#include "AlpideTraverser.h"
#include "AuxFunctions.h"

using namespace std;
using namespace DiffAlpide;

AlpideTraverser::AlpideTraverser(TTree* h101) : fTree(h101), fEvent(0), fSpill(0), fInSpill(false) {
	if(!fTree | fTree->IsZombie()) throw runtime_error("AlpideTraverser::AlpideTraverser() - Failed to construct tree");

	fMaxEntries = h101->GetEntries();
	fTree->SetBranchAddress("T_UNIQUE", &ts);
	fTree->GetEntry(0);
	fInitialTS = ts;
	tcurr = 0;
	tprev = 0;
}

uint64_t AlpideTraverser::SecondsToTimestamp(double t) {return fInitialTS + (uint64_t)(t * 1e9);}
double	 AlpideTraverser::TimestampToSeconds(uint64_t ts) {return (ts - fInitialTS)/1e9;}

template<class T> bool AlpideTraverser::GetEntry(T entry) {
	if(entry >= (T)fMaxEntries) {cerr<<"Entry > MaxEvents\n"; return false;}	
	fTree->GetEntry(entry);
	tprev = tcurr;
	tcurr  = (ts - fInitialTS)/1e9;
	tdiff  = tcurr - tprev;
	fEvent = (uint64_t)entry;
	return true;
}

TriggerType AlpideTraverser::GetTriggerType() {
	if(tdiff > kOFFSPILL) return TriggerType::kOFFSPILL;
	if(abs(tdiff - kCOOKIE) < kCOOKIE_WIDTH) return TriggerType::kCOOKIE;
	if(tdiff < kSPILL && tdiff > 0) return TriggerType::kSPILL;
	
	return TriggerType::kUNKNOWN;
}

/* This method fills the maps with info about when spill starts/ends
 * And also on which event this happens. Always remember tcurr is relative to 
 * tInitialTS (divided by 1e9)
 */

void AlpideTraverser::IdentifySpills() {
	printf("\n%s[ALPIDE]%s Identifying spills ... %s\n", KRED, KGRN, KNRM);
	
	double prevCookie = -1;
	double currCookie = -1;
	uint64_t prevCookieEntry = 0;
	uint64_t currCookieEntry = 0;
	
	fInSpill = false;
	for(uint64_t i=0; i<fMaxEntries; ++i) {
		PrintProgress(i/(float)(fMaxEntries));
		GetEntry(i);
		
		TriggerType type = GetTriggerType();

		if(fInSpill && type == TriggerType::kOFFSPILL) {
			/* Means spill ends. 
			 * During main region of the spill there cannot be
			 * triggers apart more than kOFFSPILL. */
			fInSpill = false;
			endOfSpill[fSpill] = i;
		}

		/* To figure out the 'spill-start' : 
		 * If we encounter a cookie interval which is correctly separated from previous
		 * cookie interval, then we label this (previous cookie trigger) as 'spill start'.
		 */

		if(!fInSpill && type == TriggerType::kCOOKIE) {
			prevCookie = currCookie;
			currCookie = (tcurr + tprev)/2;
			prevCookieEntry = currCookieEntry;
			currCookieEntry = i;

			if(abs(currCookie-prevCookie - kCOOKIE_DIST) < kCOOKIE_DIST_WIDTH) {
				fInSpill = true;
				++fSpill;
				// Assign a 'spill start' couple entries before the initial 'good' cookie event
				startOfSpill[fSpill] = std::max(prevCookieEntry-3, 0ul);
			}
		}
	}

	printf("\n%s[ALPIDE]%s Done with Identify spills. Found %s%d%s spills.%s\n", KRED, KGRN, KCYN, (int)startOfSpill.size(), KGRN, KNRM);
}

/* With the info from the call above, we go over a spill.
 * Since the spill 'end' is abrupt, we can clearly identify the
 * final cookie event in a spill. They *should* be periodic with
 * a period ~ kCOOKIE_DIST. 
 * 
 * A 'cookie' event is when two triggers come within [kCOOKIE +- kCOOKIE_WIDTH] window
 * Two 'cookie' events are 'good' if they are [kCOOKIE_DIST +- kCOOKIE_DIST_WIDTH] apart.
 * This 'period clock' is identifiable in both datasets! 
 *
 * From inspecting the data, it's assumed all the cookie events in the main part of the spill
 * are correctly separated and considered 'good' */

void AlpideTraverser::IdentifyCookies(int spill) {
	uint64_t endEntry   = endOfSpill.at(spill);
	uint64_t startEntry = startOfSpill.at(spill);
	// If entry not found, then std::out_of_range exception thrown
	
	printf("\n%s[ALPIDE]%s Starting spill: %s%d%s\n", KRED, KGRN, KCYN, spill, KNRM);

	std::vector<uint64_t> cookieTimestamp;
	std::vector<double> cookieTime;
	std::vector<uint64_t> cookieEntry;

	for(uint64_t i=startEntry; i<=endEntry; ++i) {
		PrintProgress((i - startEntry)/(float)(endEntry - startEntry));
		GetEntry(i);
		
		TriggerType type = GetTriggerType();
		
		if(type == TriggerType::kCOOKIE) {
			cookieTimestamp.push_back(ts);
			cookieTime.push_back(tcurr);
			cookieEntry.push_back(i);
		}
	}
	
	firstCookieTimestamp[spill] = cookieTimestamp.front();
	finalCookieTimestamp[spill] = cookieTimestamp.back();
	
	firstCookieTime[spill] = cookieTime.front();
	finalCookieTime[spill] = cookieTime.back();

	firstCookieEntry[spill] = cookieEntry.front();
	finalCookieEntry[spill] = cookieEntry.back();
	
	numberOfCookies[spill] = cookieTimestamp.size();
}

/* This method will save all the triggers that come after the main part of the spill ends.
 * They should always be kCOOKIE_DIST apart and exactly 1 per such a distance */
void AlpideTraverser::IdentifyCalibrationTrigs() {
	printf("\n%s[ALPIDE]%s - Identifying 'calibration' triggers after each of the spills ... %s\n", KRED, KGRN, KNRM);
	for(const auto [spill, startEntry] : finalCookieEntry) {
		PrintProgress(spill / (float)finalCookieEntry.size());
		vector<double> trigTime;

		uint64_t endEntry;
		try		   {endEntry = firstCookieEntry.at(spill+1);}
		catch(...) {endEntry = fMaxEntries;}

		for(uint64_t i=startEntry; i<endEntry; ++i) {
			bool b = GetEntry(i);
			if(!b) continue;

			if(abs(tdiff - kCOOKIE_DIST) < kCOOKIE_DIST_WIDTH) {
				if(trigTime.size() == 0)
					trigTime.push_back(tprev);
				trigTime.push_back(tcurr);
			}
		}
		calibrationTrigs[spill] = std::move(trigTime);
	}
}

void AlpideTraverser::WriteToFilePretty(const char* fileName) {
	std::fstream f;
	f.open(fileName, ios_base::out);
	if(!f.is_open()) throw std::invalid_argument("Wrong file name?");
	f << "Initial starting TS (64-bit): 0x" << std::hex << fInitialTS << std::dec << endl;
	f << "-----------------------------" << endl;
	for(const auto [spill, spillStart] : startOfSpill) {
		uint64_t spillEnd = endOfSpill[spill];

		GetEntry(spillStart);
		double tStart = tcurr;
		GetEntry(spillEnd);
		double tEnd   = tcurr;

		f << "Spill			    : " << spill << endl;
		f << "Spill start entry : " << spillStart << endl;
		f << "Spill start time  : " << tStart << endl;
		f << "Spill end entry   : " << spillEnd << endl;
		f << "Spill end time    : " << tEnd << endl;
		f << std::hex;
		f << "First Cookie TS          : " << firstCookieTimestamp[spill] << endl;
		f << "Final Cookie TS          : " << finalCookieTimestamp[spill] << endl;
		f << std::dec;
		f << "Number of 'good cookies' : " << numberOfCookies[spill] << endl;
		f << "-----------------------------" << endl;
	}
	f.close();
}

void AlpideTraverser:: WriteToFile(const char* fileName) {}

void AlpideTraverser::DumpContents() {
	cout << "\nInitial starting TS (64-bit): 0x" << std::hex << fInitialTS << std::dec << endl;
	cout << "-----------------------------" << endl;
	for(const auto [spill, spillStart] : startOfSpill) {
		uint64_t spillEnd = endOfSpill[spill];

		GetEntry(spillStart);
		double tStart = tcurr;
		GetEntry(spillEnd);
		double tEnd   = tcurr;

		cout << "Spill			    : " << spill << endl;
		cout << "Spill start entry : " << spillStart << endl;
		cout << "Spill start time  : " << tStart << endl;
		cout << "Spill end entry   : " << spillEnd << endl;
		cout << "Spill end time    : " << tEnd << endl;
		cout << std::hex;
		cout << "First Cookie TS          : " << firstCookieTimestamp[spill] << endl;
		cout << "Final Cookie TS          : " << finalCookieTimestamp[spill] << endl;
		cout << std::dec;
		cout << "Number of 'good cookies' : " << numberOfCookies[spill] << endl;
		cout << "-----------------------------" << endl;
	}
}

void AlpideTraverser::Go() {
	IdentifySpills();
	for(const auto [spill, spillStart] : startOfSpill)
		IdentifyCookies(spill);

	IdentifyCalibrationTrigs();
}
