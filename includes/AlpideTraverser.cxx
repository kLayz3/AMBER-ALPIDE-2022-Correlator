#include "AlpideTraverser.h"
#include "AuxFunctions.h"

using namespace std;
using namespace DiffAlpide;

AlpideTraverser::AlpideTraverser(TTree* h101) : fTree(h101), fEvent(0), fSpill(0), fInSpill(false) {
	if(!fTree | fTree->IsZombie()) throw runtime_error("AlpideTraverser::AlpideTraverser() - Failed to construct tree");

	fMaxEvents = h101->GetEntries();
	for(int x=1; x<=ALPIDE_NUM; ++x) {
		fTree->SetBranchAddress(TString::Format("ALPIDE%dT_LO", x), (uint32_t*)&t[x]);
		fTree->SetBranchAddress(TString::Format("ALPIDE%dT_HI", x), (uint32_t*)((char*)&t[x] + 4));
	}
}

bool AlpideTraverser::IsEntryValid() {
	for(int i=1; i<=ALPIDE_NUM; ++i) {
		if(t[i] == 0) return false;
	}
	return true;
}

void AlpideTraverser::FindInitialEvent() {
	for(uint64_t i=0; i<fMaxEvents; ++i) {
		fTree->GetEntry(i);
		if(IsEntryValid()) {
			fInitialTS    = t[1];
			fInitialEvent = i;
			tcurr = 0;
			tprev = 0;
			break;
		}
	}
}

uint64_t AlpideTraverser::SecondsToTimestamp(double t) {return fInitialTS + (uint64_t)(t * 1e9);}
double	 AlpideTraverser::TimestampToSeconds(uint64_t ts) {return (ts - fInitialTS)/1e9;}

template<class T> bool AlpideTraverser::GetEntry(T entry) {
	if(entry >= (T)fMaxEvents) {cerr<<"Entry > MaxEvents\n"; return false;}	
	
	fTree->GetEntry(entry);
	if(!IsEntryValid()) return false;

	tprev = tcurr;
	tcurr  = (t[1] - fInitialTS)/1e9;
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
	for(uint64_t i=0; i<fMaxEvents; ++i) {
		PrintProgress(i/(float)(fMaxEvents));
		bool b = GetEntry(i);
		if(!b) continue;
		
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
	uint64_t endEvent   = endOfSpill.at(spill);
	uint64_t startEvent = startOfSpill.at(spill);
	// If entry not found, then std::out_of_range exception thrown
	
	printf("\n%s[ALPIDE]%s Starting spill: %s%d%s\n", KRED, KGRN, KCYN, spill, KNRM);

	std::vector<uint64_t> cookieTime;
	std::vector<uint64_t> cookieEvent;

	for(uint64_t i=startEvent; i<=endEvent; ++i) {
		PrintProgress((i - startEvent)/(float)(endEvent - startEvent));
		bool b = GetEntry(i);
		if(!b) continue;
		
		TriggerType type = GetTriggerType();
		
		if(type == TriggerType::kCOOKIE) {
			cookieTime.push_back(t[1]);
			cookieEvent.push_back(i);
		}
	}
	
	firstCookieTimestamp[spill] = cookieTime.front();
	finalCookieTimestamp[spill] = cookieTime.back();
	firstCookieEntry[spill] = cookieEvent.front();
	finalCookieEntry[spill] = cookieEvent.back();
	
	numberOfCookies[spill] = cookieTime.size();
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
	FindInitialEvent();
	IdentifySpills();
	for(const auto [spill, spillStart] : startOfSpill)
		IdentifyCookies(spill);
}
