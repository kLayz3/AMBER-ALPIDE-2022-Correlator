#include "AMBERTraverser.h"
#include "AuxFunctions.h"

using namespace std;
using namespace DiffAMBER;

AMBERTraverser::AMBERTraverser(TTree* tree) : fTree(tree), fEntry(0) {
	if(!fTree | fTree->IsZombie()) throw runtime_error("AMBERTraverser::AMBERTraverser() - Failed to construct tree.\n");

	fMaxEntries = tree->GetEntries();
		
	fTree->SetBranchAddress("runNumber", &fRun);
	fTree->SetBranchAddress("spillNumber", &fSpill);
	fTree->SetBranchAddress("eventNumber", &fEventNumber);
	fTree->SetBranchAddress("eventType", &fEventType);
	fTree->SetBranchAddress("eventTime", &fEventTime);
	
	prevSpill = 1;
	currSpill = 1;

	currEvent = 0;
	prevEvent = 0;

	tcurr = -1;
	tprev = -1;
	tdiff = -1;
}

/* Returns true if a spill jump happens, else returns false */
template<class T> bool AMBERTraverser::GetEntry(T entry) {
	fTree->GetEntry(entry);
	fEntry = (uint64_t)entry;

	tprev = tcurr;
	tcurr  = fEventTime;
	tdiff  = tcurr - tprev;

	prevEvent = currEvent;
	currEvent = fEventNumber;

	prevSpill = currSpill;
	currSpill = fSpill;
	return (currSpill != prevSpill);
}

TriggerType AMBERTraverser::GetTriggerType() {
	if(tcurr < 0) return TriggerType::kUNKNOWN;

	if(tdiff > kOFFSPILL) return TriggerType::kOFFSPILL;
	if(abs(tdiff - kCOOKIE) < kCOOKIE_WIDTH) return TriggerType::kCOOKIE;
	if(tdiff < kSPILL && tdiff > 0) return TriggerType::kSPILL;
	
	return TriggerType::kUNKNOWN;
}

/* A 'cookie' event is when two triggers come within [kCOOKIE +- kCOOKIE_DIST] window
 * Two 'cookie' events are 'good' if they are [kCOOKIE_DIST +- kCOOKIE_DIST_WIDTH] apart.
 * This 'periodic clock' is identifiable in both datasets! */

void AMBERTraverser::IdentifyCookies() {
	std::vector<double>   cookieTime;
	std::vector<int>      cookieEvent;
	std::vector<uint64_t> cookieEntry;

	cookieTime.reserve(1000);
	cookieEvent.reserve(1000);
	cookieEntry.reserve(1000);

	double   prevCookieTime = -1;
	double   currCookieTime = -1;
	int      prevCookieEvent = -1;
	int      currCookieEvent = -1;
	uint64_t prevCookieEntry = 0;
	uint64_t currCookieEntry = 0;

	bool foundFirstGoodCookie = false;

	GetEntry(0); GetEntry(1); // fetch first two manually so fields are proper
	printf("\n%s[AMBER]%s Starting initial spill: %s%d%s\n", KBLUE, KGRN, KCYN, currSpill, KNRM);
	for(uint64_t i=2; i<fMaxEntries; ++i) {
		PrintProgress(i/(float)fMaxEntries);

		bool b = GetEntry(i);
		if(b) {
			printf("\n%s [AMBER]%s Done with spill: %s%d%s\n", KBLUE, KGRN, KCYN, prevSpill, KNRM);
			
			firstCookieTime[prevSpill]  = cookieTime.front();
			finalCookieTime[prevSpill]  = cookieTime.back();
			firstCookieEvent[prevSpill] = cookieEvent.front();
			finalCookieEvent[prevSpill] = cookieEvent.back();
			firstCookieEntry[prevSpill] = cookieEntry.front();
			finalCookieEntry[prevSpill] = cookieEntry.back();

			numberOfCookies[prevSpill] = cookieTime.size();

			cookieTime.clear(); 
			cookieEvent.clear();
			cookieEntry.clear();
			foundFirstGoodCookie = false;
		}
		
		TriggerType type = GetTriggerType();
		
		if(type == TriggerType::kCOOKIE) {
			prevCookieTime = currCookieTime;
			currCookieTime = (tcurr + tprev)/2;
			
			prevCookieEvent = currCookieEvent;
			currCookieEvent = currEvent;

			prevCookieEntry = currCookieEntry;
			currCookieEntry = i;
		
			/* We're only interested in cookies who come kCOOKIE_DIST apart */
			if(abs(currCookieTime-prevCookieTime - kCOOKIE_DIST) < kCOOKIE_DIST_WIDTH) {
				if(!foundFirstGoodCookie) {
					cookieTime.push_back(prevCookieTime);
					cookieEvent.push_back(prevCookieEvent);
					cookieEntry.push_back(prevCookieEntry);
					foundFirstGoodCookie = true;
				}
				
				cookieTime.push_back(currCookieTime);
				cookieEvent.push_back(currCookieEvent);
				cookieEntry.push_back(currCookieEntry);
			}
		}
	}
	
	/* Once loop is done, update the maps with final spill */
	printf("\n%s [AMBER]%s Done with spill: %s%d%s\n", KBLUE, KGRN, KCYN, prevSpill, KNRM);
	firstCookieTime[prevSpill]  = cookieTime.front();
	finalCookieTime[prevSpill]  = cookieTime.back();
	firstCookieEvent[prevSpill] = cookieEvent.front();
	finalCookieEvent[prevSpill] = cookieEvent.back();
	firstCookieEntry[prevSpill] = cookieEntry.front();
	finalCookieEntry[prevSpill] = cookieEntry.back();

	numberOfCookies[prevSpill] = cookieTime.size();
}

/* This method will save all the triggers that come after the main part of the spill ends.
 * They should always be kCOOKIE_DIST apart and exactly 1 per such a distance */
void AMBERTraverser::IdentifyCalibrationTrigs() {
	/* Go from the startEntry until the spill ends and try identify the type==8 trigs.
	 * If they're type==8 and kCOOKIE_DIST apart, save them into vector */

	printf("\n%s[AMBER]%s Identifying 'calibration triggers' after each of the spills ... %s\n", KBLUE, KGRN, KNRM);
	for(const auto [spill, startEntry] : finalCookieEntry) {
		PrintProgress(spill / (float)finalCookieEntry.size());
		vector<double> trigTime;

		for(uint64_t i=startEntry; i<fMaxEntries; ++i) {
			bool b = GetEntry(i);
			if(!b) break;

			if(fEventType == 8 && abs(tdiff - kCOOKIE_DIST) < kCOOKIE_DIST_WIDTH) {
				if(trigTime.size() == 0)
					trigTime.push_back(tprev);
				trigTime.push_back(tcurr);
			}
		}
		calibrationTrigs[spill] = std::move(trigTime);
	}
}

void AMBERTraverser::WriteToFilePretty(const char* fileName) {
	std::fstream f;
	f.open(fileName, ios_base::out);
	if(!f.is_open()) throw std::invalid_argument("AMBERTraverser::WriteToFilePretty -Wrong file name?");
	f << "Run = " << fRun << endl;
	f << " ----------------------------- " << endl;
	f << std::setprecision(9) << std::fixed;
	for(const auto [spill, cookieNum] : numberOfCookies) {
		double firstTime = firstCookieTime[spill];
		double finalTime = finalCookieTime[spill];
		int firstEvent = firstCookieEvent[spill];
		int finalEvent = finalCookieEvent[spill];
		uint64_t firstEntry  = firstCookieEntry[spill];
		uint64_t finalEntry  = finalCookieEntry[spill];

		f << "Spill                    : " << spill << endl;
		f << "First cookie entry number: " << firstEntry << endl;
		f << "First cookie event number: " << firstEvent << endl;
		f << "First cookie event time  : " << firstTime << endl;

		f << "Final cookie entry number: " << finalEntry << endl;
		f << "Final cookie event number: " << finalEvent << endl;
		f << "Final cookie event time  : " << finalTime << endl;

		f << "Number of 'good cookies' : " << numberOfCookies[spill] << endl;
		f << " ----------------------------- " << endl;
	}
}

void AMBERTraverser::WriteToFile(const char* fileName) {}

void AMBERTraverser::DumpContents() {
	cout << "\nRun = " << fRun << endl;
	cout << "-----------------------------"  << endl;
	cout << std::setprecision(9) << std::fixed;
	for(const auto [spill, cookieNum] : numberOfCookies) {
		double firstTime = firstCookieTime[spill];
		double finalTime = finalCookieTime[spill];
		int firstEvent = firstCookieEvent[spill];
		int finalEvent = finalCookieEvent[spill];
		uint64_t firstEntry  = firstCookieEntry[spill];
		uint64_t finalEntry  = finalCookieEntry[spill];

		cout << "Spill                    : " << spill << endl;
		cout << "First cookie entry number: " << firstEntry << endl;
		cout << "First cookie event number: " << firstEvent << endl;
		cout << "First cookie event time  : " << firstTime << endl;
	
		cout << "Final cookie entry number: " << finalEntry << endl;
		cout << "Final cookie event number: " << finalEvent << endl;
		cout << "Final cookie event time  : " << finalTime << endl;

		cout << "Number of 'good cookies' : " << numberOfCookies[spill] << endl;
		cout << "-----------------------------" << endl;
	}
}

void AMBERTraverser::Go() {
	/* IdentifySpills(); */
	IdentifyCookies();
	IdentifyCalibrationTrigs();
}
