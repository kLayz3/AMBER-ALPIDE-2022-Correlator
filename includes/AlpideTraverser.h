#ifndef ALPIDE_TRAVERSER_H
#define ALPIDE_TRAVERSER_H

#include "libs.hh"

namespace DiffAlpide {
	const double kOFFSPILL     = 0.002;
	const double kCOOKIE       = 0.00034;
	const double kCOOKIE_WIDTH = 0.00004;
	const double kSPILL        = 0.00025;

	const double kCOOKIE_DIST		= 0.033712;
	const double kCOOKIE_DIST_WIDTH = 0.00005;
	const double kUNKNOWN	   = -1;

	enum class TriggerType {
		kOFFSPILL,
		kSPILL,
		kCOOKIE,
		kUNKNOWN
	};
};

class AlpideTraverser {
public:
	TTree* fTree;           // Should point to h101
	uint64_t fEvent;        // Will loop up to max events;
	uint64_t fSpill;        // On which spill the traverser is.
	uint64_t fMaxEvents;    // Max Events in the ROOT file
	
	uint64_t fInitialTS;    // Time stamp of the intial event
	uint64_t fInitialEvent; // Index of initial event that holds the TS data.

	uint64_t t[ALPIDE_NUM + 1] = {0};
	
	bool fInSpill;
	double tprev;
	double tcurr;
	double tdiff;

	/* Map is used because in AMBER data spills start from index 1, not from index 0.
	 * So to avoid confusion when using potential ::size() operator, I just save everything into a map */
	std::map<int, uint64_t> firstCookieTimestamp;
	std::map<int, uint64_t> finalCookieTimestamp;
	std::map<int, uint64_t> firstCookieEntry;
	std::map<int, uint64_t> finalCookieEntry;
	std::map<int, int>    numberOfCookies;

	std::map<int, uint64_t> startOfSpill;
	std::map<int, uint64_t> endOfSpill;

	AlpideTraverser() = default;
	AlpideTraverser(TTree* h101);

	bool IsEntryValid();
	void FindInitialEvent();

	uint64_t SecondsToTimestamp(double);
	double TimestampToSeconds(uint64_t);

	template<class T> bool GetEntry(T entry);
	int NumberOfSpills() {return endOfSpill.size();}

	DiffAlpide::TriggerType GetTriggerType();

	void IdentifySpills();
	void IdentifyCookies(int spill);

	void WriteToFilePretty(const char* fileName);
	void WriteToFile(const char* fileName);
	void DumpContents();

	void Go();
};

#endif
