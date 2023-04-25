#ifndef AMBER_TRAVERSER_H
#define AMBER_TRAVERSER_H 

#include "libs.hh"

namespace DiffAMBER {
	const double kOFFSPILL     = 0.00063;
	const double kCOOKIE       = 0.0003062;
	const double kCOOKIE_WIDTH = 0.0000080;
	const double kSPILL        = 0.0000100;
	
	const double kCOOKIE_DIST		  = 0.033712;
	const double kCOOKIE_DIST_WIDTH   = 0.000011;
	const double kUNKNOWN = -1;

	enum class TriggerType {
		kOFFSPILL,
		kSPILL,
		kCOOKIE,
		kUNKNOWN
	};
};

class AMBERTraverser {
public:
	TTree* fTree;	   // Should point to TRLOTimeData
	int fRun;		   // Branch address for runNumber
	int fSpill;	       // Branch address for spillNumber
	int fEventNumber;  // Branch address for eventNumber
	int fEventType;	   // Branch address for eventType
	double fEventTime; // Branch address for eventTime
	
	uint64_t fEntry;      // Will loop up to max entries;
	uint64_t fMaxEntries; // Entry number in the ROOT file

	int prevSpill; // spillNumber of previous entry
	int currSpill; // spillNumber of current entry

	int prevEvent; // eventNumber of previous entry
	int currEvent; // eventNumber of current entry 

	double tprev; // eventTime of previous entry
	double tcurr; // eventTime of current entry
	double tdiff; // difference between them

	/* std::map is used because in AMBER data spills start from index 1, not from index 0.
	 * So to avoid confusion when using potential ::size() operator, I just save everything into a map  
	 * The key of the map is always the spillNumber */

	std::map<int, double>   firstCookieTime;
	std::map<int, double>   finalCookieTime;
	std::map<int, int>      firstCookieEvent;
	std::map<int, int>      finalCookieEvent;
	std::map<int, uint64_t> firstCookieEntry;
	std::map<int, uint64_t> finalCookieEntry;
	std::map<int, int> numberOfCookies;
	std::map<int, std::vector<double>> calibrationTrigs; // trigs type==8 after main spill part ends

	AMBERTraverser() = default;
	AMBERTraverser(TTree* tree);
	
	template<class T> bool GetEntry(T entry);

	DiffAMBER::TriggerType GetTriggerType();

	void IdentifyCookies();
	void IdentifyCalibrationTrigs();

	void WriteToFilePretty(const char* fileName);
	void WriteToFile(const char* fileName);
	void DumpContents();

	void Go();
};

#endif
