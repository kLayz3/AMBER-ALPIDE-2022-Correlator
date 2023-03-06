#define MAX_ENTRIES 200000

TFile *fIn;
TTree *tIn;
TFile* fOut;
TTree* tOut;

/* Read containers */
int32_t chunk, outChunk;
int32_t run, outRun;
int32_t spill, outSpill;
int32_t eventNumber, outEventNumber;
int32_t eventType, outEventType;
int32_t triggerMask, outTriggerMask;

double eventTime, outEventTime;

/* Output help container */
typedef std::tuple<int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, double> Entry;
typedef std::vector<Entry> SpillVector;

std::set<int> uniqueSpills;
SpillVector v;

void SortSpillVector();
void WriteSpillVector();

void sortByEventNumber() {
	string fName = "trlo.full.root";
	string fNameOut = "sorted-"+fName;

	fIn = new TFile(fName.c_str(), "r");
	if(!fIn || fIn->IsZombie()) {cerr << "File zombied.\n"; return -1;}
	tIn = dynamic_cast<TTree*>(fIn->Get("TRLOTimeData"));
	if(!tIn || tIn->IsZombie()) {cerr << "Tree zombied.\n"; return -1;}
	
	tIn->SetBranchAddress("chunkNumber", &chunk);
	tIn->SetBranchAddress("runNumber", &run);
	tIn->SetBranchAddress("spillNumber", &spill);
	tIn->SetBranchAddress("eventNumber", &eventNumber);
	tIn->SetBranchAddress("eventType", &eventType);
	tIn->SetBranchAddress("triggerMask", &triggerMask);
	tIn->SetBranchAddress("eventTime", &eventTime);
	
	fOut = new TFile(fNameOut.c_str(), "RECREATE");
	tOut = new TTree("TRLOTimeData","TRLOTimeData");

	tOut->Branch("chunkNumber", &outChunk);
	tOut->Branch("runNumber", &outRun);
	tOut->Branch("spillNumber", &outSpill);
	tOut->Branch("eventNumber", &outEventNumber);
	tOut->Branch("eventType", &outEventType);
	tOut->Branch("triggerMask", &outTriggerMask);
	tOut->Branch("eventTime", &outEventTime);
	
	/* First identify the spills */
	uniqueSpills.clear();
	for(uint64_t i=0; i<tIn->GetEntries(); ++i) {
		tIn->GetEntry(i);
		if(uniqueSpills.find(spill) == uniqueSpills.end()) 
			uniqueSpills.insert(spill);
	}

	v.reserve(MAX_ENTRIES);

	/* Runtime complexity is O(spill x entry number), however multiple passes through 
	 * all events make it more memory efficient. It can be done in O(entry number) runtime 
	 * but then we cache the whole ROOT file. Possibly RAM overflowing .. */

	for(auto referentSpill : uniqueSpills) {
		v.clear();

		for(uint64_t i=0; i<tIn->GetEntries(); ++i) {
			tIn->GetEntry(i);
			if(spill != referentSpill) continue; 

			v.emplace_back(std::make_tuple(
						chunk,
						run,
						spill,
						eventNumber,
						eventType,
						triggerMask,
						eventTime
						));
		}
		SortSpillVector();
		WriteSpillVector();
	}

	v.clear();
	v.shrink_to_fit();
	fOut->Write();
	fOut->Close();
	fIn->Close();

	printf("\nWriting completed, exiting...\n");
}

void SortSpillVector() {
	if(v.size() == 0) return;
	std::sort(v.begin(), v.end(), [](Entry& l, Entry& r) -> bool 
			{
				return (std::get<3>(l) < std::get<3>(r)) ? 1 : 0;
			});
}

void WriteSpillVector() {	
	if(v.size() == 0) return;
	for(Entry& entry : v) {
		std::tie(outChunk, 
				 outRun,
				 outSpill,
				 outEventNumber,
				 outEventType,
				 outTriggerMask,
				 outEventTime
				 ) = entry;
		tOut->Fill();
	}
}

