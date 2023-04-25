#include "/home/frsuser/alpide-DAQ/analysis/includes/AuxFunctions.h"
#define ALPIDE_NUM 6

TFile *fIn;
TTree *tIn;
TFile* fOut;
TTree* tOut;

/* Read containers */
uint64_t t[ALPIDE_NUM+1];
int IsInSpill(double);

void sortAlpideFile() {
	string fName = "alpide-run-298967.root";
	string fNameOut = "sorted-"+fName;

	fIn = new TFile(fName.c_str(), "r");
	if(!fIn || fIn->IsZombie()) {cerr << "File zombied.\n"; return -1;}
	tIn = dynamic_cast<TTree*>(fIn->Get("h101"));
	if(!tIn || tIn->IsZombie()) {cerr << "Tree zombied.\n"; return -1;}

	t[0] = 0;
	for(int i=1; i<=ALPIDE_NUM; ++i) {
		tIn->SetBranchAddress(TString::Format("ALPIDE%dT_LO", i), (uint32_t*)&t[i]);
		tIn->SetBranchAddress(TString::Format("ALPIDE%dT_HI", i), (uint32_t*)((char*)&t[i] + 4));
	}

	uint64_t lastTs{0};
	uint64_t currTs{0};
	uint64_t cutoff = 1000;

	uint64_t eventsInSpill = 0;
	uint64_t eventsOutSpill = 0;

	uint64_t fragmentedEventsOutSpill = 0;
	uint64_t fragmentedEventsInSpill = 0;
	uint64_t ultraBadEvents = 0;

	fOut = new TFile(fNameOut.c_str(), "RECREATE");
	tOut = new TTree("h101", "h101");
	tOut->Branch("T_UNIQUE", &currTs, "T_UNIQUE/l");

	uint64_t nentries = tIn->GetEntries();
	tIn->GetEntry(0);
	uint64_t offset = *std::max_element(t+1, t+1+ALPIDE_NUM);

	for(auto i=0; i<nentries; ++i) {
		if(i%100 == 0) PrintProgress((float)i / nentries);
		tIn->GetEntry(i);
		currTs = *std::max_element(t+1, t+1+ALPIDE_NUM);
		if(currTs == 0) ++ultraBadEvents;

		double x = (currTs-offset)/1e9;
		int isIn = IsInSpill(x);
		if(isIn == -1) eventsOutSpill++;
		else eventsInSpill++;

		if(currTs - lastTs > cutoff) {
			tOut->Fill();
		}
		else {
			if(isIn == -1) fragmentedEventsOutSpill++;
			else fragmentedEventsInSpill++;
		}

		lastTs = currTs;
	}

	tOut->Write();
	fIn->Close();
	fOut->Close();

	cout << "\nDone...\n";
	printf("Found a total of %lu events in-spill.\n", eventsInSpill);
	printf("Found a total of %lu events out-spill.\n", eventsOutSpill);
	printf("Found a total of %lu fragmented entries (in-spill).\n", fragmentedEventsInSpill); 
	printf("Found a total of %lu fragmented entries (out-spill).\n", fragmentedEventsOutSpill); 
	printf("Found a total of %lu ultra-bad entries.\n", ultraBadEvents);

	delete fIn;
	delete fOut;
}

int IsInSpill(double x) {
	if(x >  5.6 && x < 9.75) return 1;
	if(x > 27.1 && x < 31.2) return 2;
	if(x > 41.5 && x < 45.7) return 3;
	if(x > 63.1 && x < 67.3) return 4;
	if(x > 84.68 && x < 88.9) return 5;
	if(x > 99.1 && x < 103.3) return 6;
	if(x > 120.7 && x < 124.9) return 7;
	if(x > 142.4 && x < 146.5) return 8;
	if(x > 156.8 && x < 160.9) return 9;

	return -1;
}
