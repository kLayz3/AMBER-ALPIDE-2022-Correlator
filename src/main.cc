#include "AuxFunctions.h"
#include "CMDLineParser.h"
#include "AlpideTraverser.h"
#include "AMBERTraverser.h"
#include "libs.hh"

using namespace std;
extern const std::string _help;

int main(int argc, char* argv[]) {
	if(IsCmdArg("help", argc, argv)) {cout << _help << endl; return 0;}

	string AlpideFileName;
	string AMBERFileName;

	if(!ParseCmdLine("alpide", AlpideFileName, argc,argv)) {
			cerr << "No Alpide file specified!\n";
		cerr << _help << endl; return -1;
	}
	if(!ParseCmdLine("amber", AMBERFileName, argc,argv)) {
		cerr << "No AMBER file specified!\n";
		cerr << _help << endl; return -1;
	}
	
	TFile* f_Alpide = new TFile(AlpideFileName.c_str(), "READ");
	TFile* f_AMBER  = new TFile(AMBERFileName.c_str(), "READ");
	if(!f_Alpide || f_Alpide->IsZombie()) {cerr << "Can' open Alpide rootfile.\n" << endl << _help; return -1;} 
	if(!f_AMBER || f_AMBER->IsZombie()) {cerr << "Can' open AMBER rootfile.\n" << endl << _help; return -1;} 
	
	TTree* t_Alpide = dynamic_cast<TTree*>(f_Alpide->Get("h101"));
	TTree* t_AMBER  = dynamic_cast<TTree*>(f_AMBER->Get("TRLOTimeData"));
	if(!t_Alpide || t_Alpide->IsZombie()) {cerr << "Can' open Alpide TTree.\n" << endl << _help; return -1;} 
	if(!t_AMBER || t_AMBER->IsZombie()) {cerr << "Can' open AMBER TTree.\n" << endl << _help; return -1;} 
	
	AlpideTraverser a(t_Alpide);
	AMBERTraverser A(t_AMBER);

	printf("Starting the traversal!\n");
	
	a.Go();
	a.DumpContents();

	A.Go();
	A.DumpContents();

	map<int,uint64_t> alpideMap = a.finalCookieTimestamp;
	map<int,double>   amberMap  = A.finalCookieTime;

	assert(alpideMap.size() == amberMap.size());

	printf("%s\nFound the corresponding correlated events: %s\n", KGRN, KNRM);
	printf("\n%sRun Number = %s%d%s\n", KCYN, KMAG, A.fRun, KNRM);
	for(const auto [spill, ts] : alpideMap) {
		printf("Spill %s%d%s : ALPIDE: %s0x%lx%s <==> AMBER: %s%.8f%s\n", 
				KCYN, spill, KNRM,
				KYLW, ts, KNRM,
				KBLUE, amberMap[spill], KNRM);	
	}

	f_Alpide->Close();
	f_AMBER->Close();

	cout << endl; return 0;
}
