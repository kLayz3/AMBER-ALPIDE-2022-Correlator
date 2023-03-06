#pragma once

#include <bits/stdc++.h>

#include "Riostream.h"
#include "TApplication.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCut.h"
#include "TCutG.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TEllipse.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLorentzVector.h"
#include "TMarker.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMultiDimFit.h"
#include "TPaletteAxis.h"
#include "TPRegexp.h"
#include "TPrincipal.h"
#include "TProfile.h"
#include "TRandom1.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TText.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TVectorT.h"
#include "Rtypes.h"

#define LEN(x) (sizeof x / sizeof *x)
#define timeNow() std::chrono::high_resolution_clock::now()

#define KMAG  "\x1B[35m"
#define KRED  "\x1B[31m"
#define KYLW  "\x1B[33m"
#define KBLUE "\x1B[34m"
#define KNRM  "\x1B[0m"
#define KGRN  "\x1B[32m"
#define KCYN  "\x1B[36m"

inline extern const std::string _help = 
"Usage:  ./main --amber=FILE_AMBER.root --alpide=FILE_ALPIDE.root\n\
\n\
Both ROOT files must be from the same run and AMBER file\n\
must contain 'usual' TRLOTimeData tree containing runNumber, spillNumber, eventNumber, eventTime, ... branches.\n\
ALPIDE file must be raw, without any clustering applied. \n\
\n\
";

#define ALPIDE_NUM 6
