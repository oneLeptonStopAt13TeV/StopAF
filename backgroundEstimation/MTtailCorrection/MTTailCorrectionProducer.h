//
// ########################################################################
// #  Thanks to Wouter Verkerke, awesome stat teacher, for original code  #
// #                                                                      #
// # 'ORGANIZATION AND SIMULTANEOUS FITS' RooFit tutorial macro #501      #
// # Using simultaneous p.d.f.s to describe simultaneous fits to multiple #
// # datasets                                                             #
// #                                                                      #
// # 07/2008                                                              #
// ########################################################################

// ROOT headers

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TRandom.h"
#include "TMath.h"

using namespace RooFit;

// STL headers

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

// Sonic screwdriver headers

#include "interface/Table.h"
#include "interface/Figure.h"

using namespace theDoctor;

// ####################################
// # MT tail correction configuration #
// ####################################

#define INPUT_FOLDER  "./plots/"
#define OUTPUT_FOLDER "./results/"

#define OBSERVABLE_FOR_FIT "Mlb_leadb"

#define PROCESS_NAME_TT_1L "TT_1l"
#define PROCESS_NAME_TT_2L "TT_2l"
#define PROCESS_NAME_WJETS "WJets"
#define PROCESS_NAME_RARE "OTHER"

//#define CHANNEL_NAME singleLepton
#define CHANNEL_NAME "muon"


TRandom* randomnessGenerator;

// Uncertainty of the template fit method itself,
// coming from test with MCstat, JES, algorithms, ...

#define TEMPLATE_FIT_METHOD_UNCERTAINTY 0.2



// ########################################
// # Regions for cut-and-count correction #
// ########################################

vector<string> listIndividualCuts =
{
    "MT_80",
    "MT_90",
    "MT_100",
    "MT_110"
    //"MT_120"
};

vector<string> listIndividualCuts_MTtail =
{
    "CR0b_presel_MTtail_80",
    "CR0b_presel_MTtail_90",
    "CR0b_presel_MTtail_100",
    "CR0b_presel_MTtail_110"
    //"CR0b_presel_MTtail_120"
};

vector<string>  listIndividualCuts_MTpeak =
{
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak",
    "CR0b_presel_MTpeak"
    //"CR0b_presel_MTpeak"
};

vector<string> listCutAndCounts =
{
    //"preselection"
};

map<string,vector<string> > listCutAndCounts_cuts;

void initCutAndCountCuts()
{
    //listCutAndCounts_cuts["preselection"]                     = { "MT_100"                     };
}
