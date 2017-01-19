#include <vector>
#include <iostream>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"

#define USE_TTZ
//#define USE_TTZNLO
//#define USE_WZ
//#define USE_ZZ

#define USE_VAR_BASELINE
#define USE_LEP1
#define USE_LEP2
#define USE_JETS
#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_GLOBAL_VAR

#include "../../common/TFFactory.h"
#include "../../Selection/moriondPlots.h"

using namespace std;


// ----------------------------------------------
// Should be called only here because many
// struct and fuctions have to be declare first
// ----------------------------------------------
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"

uint32_t counter = 0;
string empty = "";
string storedDataset;
TH2D *h2 = NULL;
TAxis *xaxis = NULL;
TAxis *yaxis = NULL;
bool checkNegativeYields = false;
uint32_t nthentry = 0;
string outputName = "";


float getWeight(string currentProcessType, float lumi);
float reweightTop(float topPt, float atopPt);
map< pair<uint32_t,uint32_t>, string > scanMap;

bool lepChannel() 
{ 
    return true; 
}
    
ofstream statnames("statNamesPlots.txt");
//Add this as a global variable 

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    //babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11/NotSkimmed";
    babyTuplePath = "/opt/sbg/data/data1/cms/echabert/Stop1lSharedBabies/isuarez_v11";
    totalNumberOfWorkers = 1;

    vector<float> METBins1 = {250,350,450,600};
    vector<float> METBins2 = {250,350,450,550,650};
    vector<float> METBins3 = {250,350,450,550};
    vector<float> topnessB = {-20,0,10,20};

    AddVariable("MET", "E_{T}^{miss}",  "GeV", 4 ,250, 650,  &(myEvent.pfmet), "noUnderflowInFirstBin,logY");
    AddVariable("METAB", "E_{T}^{miss}",  "GeV", (int) (METBins1.size()-1), METBins1.data(),  &(myEvent.pfmet), "noUnderflowInFirstBin,logY");
    AddVariable("METCD", "E_{T}^{miss}",  "GeV", (int) (METBins2.size()-1), METBins2.data(),  &(myEvent.pfmet),  "noUnderflowInFirstBin,logY");
    AddVariable("MET3EFGHI", "E_{T}^{miss}",  "GeV", (int) (METBins3.size()-1), METBins3.data(),  &(myEvent.pfmet), "noUnderflowInFirstBin,logY");
    AddVariable("Mlb", "M_{lb}", "GeV", 2 , 0, 350,  &(myEvent.Mlb), "noUnderflowInFirstBin,logY");
    AddVariable("topnessMod", "t_{mod}", "GeV", (int)(topnessB.size()-1),topnessB.data(),  &(myEvent.topnessMod), "logY");
    AddVariable("Njets", "N_{jets}", "",5,2,7,  &(myEvent.ngoodjets));

    #ifdef USE_TTZ
    AddProcessClass("ttZ", "ttZ", "background", kWhite);
    	AddDataset("ttZJets_13TeV_madgraphMLM","ttZ",0,0);
/*    	AddDataset("ttZJets_13TeV_madgraphMLM_1","ttZ",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM_2","ttZ",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM_3","ttZ",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM_4","ttZ",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM_5","ttZ",0,0);
    	AddDataset("ttZJets_13TeV_madgraphMLM_6","ttZ",0,0);*/
        outputName = "METplotsttZ";
    #endif

    #ifdef USE_TTZNLO
    AddProcessClass("ttZ", "ttZ", "background", kWhite);
        AddDataset("TTZToLLNuNu_M-10_amcnlo_pythia8_25ns","ttZ",0,0);
    //outputName = "yieldZnunuMorTTZNLO";
    outputName = "METplotsttZNLO";
    #endif
    
	
    #ifdef USE_ZZWZ
    AddProcessClass("ZZWZ", "ZZWZ", "background", kBlue);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","ZZWZ",0,0);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZWZ",0,0);
        outputName = "METplotsZZWZ";
    #endif

    #ifdef USE_ZZ
    AddProcessClass("ZZ", "ZZ", "background", kBlue);
    	AddDataset("ZZTo2Q2Nu_amcnlo_pythia8_25ns","ZZ",0,0);
        outputName = "METplotsZZ";
    #endif

    #ifdef USE_WZ
    AddProcessClass("WZ", "WZ", "background", kBlue);
    	AddDataset("WZTo1L3Nu_amcnlo_pythia8_25ns","WZ",0,0);
        outputName = "METplotsWZ";
    #endif

	//regions

AddRegion("SR1l_AB_250lessMETlessInf","SR1l_AB_250lessMETlessInf",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfLSFdown","SR1l_AB_250lessMETlessInfLSFdown",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfLSFup","SR1l_AB_250lessMETlessInfLSFup",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfBTlightDown","SR1l_AB_250lessMETlessInfBTlightDown",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfBTlightUp","SR1l_AB_250lessMETlessInfBTlightUp",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfBTheavyDown","SR1l_AB_250lessMETlessInfBTheavyDown",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfBTheavyUp","SR1l_AB_250lessMETlessInfBTheavyUp",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfPUdown","SR1l_AB_250lessMETlessInfPUdown",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfPUup","SR1l_AB_250lessMETlessInfPUup",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfPDFdown","SR1l_AB_250lessMETlessInfPDFdown",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfPDFup","SR1l_AB_250lessMETlessInfPDFup",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfalphaSdown","SR1l_AB_250lessMETlessInfalphaSdown",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfalphaSup","SR1l_AB_250lessMETlessInfalphaSup",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfQ2down","SR1l_AB_250lessMETlessInfQ2down",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_AB_250lessMETlessInfQ2up","SR1l_AB_250lessMETlessInfQ2up",&SR1l_AB_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInf","SR1l_CD_250lessMETlessInf",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfLSFdown","SR1l_CD_250lessMETlessInfLSFdown",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfLSFup","SR1l_CD_250lessMETlessInfLSFup",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfBTlightDown","SR1l_CD_250lessMETlessInfBTlightDown",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfBTlightUp","SR1l_CD_250lessMETlessInfBTlightUp",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfBTheavyDown","SR1l_CD_250lessMETlessInfBTheavyDown",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfBTheavyUp","SR1l_CD_250lessMETlessInfBTheavyUp",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfPUdown","SR1l_CD_250lessMETlessInfPUdown",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfPUup","SR1l_CD_250lessMETlessInfPUup",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfPDFdown","SR1l_CD_250lessMETlessInfPDFdown",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfPDFup","SR1l_CD_250lessMETlessInfPDFup",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfalphaSdown","SR1l_CD_250lessMETlessInfalphaSdown",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfalphaSup","SR1l_CD_250lessMETlessInfalphaSup",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfQ2down","SR1l_CD_250lessMETlessInfQ2down",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_CD_250lessMETlessInfQ2up","SR1l_CD_250lessMETlessInfQ2up",&SR1l_CD_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInf","SR1l_EFGH_250lessMETlessInf",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfLSFdown","SR1l_EFGH_250lessMETlessInfLSFdown",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfLSFup","SR1l_EFGH_250lessMETlessInfLSFup",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfBTlightDown","SR1l_EFGH_250lessMETlessInfBTlightDown",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfBTlightUp","SR1l_EFGH_250lessMETlessInfBTlightUp",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfBTheavyDown","SR1l_EFGH_250lessMETlessInfBTheavyDown",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfBTheavyUp","SR1l_EFGH_250lessMETlessInfBTheavyUp",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfPUdown","SR1l_EFGH_250lessMETlessInfPUdown",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfPUup","SR1l_EFGH_250lessMETlessInfPUup",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfPDFdown","SR1l_EFGH_250lessMETlessInfPDFdown",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfPDFup","SR1l_EFGH_250lessMETlessInfPDFup",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfalphaSdown","SR1l_EFGH_250lessMETlessInfalphaSdown",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfalphaSup","SR1l_EFGH_250lessMETlessInfalphaSup",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfQ2down","SR1l_EFGH_250lessMETlessInfQ2down",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_EFGH_250lessMETlessInfQ2up","SR1l_EFGH_250lessMETlessInfQ2up",&SR1l_EFGH_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInf","SR1l_I_250lessMETlessInf",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfLSFdown","SR1l_I_250lessMETlessInfLSFdown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfLSFup","SR1l_I_250lessMETlessInfLSFup",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfBTlightDown","SR1l_I_250lessMETlessInfBTlightDown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfBTlightUp","SR1l_I_250lessMETlessInfBTlightUp",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfBTheavyDown","SR1l_I_250lessMETlessInfBTheavyDown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfBTheavyUp","SR1l_I_250lessMETlessInfBTheavyUp",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfPUdown","SR1l_I_250lessMETlessInfPUdown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfPUup","SR1l_I_250lessMETlessInfPUup",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfPDFdown","SR1l_I_250lessMETlessInfPDFdown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfPDFup","SR1l_I_250lessMETlessInfPDFup",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfalphaSdown","SR1l_I_250lessMETlessInfalphaSdown",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfalphaSup","SR1l_I_250lessMETlessInfalphaSup",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfQ2down","SR1l_I_250lessMETlessInfQ2down",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_I_250lessMETlessInfQ2up","SR1l_I_250lessMETlessInfQ2up",&SR1l_I_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInf","SR1l_CDEFGH_250lessMETlessInf",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfLSFdown","SR1l_CDEFGH_250lessMETlessInfLSFdown",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfLSFup","SR1l_CDEFGH_250lessMETlessInfLSFup",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfBTlightDown","SR1l_CDEFGH_250lessMETlessInfBTlightDown",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfBTlightUp","SR1l_CDEFGH_250lessMETlessInfBTlightUp",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfBTheavyDown","SR1l_CDEFGH_250lessMETlessInfBTheavyDown",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfBTheavyUp","SR1l_CDEFGH_250lessMETlessInfBTheavyUp",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfPUdown","SR1l_CDEFGH_250lessMETlessInfPUdown",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfPUup","SR1l_CDEFGH_250lessMETlessInfPUup",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfPDFdown","SR1l_CDEFGH_250lessMETlessInfPDFdown",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfPDFup","SR1l_CDEFGH_250lessMETlessInfPDFup",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfalphaSdown","SR1l_CDEFGH_250lessMETlessInfalphaSdown",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfalphaSup","SR1l_CDEFGH_250lessMETlessInfalphaSup",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfQ2down","SR1l_CDEFGH_250lessMETlessInfQ2down",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_CDEFGH_250lessMETlessInfQ2up","SR1l_CDEFGH_250lessMETlessInfQ2up",&SR1l_CDEFGH_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInf","SR1l_NJlowTM_250lessMETlessInf",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfLSFdown","SR1l_NJlowTM_250lessMETlessInfLSFdown",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfLSFup","SR1l_NJlowTM_250lessMETlessInfLSFup",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfBTlightDown","SR1l_NJlowTM_250lessMETlessInfBTlightDown",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfBTlightUp","SR1l_NJlowTM_250lessMETlessInfBTlightUp",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfBTheavyDown","SR1l_NJlowTM_250lessMETlessInfBTheavyDown",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfBTheavyUp","SR1l_NJlowTM_250lessMETlessInfBTheavyUp",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfPUdown","SR1l_NJlowTM_250lessMETlessInfPUdown",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfPUup","SR1l_NJlowTM_250lessMETlessInfPUup",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfPDFdown","SR1l_NJlowTM_250lessMETlessInfPDFdown",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfPDFup","SR1l_NJlowTM_250lessMETlessInfPDFup",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfalphaSdown","SR1l_NJlowTM_250lessMETlessInfalphaSdown",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfalphaSup","SR1l_NJlowTM_250lessMETlessInfalphaSup",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfQ2down","SR1l_NJlowTM_250lessMETlessInfQ2down",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJlowTM_250lessMETlessInfQ2up","SR1l_NJlowTM_250lessMETlessInfQ2up",&SR1l_NJlowTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInf","SR1l_NJmidTM_250lessMETlessInf",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfLSFdown","SR1l_NJmidTM_250lessMETlessInfLSFdown",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfLSFup","SR1l_NJmidTM_250lessMETlessInfLSFup",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfBTlightDown","SR1l_NJmidTM_250lessMETlessInfBTlightDown",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfBTlightUp","SR1l_NJmidTM_250lessMETlessInfBTlightUp",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfBTheavyDown","SR1l_NJmidTM_250lessMETlessInfBTheavyDown",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfBTheavyUp","SR1l_NJmidTM_250lessMETlessInfBTheavyUp",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfPUdown","SR1l_NJmidTM_250lessMETlessInfPUdown",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfPUup","SR1l_NJmidTM_250lessMETlessInfPUup",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfPDFdown","SR1l_NJmidTM_250lessMETlessInfPDFdown",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfPDFup","SR1l_NJmidTM_250lessMETlessInfPDFup",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfalphaSdown","SR1l_NJmidTM_250lessMETlessInfalphaSdown",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfalphaSup","SR1l_NJmidTM_250lessMETlessInfalphaSup",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfQ2down","SR1l_NJmidTM_250lessMETlessInfQ2down",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJmidTM_250lessMETlessInfQ2up","SR1l_NJmidTM_250lessMETlessInfQ2up",&SR1l_NJmidTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInf","SR1l_NJhighTM_250lessMETlessInf",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfLSFdown","SR1l_NJhighTM_250lessMETlessInfLSFdown",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfLSFup","SR1l_NJhighTM_250lessMETlessInfLSFup",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfBTlightDown","SR1l_NJhighTM_250lessMETlessInfBTlightDown",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfBTlightUp","SR1l_NJhighTM_250lessMETlessInfBTlightUp",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfBTheavyDown","SR1l_NJhighTM_250lessMETlessInfBTheavyDown",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfBTheavyUp","SR1l_NJhighTM_250lessMETlessInfBTheavyUp",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfPUdown","SR1l_NJhighTM_250lessMETlessInfPUdown",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfPUup","SR1l_NJhighTM_250lessMETlessInfPUup",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfPDFdown","SR1l_NJhighTM_250lessMETlessInfPDFdown",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfPDFup","SR1l_NJhighTM_250lessMETlessInfPDFup",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfalphaSdown","SR1l_NJhighTM_250lessMETlessInfalphaSdown",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfalphaSup","SR1l_NJhighTM_250lessMETlessInfalphaSup",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfQ2down","SR1l_NJhighTM_250lessMETlessInfQ2down",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJhighTM_250lessMETlessInfQ2up","SR1l_NJhighTM_250lessMETlessInfQ2up",&SR1l_NJhighTM_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInf","SR1l_NJ_250lessMETlessInf",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfLSFdown","SR1l_NJ_250lessMETlessInfLSFdown",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfLSFup","SR1l_NJ_250lessMETlessInfLSFup",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfBTlightDown","SR1l_NJ_250lessMETlessInfBTlightDown",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfBTlightUp","SR1l_NJ_250lessMETlessInfBTlightUp",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfBTheavyDown","SR1l_NJ_250lessMETlessInfBTheavyDown",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfBTheavyUp","SR1l_NJ_250lessMETlessInfBTheavyUp",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfPUdown","SR1l_NJ_250lessMETlessInfPUdown",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfPUup","SR1l_NJ_250lessMETlessInfPUup",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfPDFdown","SR1l_NJ_250lessMETlessInfPDFdown",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfPDFup","SR1l_NJ_250lessMETlessInfPDFup",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfalphaSdown","SR1l_NJ_250lessMETlessInfalphaSdown",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfalphaSup","SR1l_NJ_250lessMETlessInfalphaSup",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfQ2down","SR1l_NJ_250lessMETlessInfQ2down",&SR1l_NJ_250lessMETlessInf);
AddRegion("SR1l_NJ_250lessMETlessInfQ2up","SR1l_NJ_250lessMETlessInfQ2up",&SR1l_NJ_250lessMETlessInf);
    // ------------------
    // Channels
    // ------------------
    
    AddChannel("lepChannel","lepChannel", &lepChannel);

    SetLumi(36.46);

    Create1DHistos();

    fillYieldsVector();
    WriteXMLConfig(); 
}

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
    counter++;
    nthentry++;

    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);


    checkNegativeYields = false;

    if(nthentry == myEvent.nentries)
    {
        nthentry=0;
        cout << "seting true to check negative " << endl;
        checkNegativeYields = true;
    }
    vector<string> classLabels;
    GetProcessClassLabelList(&classLabels);

    myEvent.trigger = CheckTrigger( myEvent.is_data, currentDataset);
    if( ( currentProcessClass == "ttZ" || currentProcessClass == "ZZWZ" || currentProcessClass == "ZZ"|| currentProcessClass == "WZ" ) && !(myEvent.isZtoNuNu) )
         currentProcessClass = "";

    float weightLumi = getWeight(currentProcessType, GetLumi()); //@MJ@ TODO cross section form file?!

    //@MJ@ TODO I hate myself for this, but no better solution foud
    //computation of up/down weights
    //NOTICE, important is to fill only weight histo and only have 1 process class

    float weight_pt = 1;
    if( currentProcessClass == "ttZ")
        weight_pt = reweightTop(myEvent.top_pt, myEvent.atop_pt);
 
    vector<float> weightV;
    vector<float> weightVTot;
    weightV.clear();
    float nEvents =  myEvent.wNormalization.at(22);
    float SF = 1.37;
    //for number of SR
    for(uint32_t SR=0; SR<9; SR++) //@MJ@ TODO nr of sig regions changes
    {
        float w = 0;
        //normal
        weightV.push_back(weightLumi*SF); //@MJ@ TODO do not forget SF
        //LSFdown
        if(counter == 1) statnames << "lepSFDN" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF_down*( nEvents / myEvent.wNormalization.at(30) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ));
        weightV.push_back(w);
        //LSFup
        if(counter == 1) statnames << "lepSFUP" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF_up*( nEvents / myEvent.wNormalization.at(29) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ));
        weightV.push_back(w);
        //BTlightdown
        if(counter == 1) statnames << "btagLightDN" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_light_DN*( nEvents / myEvent.wNormalization.at(18) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTlightup
        if(counter == 1) statnames << "btagLightUP" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_light_UP*( nEvents / myEvent.wNormalization.at(16) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTheabydown
        if(counter == 1) statnames << "btagHeavyDN" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_heavy_DN*( nEvents / myEvent.wNormalization.at(17) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //BTheavyup
        if(counter == 1) statnames << "btagHeavyUP" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_btagsf_heavy_UP*( nEvents / myEvent.wNormalization.at(15) ) * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //PUdown
        if(counter == 1) statnames << "PUdown" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) * myEvent.weight_PUdown; //@MJ@ TODO PU without any normalization?!
        weightV.push_back(w);
        //PUup
        if(counter == 1) statnames << "PUup"<< endl;
        w = SF * GetLumi() *  myEvent.scale1fb  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) * myEvent.weight_PUup;
        weightV.push_back(w);
        //PDFdown
        if(counter == 1) statnames << "pdfDN" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs((myEvent.pdf_down_weight/myEvent.genweights->at(0)) * (  myEvent.wNormalization.at(1) / myEvent.wNormalization.at(11) ))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        //cout << "pdfd  " << pdfd << " w " << w << endl;
        weightV.push_back(w);
        //PDFup
        if(counter == 1) statnames << "pdfUP" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )  * abs((myEvent.pdf_up_weight/myEvent.genweights->at(0)) * (  myEvent.wNormalization.at(1)/ myEvent.wNormalization.at(10) ))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        weightV.push_back(w);
        //alphaSdown
        if(counter == 1) statnames << "alphaSDN" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * abs(myEvent.weight_alphas_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(13)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) ; //TODO
        weightV.push_back(w);
        //alphaSup
        if(counter == 1) statnames << "alphaSUP" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs(myEvent.weight_alphas_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(12)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ) ; //TODO
        weightV.push_back(w);
        //Q2down
        if(counter == 1) statnames << "Q2DN" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb  * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) )* abs(myEvent.weight_q2_down*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(9)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) )  ; //TODO
        weightV.push_back(w);
        //Q2up
        if(counter == 1) statnames << "Q2UP" << endl;
        w = SF * GetLumi() *  myEvent.scale1fb *  myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28)) * abs(myEvent.weight_q2_up*( myEvent.wNormalization.at(1) / myEvent.wNormalization.at(5)))* myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) ); //TODO
        weightV.push_back(w);
    }

    vector<string> theReg;
    GetRegionTagList(&theReg);
    vector<string> theVar;
    GetVariablesTagList(&theVar);
    cout << "size of variables " << theVar.size() << endl;
    
    for(uint32_t v=0; v<theVar.size(); v++)
    {
        for(uint32_t s=0; s< weightV.size(); s++)
        weightVTot.push_back(weightV.at(s));
    }


    if( weightVTot.size() != theReg.size()*theVar.size())
        throw std::runtime_error("vector of weights does not have same size as regions, the weights will not be correctly assesed");

   //cout << "scale 1fb " << myEvent.scale1fb << " kfactor " << myEvent.kfactor << " xsec " << myEvent.crossSection << " 1/weight " << (myEvent.scale1fb)/(myEvent.kfactor*myEvent.crossSection) << " pdf down " << abs(myEvent.pdf_down_weight* ( nEvents / myEvent.wNormalization.at(11) )) << " pdf up " << abs(myEvent.pdf_up_weight*( nEvents / myEvent.wNormalization.at(10) )) << " gen[0]" << myEvent.genweights->at(0) << endl; 
 
    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    AutoFillProcessClass(currentProcessClass, weight, checkNegativeYields, weightVTot, true);


    if(counter % 10000 == 0)
    {
        cout << counter << endl;
    }

}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    statnames.close();

    SetGlobalStringOption("Plot", "infoTopRight","#bf{CMS} #it{Preliminary}" );
    SetGlobalStringOption("Plot", "infoTopLeft","36.5^{-1} (13 TeV)" );

    // Make and write the plots
    SchedulePlots("1DSuperimposedNoNorm");
    SchedulePlots("1DDataMCComparison");

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./"+outputName+"plotsRatios/");

    // ######################
    //  Tables and other stuff
    // ######################

vector<string> totYield = { "SR1l_AB_250lessMETlessInf" , "SR1l_AB_250lessMETlessInfLSFdown" , "SR1l_AB_250lessMETlessInfLSFup" , "SR1l_AB_250lessMETlessInfBTlightDown" , "SR1l_AB_250lessMETlessInfBTlightUp" , "SR1l_AB_250lessMETlessInfBTheavyDown" , "SR1l_AB_250lessMETlessInfBTheavyUp" , "SR1l_AB_250lessMETlessInfPUdown" , "SR1l_AB_250lessMETlessInfPUup" , "SR1l_AB_250lessMETlessInfPDFdown" , "SR1l_AB_250lessMETlessInfPDFup" , "SR1l_AB_250lessMETlessInfalphaSdown" , "SR1l_AB_250lessMETlessInfalphaSup" , "SR1l_AB_250lessMETlessInfQ2down" , "SR1l_AB_250lessMETlessInfQ2up" , "SR1l_CD_250lessMETlessInf" , "SR1l_CD_250lessMETlessInfLSFdown" , "SR1l_CD_250lessMETlessInfLSFup" , "SR1l_CD_250lessMETlessInfBTlightDown" , "SR1l_CD_250lessMETlessInfBTlightUp" , "SR1l_CD_250lessMETlessInfBTheavyDown" , "SR1l_CD_250lessMETlessInfBTheavyUp" , "SR1l_CD_250lessMETlessInfPUdown" , "SR1l_CD_250lessMETlessInfPUup" , "SR1l_CD_250lessMETlessInfPDFdown" , "SR1l_CD_250lessMETlessInfPDFup" , "SR1l_CD_250lessMETlessInfalphaSdown" , "SR1l_CD_250lessMETlessInfalphaSup" , "SR1l_CD_250lessMETlessInfQ2down" , "SR1l_CD_250lessMETlessInfQ2up" , "SR1l_EFGH_250lessMETlessInf" , "SR1l_EFGH_250lessMETlessInfLSFdown" , "SR1l_EFGH_250lessMETlessInfLSFup" , "SR1l_EFGH_250lessMETlessInfBTlightDown" , "SR1l_EFGH_250lessMETlessInfBTlightUp" , "SR1l_EFGH_250lessMETlessInfBTheavyDown" , "SR1l_EFGH_250lessMETlessInfBTheavyUp" , "SR1l_EFGH_250lessMETlessInfPUdown" , "SR1l_EFGH_250lessMETlessInfPUup" , "SR1l_EFGH_250lessMETlessInfPDFdown" , "SR1l_EFGH_250lessMETlessInfPDFup" , "SR1l_EFGH_250lessMETlessInfalphaSdown" , "SR1l_EFGH_250lessMETlessInfalphaSup" , "SR1l_EFGH_250lessMETlessInfQ2down" , "SR1l_EFGH_250lessMETlessInfQ2up" , "SR1l_I_250lessMETlessInf" , "SR1l_I_250lessMETlessInfLSFdown" , "SR1l_I_250lessMETlessInfLSFup" , "SR1l_I_250lessMETlessInfBTlightDown" , "SR1l_I_250lessMETlessInfBTlightUp" , "SR1l_I_250lessMETlessInfBTheavyDown" , "SR1l_I_250lessMETlessInfBTheavyUp" , "SR1l_I_250lessMETlessInfPUdown" , "SR1l_I_250lessMETlessInfPUup" , "SR1l_I_250lessMETlessInfPDFdown" , "SR1l_I_250lessMETlessInfPDFup" , "SR1l_I_250lessMETlessInfalphaSdown" , "SR1l_I_250lessMETlessInfalphaSup" , "SR1l_I_250lessMETlessInfQ2down" , "SR1l_I_250lessMETlessInfQ2up" , "SR1l_CDEFGH_250lessMETlessInf" , "SR1l_CDEFGH_250lessMETlessInfLSFdown" , "SR1l_CDEFGH_250lessMETlessInfLSFup" , "SR1l_CDEFGH_250lessMETlessInfBTlightDown" , "SR1l_CDEFGH_250lessMETlessInfBTlightUp" , "SR1l_CDEFGH_250lessMETlessInfBTheavyDown" , "SR1l_CDEFGH_250lessMETlessInfBTheavyUp" , "SR1l_CDEFGH_250lessMETlessInfPUdown" , "SR1l_CDEFGH_250lessMETlessInfPUup" , "SR1l_CDEFGH_250lessMETlessInfPDFdown" , "SR1l_CDEFGH_250lessMETlessInfPDFup" , "SR1l_CDEFGH_250lessMETlessInfalphaSdown" , "SR1l_CDEFGH_250lessMETlessInfalphaSup" , "SR1l_CDEFGH_250lessMETlessInfQ2down" , "SR1l_CDEFGH_250lessMETlessInfQ2up" , "SR1l_NJlowTM_250lessMETlessInf" , "SR1l_NJlowTM_250lessMETlessInfLSFdown" , "SR1l_NJlowTM_250lessMETlessInfLSFup" , "SR1l_NJlowTM_250lessMETlessInfBTlightDown" , "SR1l_NJlowTM_250lessMETlessInfBTlightUp" , "SR1l_NJlowTM_250lessMETlessInfBTheavyDown" , "SR1l_NJlowTM_250lessMETlessInfBTheavyUp" , "SR1l_NJlowTM_250lessMETlessInfPUdown" , "SR1l_NJlowTM_250lessMETlessInfPUup" , "SR1l_NJlowTM_250lessMETlessInfPDFdown" , "SR1l_NJlowTM_250lessMETlessInfPDFup" , "SR1l_NJlowTM_250lessMETlessInfalphaSdown" , "SR1l_NJlowTM_250lessMETlessInfalphaSup" , "SR1l_NJlowTM_250lessMETlessInfQ2down" , "SR1l_NJlowTM_250lessMETlessInfQ2up" , "SR1l_NJmidTM_250lessMETlessInf" , "SR1l_NJmidTM_250lessMETlessInfLSFdown" , "SR1l_NJmidTM_250lessMETlessInfLSFup" , "SR1l_NJmidTM_250lessMETlessInfBTlightDown" , "SR1l_NJmidTM_250lessMETlessInfBTlightUp" , "SR1l_NJmidTM_250lessMETlessInfBTheavyDown" , "SR1l_NJmidTM_250lessMETlessInfBTheavyUp" , "SR1l_NJmidTM_250lessMETlessInfPUdown" , "SR1l_NJmidTM_250lessMETlessInfPUup" , "SR1l_NJmidTM_250lessMETlessInfPDFdown" , "SR1l_NJmidTM_250lessMETlessInfPDFup" , "SR1l_NJmidTM_250lessMETlessInfalphaSdown" , "SR1l_NJmidTM_250lessMETlessInfalphaSup" , "SR1l_NJmidTM_250lessMETlessInfQ2down" , "SR1l_NJmidTM_250lessMETlessInfQ2up" , "SR1l_NJhighTM_250lessMETlessInf" , "SR1l_NJhighTM_250lessMETlessInfLSFdown" , "SR1l_NJhighTM_250lessMETlessInfLSFup" , "SR1l_NJhighTM_250lessMETlessInfBTlightDown" , "SR1l_NJhighTM_250lessMETlessInfBTlightUp" , "SR1l_NJhighTM_250lessMETlessInfBTheavyDown" , "SR1l_NJhighTM_250lessMETlessInfBTheavyUp" , "SR1l_NJhighTM_250lessMETlessInfPUdown" , "SR1l_NJhighTM_250lessMETlessInfPUup" , "SR1l_NJhighTM_250lessMETlessInfPDFdown" , "SR1l_NJhighTM_250lessMETlessInfPDFup" , "SR1l_NJhighTM_250lessMETlessInfalphaSdown" , "SR1l_NJhighTM_250lessMETlessInfalphaSup" , "SR1l_NJhighTM_250lessMETlessInfQ2down" , "SR1l_NJhighTM_250lessMETlessInfQ2up" , "SR1l_NJ_250lessMETlessInf" , "SR1l_NJ_250lessMETlessInfLSFdown" , "SR1l_NJ_250lessMETlessInfLSFup" , "SR1l_NJ_250lessMETlessInfBTlightDown" , "SR1l_NJ_250lessMETlessInfBTlightUp" , "SR1l_NJ_250lessMETlessInfBTheavyDown" , "SR1l_NJ_250lessMETlessInfBTheavyUp" , "SR1l_NJ_250lessMETlessInfPUdown" , "SR1l_NJ_250lessMETlessInfPUup" , "SR1l_NJ_250lessMETlessInfPDFdown" , "SR1l_NJ_250lessMETlessInfPDFup" , "SR1l_NJ_250lessMETlessInfalphaSdown" , "SR1l_NJ_250lessMETlessInfalphaSup" , "SR1l_NJ_250lessMETlessInfQ2down" , "SR1l_NJ_250lessMETlessInfQ2up"};


    TableDataMC(this, totYield,"lepChannel",  "includeSignal" ).Print(outputName+".tab", 6);
    TableDataMC(this, totYield,"lepChannel", "includeSignal" ).PrintLatex(outputName+".tex", 6);

    ofstream sigfile("signalRegMorPlots.txt");
    if (sigfile.is_open())
    {
        for(uint32_t r=0; r<totYield.size(); r++)
        {
            sigfile << totYield.at(r) << endl;
        }
            sigfile.close();
    }


    cout << "end of processing" << endl;
 }


    float getWeight(string currentProcessType, float lumi)
    {
        float nEvents =  myEvent.wNormalization.at(22);
        float all_weights = lumi*  myEvent.scale1fb * myEvent.weight_lepSF*( nEvents / myEvent.wNormalization.at(28) ) * myEvent.weight_vetoLepSF*( nEvents / myEvent.wNormalization.at(31) );
        if(currentProcessType == "signal")
             throw std::runtime_error("weight for signal still waitning to be implemented!");
        return all_weights;
    }

    
    float reweightTop(float topPt, float atopPt)
    {
        return sqrt(exp(0.0615-(0.0005*topPt))*exp(0.0615-(0.0005*atopPt))); 
    }
