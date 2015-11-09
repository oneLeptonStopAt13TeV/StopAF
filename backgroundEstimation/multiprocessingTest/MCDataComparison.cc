#include <vector>
#include <iostream>
using namespace std;

#include "TTree.h"
//babyEvent myEvent;




//-----------------------------------
//-- Load only necessary branches
//-----------------------------------

//#define USE_GEN_INFO
//#define USE_GEN_INFO_EXT
//#define USE_SKIMMING_VAR
#define USE_JETS
//#define USE_JETS_EXT
//#define USE_GEN_INFO
//#define USE_GEN_INFO_EXT
#define USE_LEP1
//#define USE_LEP1_EXT
#define USE_LEP2
//#define USE_LEP2_EXT
#define USE_JETS
//#define USE_JETS_EXT
#define USE_PV
#define USE_WEIGHTS
#define USE_VAR_BASELINE
#define USE_GLOBAL_VAR
//#define USE_OLD_VAR


//-----------------------------------
// -- load Babytuple Format --
// -- Need first to activate the branches
//-----------------------------------
//#include "../common/BabyTupleFormat.h"
#include "../../common/Reader_CommonFormat.h"

#include "../../common/selectionDefinitions13TeV.h"


//---------------------------------------
// -- Add variables computed on the flyt
// -- Need first to load the format
//--------------------------------------
//#include "OnTheFlyVariables.h"

// Should be called after selectionDefinition 'cause it need myEvent varible ...
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"



bool muonChannelSelector() { return true; }

// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");
    cout << "non correct version" << endl;
    //babyTuplePath = "./";
    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Nov1st_missingInfo/";
    
    totalNumberOfWorkers = 8;

    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 300,  &(myEvent.mt_met_lep));
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 300,  &(myEvent.pfmet));
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W));
    AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets));
    AddVariable("Nbjets", "N_{b-jets}",  "", 6,   0, 5,  &(myEvent.ngoodbtags));
    //AddVariable("SRBins", "Signal Region box",  "", 11,   -5, 5,  &(onTheFlyVariables.SRbins));
    
    
    //AddVariable("x", "whatever2", "", 20, -40,  40,  &(onTheFlyVariables.someCoolVariable));
    // ...

    AddProcessClass("SingleTop", "Single top", "background", kGreen);
    	AddDataset("ST_s", "SingleTop" );
    	AddDataset("ST_t-atop", "SingleTop" );
    	AddDataset("ST_t-top", "SingleTop" );
    	AddDataset("ST_tW-atop", "SingleTop" );
    	AddDataset("ST_tW-top", "SingleTop" );

    AddProcessClass("ttV", "tt+boson", "background", kMagenta);
    	AddDataset("TTW_ln", "ttV");
    	AddDataset("TTW_qq", "ttV");
    	AddDataset("TTZ_ll", "ttV");
    	AddDataset("TTZ_qq", "ttV");

    AddProcessClass("VV", "di-bosons", "background", kMagenta+1);
    	AddDataset("WW_aMC", "VV");
    	AddDataset("ZZ", "VV");
    	AddDataset("WZ", "VV");
    
    AddProcessClass("TTJets", "TT+jets", "background", kRed);
    	AddDataset("TTjets_M5", "TTJets");
   
    AddProcessClass("WJets", "W+jets", "background", kBlue);
    	AddDataset("Wjets_aMC", "WJets");
    
    AddProcessClass("DATA", "DATA", "data", kRed);
    	//AddDataset("SingleMuon_RunD", "DATA");
    	//AddDataset("SingleElectron_RunD", "DATA");
    	AddDataset("SingleMuon_RunD_prompt", "DATA");
    	AddDataset("SingleElectron_RunD_prompt", "DATA");
    
    // ...

    AddRegion("preselection", "Preselection", &goesInPreselection);
    AddRegion("BaselineSelection", "BaselineSelection", &goesInBaselineSearchSR);
    // ...

    AddChannel("def", "#mu channel", &muonChannelSelector);
    AddChannel("muon", "#mu channel", &goesInSingleMuonChannel);
    AddChannel("electron", "electron channel", &goesInSingleElecChannel);
    AddChannel("electronEC", "electron (EC) channel", &goesInSingleElecECChannel);
    AddChannel("singleLepton", "single lepton channel", &goesInSingleLeptonChannel);
    // ...

    //SetLumi(1226);
    SetLumi(578.3);

    Create1DHistos();
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{

   /*
   cout << "testing run" << endl;
   cout << "nr. of selected leptons: " << myEvent.numberOfSelectedLeptons << endl;
   cout << "nr. of init events: " << myEvent.totalNumberOfInitialEvent << endl;


    cout<<myEvent.mt_met_lep<<endl;
    cout<<"sel: "<<goesInPreselection()<<endl;
    */
    myEvent.crossSection  = CrossSection(currentDataset);
    //myEvent.totalNumberOfInitialEvent = 
    /*
    cout<<"xsextion : "<<myEvent.crossSection<<endl;
    cout<<"tot weight: "<<myEvent.totalNumberOfInitialEvent<<endl;
    cout<<"weight: "<<myEvent.mc_weight<<endl;
    cout<<"met: "<<myEvent.pfmet<<endl;
    */
    // Compute on the fly variables if needed

    //ComputeOnTheFlyVariables();

    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    // Compute weight for current event

    double weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;

    double weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;
    //if (currentProcessClass == "WJets") weight *= 2.5;

    // Fill this event in the histo collections
    //cout<<"weight = "<<weight<<" "<<myEvent.crossSection<<" "<<GetLumi()<<" "<< myEvent.mc_weight<<" "<<myEvent.totalNumberOfInitialEvent<<endl;
    AutoFillProcessClass(currentProcessClass, weight);
    //AutoFillProcessClass(currentProcessClass, 0.001);
}

// ################################################################

void BabyScrewdriver::PostProcessingStep()
{
    // ######################
    //  Plot configuration and production
    // ######################

    // Schedule plots

    SchedulePlots("1DSuperimposed");
    SchedulePlots("1DStack");
    SchedulePlots("1DDataMCComparison");
    
    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 8 TeV");

    SetGlobalBoolOption("Plot", "exportPdf", true);
    SetGlobalBoolOption("Plot", "exportEps", false);
    SetGlobalBoolOption("Plot", "exportPng", false);

    // Make and write the plots

    cout << endl;
    cout << "   > Making plots..." << endl;
    MakePlots();
    cout << "   > Saving plots..." << endl;
    WritePlots("./plots/");

    // ######################
    //  Tables and other stuff
    // ######################

    TableDataMC(this,{"preselection"},"muon", "").Print("table.tab",4);

}



