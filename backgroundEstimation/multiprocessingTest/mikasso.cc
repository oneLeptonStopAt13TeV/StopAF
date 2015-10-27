#include <vector>
#include <iostream>
using namespace std;

#include "TTree.h"
//babyEvent myEvent;




//-----------------------------------
//-- Load only necessary branches
//-----------------------------------

#define USE_GEN_INFO
#define USE_GEN_INFO_EXT
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
#define USE_OLD_VAR
#define USE_AK8_JETS


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
#include "OnTheFlyVariables.h"

// Should be called after selectionDefinition 'cause it need myEvent varible ...
#include "../../sonicScrewdriver/interface/BabyScrewdriver.h"



bool muonChannelSelector() { return true; }
uint64_t counting;
// ################################################################

void BabyScrewdriver::Init()
{
    PrintBoxedMessage("Initializing babyScrewdriver");

    babyTuplePath = "/opt/sbg/scratch1/cms/echabert/store/babyTuples/Oct21_postpatch8";
    //babyTuplePath = "./";
    
    totalNumberOfWorkers = 1;



    AddVariable("MT", "M_T [GeV]",  "", 50,   0, 300,  &(myEvent.mt_met_lep));
    AddVariable("MET", "MET [GeV]",  "", 50,   0, 300,  &(myEvent.pfmet));
    AddVariable("MT2W", "MT2W",  "", 50,   0, 300,  &(myEvent.MT2W));
    AddVariable("Njets", "N_{jets}",  "", 11,   0, 10,  &(myEvent.ngoodjets));
    AddVariable("Nbjets", "N_{b-jets}",  "", 6,   0, 5,  &(myEvent.ngoodbtags));
    AddVariable("SRBins", "Signal Region box",  "", 11,   -5, 5,  &(onTheFlyVariables.SRbins));
    AddVariable("genWPt", "Pt of gen W",  "", 500,   0, 1000,  &(onTheFlyVariables.m_genWPt));
    AddVariable("genTopPt", "Pt of gen top quark",  "", 500,   0, 1000,  &(onTheFlyVariables.m_genTopPt));
    AddVariable("genWdR", "#delta R related to gen W",  "", 100,   0, 4,  &(onTheFlyVariables.m_genWdR));
    AddVariable("genTopdR", "#deltar R related to gen top quark",  "", 100,   0, 4,  &(onTheFlyVariables.m_genTopdR));
    AddVariable("recWPt", "Pt of rec W",  "", 500,   0, 500,  &(onTheFlyVariables.m_recWPt));
    AddVariable("ak8WMass", "mass of W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WMass));
    AddVariable("ak8WPrunedMass", "pruned mass of W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WPrunedMass));
    AddVariable("ak8WTrimmedMass", "trimmed mass of W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WTrimmedMass));
    AddVariable("ak8WSoftDropMass", "softdrop mass of W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WSoftDropMass));
    AddVariable("ak8recoWPt", "pt of W ak8 jet",  "", 200,   0, 1000,  &(onTheFlyVariables.m_ak8recoWPt));
    AddVariable("ak8WFakeMass", "mass of fake W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WFakeMass));
    AddVariable("ak8WFakePrunedMass", "pruned mass of fake W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WFakePrunedMass));
    AddVariable("ak8WFakeTrimmedMass", "trimmed mass of fake W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WFakeTrimmedMass));
    AddVariable("ak8WFakeSoftDropMass", "softdrop mass of fake W ak8 jet",  "", 100,   0, 400,  &(onTheFlyVariables.m_ak8WFakeSoftDropMass));
    AddVariable("ak8recoWFakePt", "pt of fake W ak8 jet",  "", 200,   0, 1000,  &(onTheFlyVariables.m_ak8recoWFakePt));
       
    // ...

    AddProcessClass("TTJets-FXFX", "TTJets-FXFX", "background", kPink);
        AddDataset("TTJets-FXFX", "TTJets-FXFX", 0, 0);
    // ...

    AddRegion("preselection", "Preselection", &goesInPreselection);
    AddRegion("BaselineSelection", "BaselineSelection", &goesInBaselineSearchSR);
    // ...

    AddChannel("muon", "#mu channel", &muonChannelSelector);
    // ...

    SetLumi(2170.);  //@MJ@ TODO correct lumi

    Create1DHistos();
    Add2DHisto("genWPt","genWdR");
    Add2DHisto("genTopPt","genTopdR");
}

// ################################################################

void BabyScrewdriver::ActionForEachEvent(string currentDataset)
{
     
    //cout<<myEvent.mt_met_lep<<endl;
    //cout << " sel: "<< goesInPreselection() << endl;
    myEvent.crossSection  = CrossSection(currentDataset);
    //cout<<"xsextion : "<<myEvent.crossSection<<endl;
    //cout<<"tot weight: "<<myEvent.totalNumberOfInitialEvent<<endl;
    //cout<<"weight: "<<myEvent.mc_weight<<endl;
    // Compute on the fly variables if needed

    ComputeOnTheFlyVariables();
    
    //find generated particles and their properties
    //Initial clean up
    onTheFlyVariables.m_genWPt = 0;
    onTheFlyVariables.m_genWdR = 0;
    onTheFlyVariables.m_genTopPt = 0;
    onTheFlyVariables.m_genTopdR = 0;
    
    //find W decaying to two quarks
    findGenParticleProps(24, &(onTheFlyVariables.m_genWPt), &(onTheFlyVariables.m_genWdR)); //W
    //cout << "1) back in main; W: gen pt is: " << onTheFlyVariables.m_genWPt << " gen dR is: " <<  onTheFlyVariables.m_genWdR << " end of one event" << endl;
    
    //find t decaying to W and quark
    findGenParticleProps(6, &(onTheFlyVariables.m_genTopPt), &(onTheFlyVariables.m_genTopdR)); //top
    //cout << "2) back in main; top:  gen pt is: " << onTheFlyVariables.m_genTopPt << " gen dR is: " <<  onTheFlyVariables.m_genTopdR << " end of one event" << endl;

    //count efficiency of W tagging
    countEfficiency(24);
    //cout << "W ak8 pt" << onTheFlyVariables.m_ak8recoWPt << "ake ak8 W pt" << onTheFlyVariables.m_ak8recoWFakePt << endl;

    // Determine which processClass to fill
    // (in the most trivial case, only call GetProcessClass(currentDataset),
    // but you might want to split a dataset according to
    // the number of generated leptons, for instance)
    string currentProcessClass = GetProcessClass(currentDataset);
    string currentProcessType  = GetProcessClassType(currentProcessClass);

    // Compute weight for current event

    float weightLumi = myEvent.crossSection * GetLumi() * myEvent.mc_weight / myEvent.totalNumberOfInitialEvent;

    float weight     = weightLumi;
    if (currentProcessType == "data") weight = 1.0;

    // Fill this event in the histo collections

    AutoFillProcessClass(currentProcessClass, weight);
    //AutoFillProcessClass(currentProcessClass, 0.001);

    counting++;
    if(counting % 1000 == 0)
    {
        cout << counting << endl;
    }
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
    SchedulePlots("2D");

    // Config plots

    SetGlobalStringOption("Plot", "infoTopRight", "CMS Simulation");
    SetGlobalStringOption("Plot", "infoTopLeft",  "#sqrt{s} = 13 TeV");

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



