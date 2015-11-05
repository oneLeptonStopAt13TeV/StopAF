#define MT_CUT    150
#define MET_CUT   250
#define MT_PSCUT    100
#define MET_PSCUT   80
#define NJET_CUT  4
#define NBJET_CUT 1
#define NLEP_CUT  1
#define MET_CUTLL 50
#define MTW2_CUT 200
#define MLB_CUT 175
#define ETMISS_BOUND_250 250
#define ETMISS_BOUND_300 300
#define ETMISS_BOUND_325 325
#define ETMISS_BOUND_375 375
#define ETMISS_BOUND_400 400

// Not sure that it is a good idea to include this here,
// since one often needs to use a modified format because
// of skimming or use of tiny tuples...

//#include "Reader.h"
//#include "Reader_final.h"  // Has the extended BDT info defined
#include "Reader_CommonFormat.h"  // Has the extended BDT info defined
#include "CrossSection.h" // used to compute the weight

#include <iostream>
#include <fstream>

using namespace std;

// NB : When you call any of the following functions,
// these three variables must be filled with the current
// event being taken care of, the name of the sample and
// the sample type (background, signal or data)

babyEvent myEvent;
string sampleName;
string sampleType;

// MT cuts definitions
// ###################

bool goesInMTpeak()     { if ((myEvent.mt_met_lep > 50) && (myEvent.mt_met_lep < 80)) return true; else return false; }
bool goesInMTtail()     { if (myEvent.mt_met_lep > MT_CUT)                    return true; else return false; }
bool goesInMTinverted() { if (myEvent.mt_met_lep < MT_CUT)                    return true; else return false; }

// Control region definitions
// ##########################

bool goesInPreVetoSelection()
{

    if (myEvent.pfmet < MET_PSCUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    //if (myEvent.ngoodbtags < NBJET_CUT)  return false;

    return true;
}

bool goesInPreVetoSelectionMTtail()     { return (goesInPreVetoSelection() && goesInMTtail());     }
bool goesInPreVetoSelectionMTpeak()     { return (goesInPreVetoSelection() && goesInMTpeak());     }
bool goesInPreVetoSelectionMTinverted() { return (goesInPreVetoSelection() && goesInMTinverted()); }


bool goesInPreselectionNoVeto()
{
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;

    return true;
}

bool goesInPreselectionNoVetoNoMetCut()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;

    return true;
}

bool goesInPreselection()
{
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

    return true;
}

bool goesInBaselineSearchSR() {
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    
    if ((!myEvent.PassTrackVeto) || (!myEvent.PassTauVeto)) return false;

return (goesInPreselection() && goesInMTtail() );}
bool goesInBaselineSearchSR2b() { return (goesInPreselection() && goesInMTtail() && myEvent.ngoodbtags >=2 );}
//bool goesInLargeDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.MT2W > 200 && myEvent.dphi_Wlep> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}
//bool goesInSmallDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.dphi_Wlep> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}

/*
bool goesInLargeDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.MT2W > 200 && myEvent.minDPhi_jmet> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}
bool goesInSmallDMSR() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.minDPhi_jmet> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 200);}
bool goesInSmallDMSR300() { return (goesInPreselection() &&  myEvent.mt_met_lep > 150 && myEvent.minDPhi_jmet> 0.8 && myEvent.hadronic_top_chi2 < 10 && myEvent.pfmet > 300);}
bool goesInSmallDMSR2b() { return (goesInSmallDMSR() && myEvent.ngoodbtags>=2);}
*/

// Categorization of Events
// ########################

bool NoAk8Jets()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() != 0) return false;

    return true;
}

bool AtLeastOneAk8Jet()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() == 0) return false;

    return true;
}

bool NoAk10Jets()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak10pfjets_mass.size() != 0) return false; //@MJ@ TODO uncomment this one day

    return true;
}

bool AtLeastOneAk10Jet()
{
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak10pfjets_mass.size() == 0) return false; //@MJ@ TODO uncomment this one day

    return true;
}

bool ThreeJetsInEventOneOfThemAk8()
{
    if (myEvent.ngoodjets < 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    if (myEvent.ak8pfjets_mass.size() != 0) return false;

    return true;
}

bool ThreeJetsInEventOneOfThemAk10()
{
    if (myEvent.ngoodjets < 3)  return false;
    if (myEvent.MT2W < MTW2_CUT)  return false;
    //if (myEvent.ak10pfjets_mass.size() != 0) return false; //@MJ@ TODO uncomment this one day

    return true;
}

// New signal regions
// ###################

bool SRbin1(ofstream &outputFile)
{
    static uint32_t j;
    j++;

    if(j == 1)
    {
        outputFile << "Signal region bin 1" << endl;
        outputFile << "----------------------------" << endl;
        outputFile << "myEvent.ngoodleps == " << NLEP_CUT << endl;
        outputFile << "myEvent.ngoodjets > " << NJET_CUT << endl;
        outputFile << "myEvent.ngoodbtags > " << NBJET_CUT << endl;
        outputFile << "myEvent.Mlb < " << MLB_CUT << endl;
        outputFile << "myEvent.MT2W < " << MTW2_CUT << endl;
        //outputFile <<  ETMISS_BOUND_250 " < myEvent.ETmiss <
    } 
    
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.ETmiss < ETMISS_BOUND_250)  return false;
    if (myEvent.ETmiss > ETMISS_BOUND_325)  return false;

    return true;
}

bool SRbin2()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.ETmiss <= ETMISS_BOUND_325)  return false;

    return true;
}

bool SRbin3()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ETmiss < ETMISS_BOUND_250)  return false;
    if (myEvent.ETmiss > ETMISS_BOUND_375)  return false;

    return true;
}

bool SRbin4()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb > MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ETmiss <= ETMISS_BOUND_375)  return false;

    return true;
}

bool SRbin5()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W > MTW2_CUT)  return false;
    if (myEvent.ETmiss <= ETMISS_BOUND_250)  return false;

    return true;
}

bool SRbin6()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ETmiss < ETMISS_BOUND_250)  return false;
    if (myEvent.ETmiss > ETMISS_BOUND_300)  return false;

    return true;
}

bool SRbin7()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ETmiss <= ETMISS_BOUND_300)  return false;
    if (myEvent.ETmiss > ETMISS_BOUND_400)  return false;

    return true;
}

bool SRbin8()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.Mlb <= MLB_CUT)  return false;
    if (myEvent.MT2W <= MTW2_CUT)  return false;
    if (myEvent.ETmiss <= ETMISS_BOUND_400)  return false;

    return true;
}

bool goesInPreselectionNoBVeto()
{
    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    //if ((!myEvent.isolatedTrackVeto) || (!myEvent.tauVeto)) return false;

    return true;
}


bool goesInPreselectionMTtailNoBeto(){return (goesInPreselectionNoBVeto() && goesInMTtail());}
bool goesInPreselectionMTtail()     { return (goesInPreselection() && goesInMTtail());     }
bool goesInPreselectionMTpeak()     { return (goesInPreselection() && goesInMTpeak());     }
bool goesInPreselectionMTinverted() { return (goesInPreselection() && goesInMTinverted()); }

bool goesIn0BtagControlRegion()
{

    if (myEvent.pfmet < MET_CUT) return false;
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags != 0 )  return false;
    //if ((!myEvent.isolatedTrackVeto) || (!myEvent.tauVeto)) return false;

    return true;
}

bool goesIn0BtagControlRegionMTtail()     { return (goesIn0BtagControlRegion() && goesInMTtail());     }
bool goesIn0BtagControlRegionMTpeak()     { return (goesIn0BtagControlRegion() && goesInMTpeak());     }
bool goesIn0BtagControlRegionMTinverted() { return (goesIn0BtagControlRegion() && goesInMTinverted()); }

bool goesInDileptonControlRegion(short int nJetCut = -1)
{
    if (myEvent.ngoodleps != 2) return false;

    if (nJetCut == -1) { if (myEvent.ngoodjets < 1) return false; }
    else if (nJetCut == 2)  { if ((myEvent.ngoodjets != 1) && (myEvent.ngoodjets != 2)) return false; }
    else if (nJetCut == 3)  { if (myEvent.ngoodjets != 3) return false; }
    else if (nJetCut == 4)  { if (myEvent.ngoodjets <  4) return false; }

    if (myEvent.ngoodbtags < NBJET_CUT)   return false;
    if (myEvent.pfmet   < MET_CUTLL)   return false;

    // Remove same-sign events
    if ((myEvent.lep1_pdgid < 0) && (myEvent.lep2_pdgid < 0)) return false;
    if ((myEvent.lep1_pdgid > 0) && (myEvent.lep2_pdgid > 0)) return false;

    // Remove Z mass peak
    if (fabs((myEvent.leadingLepton + myEvent.secondLepton).M() - 91) < 15) return false;

    return true;
}

bool goesInVetosControlRegion()
{
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    if (myEvent.ngoodjets < NJET_CUT)  return false;
    if (myEvent.ngoodbtags < NBJET_CUT)  return false;
    if (myEvent.pfmet   < MET_CUT) return false;

    // Apply reversed vetos
    //if ((myEvent.isolatedTrackVeto) && (myEvent.tauVeto)) return false;

    return true;
}

bool goesInDileptonControlRegionMTtail()     { return (goesInDileptonControlRegion() && goesInMTtail());     }
bool goesInDileptonControlRegionMTpeak()     { return (goesInDileptonControlRegion() && goesInMTpeak());     }
bool goesInDileptonControlRegionMTinverted() { return (goesInDileptonControlRegion() && goesInMTinverted()); }

bool goesInVetoControlRegionMTtail()     { return (goesInVetosControlRegion() && goesInMTtail());     }
bool goesInVetoControlRegionMTpeak()     { return (goesInVetosControlRegion() && goesInMTpeak());     }
bool goesInVetoControlRegionMTinverted() { return (goesInVetosControlRegion() && goesInMTinverted()); }

// Single-lepton channels definitions
// ##################################

bool goesInSingleElecChannel()
{
    // Keep only events with ngoodleps == 1
    if (myEvent.ngoodleps != NLEP_CUT) return false;
    // For data, keep only events from SingleElec dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "SingleElec") || (!myEvent.triggerElec)) return false;
    }*/

    // Remove electrons with pT < 30 GeV
    //if (myEvent.leadingLepton.Pt() < 30)  return false;

    // Keep only events with an electron as leading lepton
    return (abs(myEvent.lep1_pdgid) == 11);
}


bool goesInSingleMuonChannel()
{
    // Keep only events with ngoodleps == 1
    if (myEvent.ngoodleps != 1) return false;
    // For data, keep only events from SingleMuon dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "SingleMuon") || ((!myEvent.triggerMuon) && (!myEvent.xtriggerMuon))) return false;

        // Take care of the splitting due to x-trigger
        if ((myEvent.leadingLepton.Pt() >= 26) && (!myEvent.triggerMuon))  return false;
        if ((myEvent.leadingLepton.Pt() <  26) && (!myEvent.xtriggerMuon)) return false;
    }
    */

    // Keep only events with a muon as leading lepton
    return (abs(myEvent.lep1_pdgid) == 13);
}


bool goesInSingleLeptonChannel()
{
    // Single-lepton channel is the union of e-channel + mu-channel
    bool test = (goesInSingleElecChannel()
              || goesInSingleMuonChannel());
    return test;
}



// Double-lepton channels definitions
// ##################################

bool goesInDoubleElecChannel()
{
    // Keep only events with ngoodleps == 2
    if (myEvent.ngoodleps != 2) return false;
    // For data, keep only events from DoubleElec dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "DoubleElec") || (!myEvent.triggerDoubleElec)) return false;
    }
    */
    // Keep only events with two electrons
    return ((abs(myEvent.lep1_pdgid) == 11)
         && (abs(myEvent.lep2_pdgid)  == 11));
}

bool goesInDoubleMuonChannel()
{
    // Keep only events with ngoodleps == 2
    if (myEvent.ngoodleps != 2) return false;
    // For data, keep only events from DoubleMuon dataset that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "DoubleMuon") || (!myEvent.triggerDoubleMuon)) return false;
    }
    */
    // Keep only events with two muons
    return ((abs(myEvent.lep1_pdgid) == 13)
         && (abs(myEvent.lep2_pdgid)  == 13));
}

bool goesInMuonElecChannel()
{
    // Keep only events with ngoodleps == 2
    if (myEvent.ngoodleps != 2) return false;
    // For data, keep only events from SingleMuon channel that fired the trigger
    /*
    if (sampleType == "data")
    {
        if ((sampleName != "MuEl") || (!myEvent.triggerMuonElec)) return false;
    }
    */
    // Keep only events with an electron and a muon
    return   (((abs(myEvent.lep1_pdgid) == 13)
            && (abs(myEvent.lep2_pdgid)  == 11))
       ||     ((abs(myEvent.lep1_pdgid) == 11)
            && (abs(myEvent.lep2_pdgid)  == 13)));
}



bool goesInDoubleLeptonChannel()
{
    // Double-lepton channel is the union of ee, mumu and emu channels
    bool test = (goesInDoubleElecChannel()
              || goesInDoubleMuonChannel()
              ||   goesInMuonElecChannel()  );
    return test;
}



// Weighting information
// #####################

float getLumi()
{
         if (goesInSingleElecChannel())  return 19508.0;
    else if (goesInSingleMuonChannel())  return 19514.0;
    else if (goesInMuonElecChannel  ())  return 19513.0;
    else if (goesInDoubleMuonChannel())  return 19504.0;
    else if (goesInDoubleElecChannel())  return 19517.0;
    else                                 return 0.0;
}

float getWeight(const string& dataset, int nofEventsInTree, float sf_fracEvent = 1) //sf_fracEvent = used if we read only a part of the events
{
    float weight = 1.0;
    if (sampleType == "data") return weight;

    // Get the lumi
    //float lumi = getLumi();
    float lumi = 10000.;

    // Normalize to cross section times lumi
    int nofInitEvent = myEvent.totalNumberOfInitialEvent;
    /*
    cout<<"nofIniEvent = "<<nofInitEvent<<endl;
    cout<<"xs = "<<CrossSection(dataset)<<endl;
    cout<<"lumi = "<<lumi<<" sf = "<<sf_fracEvent<<endl;
    */
    //nofEventsInTree can be different from nofIniEvent due to skimming and preselection
    if(nofInitEvent!=0)
    	//weight = CrossSection(dataset) * lumi * sf_fracEvent *  nofEventsInTree/ nofInitEvent;
    	weight = (float) CrossSection(dataset) * lumi / sf_fracEvent / nofInitEvent;
    else 
    	weight = 0;
    
    // Normalize to cross section times lumi
    //weight *= myEvent.weightCrossSection * lumi;

    /*
    // Apply trigger efficiency weights for singleLepton channels
    if (myEvent.ngoodleps == 1)
    {
        weight *= myEvent.weightTriggerEfficiency;
    }
    // Apply trigger efficiency weights for doubleLepton channels
    else if (myEvent.ngoodleps == 2)
    {
        if (goesInDoubleElecChannel()) weight *= 0.96;
        if (goesInDoubleMuonChannel()) weight *= 0.88;
        if (goesInMuonElecChannel())   weight *= 0.93;
    }

    // Apply lepton ID efficiency and isolation scale factor to singleLepton channels
    if (myEvent.ngoodleps == 1)
    {
        weight *= myEvent.leadingLeptonIdEfficiency * myEvent.leadingLeptonIsoScaleFactor;

    }
    // TODO not sure about this, to be confirmed
    else if (myEvent.ngoodleps == 2)
    {
        weight *= myEvent.leadingLeptonIdEfficiency * myEvent.leadingLeptonIsoScaleFactor;
        weight *= myEvent.secondLeptonIdEfficiency  * myEvent.secondLeptonIsoScaleFactor;
    }

    // Apply pile-up weight
    // TODO : Do we confirm we'll use also this for signal ?
        weight *= myEvent.weightPileUp;


    if (sampleType == "signal")
    {
        // For signal, apply ISR reweighting
        weight *= myEvent.weightISRmodeling;

        // Check if event has been used in BDT training
        if (myEvent.isUsedInBDT == 0) weight *= 2;
        else                          weight  = 0;
    }


    // For ttbar only, apply topPt reweighting
    if (sampleName.find("ttbar") != string::npos)
        weight *= myEvent.weightTopPt;

    */
    return weight;
}


