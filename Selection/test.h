#include "Reader_CommonFormat.h"

//global variable
babyEvent myEvent;
bool trigged;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

void CheckTrigger(bool data, string dataset){
        if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");

	bool passMETMHTTrigger = false;
        bool passElTrigger = false;
        bool passMuTrigger = false;
	
	trigged = false

        if(data)
        {
		// METMH, electron or muon trigger 
		//for(unsigned int i=0; i<myEvent.trigger_name.size();i++)
		//{
			// -- MET MHT trigger
			if(myEvent.HLT_PFMET170 == true && dataset.find("MET")!=std::string::npos)
                                passMETMHTTrigger = true;
			if(myEvent.HLT_PFMET100_PFMHT100_IDTight  == true && dataset.find("MET")!=std::string::npos)
				passMETMHTTrigger = true;
			// -- Electron trigger
			if(myEvent.HLT_Ele25_eta2p1_WPLoose  == true && dataset.find("SE")!=std::string::npos)
			        if(abs(myEvent.lep1_pdgid) == 11) passElTrigger = true;
			if(myEvent.HLT_Ele27_eta2p1_WPLoose  == true && dataset.find("SE")!=std::string::npos )		
				if(abs(myEvent.lep1_pdgid) == 11) passElTrigger = true;
			
			// -- Muon trigger
			if(myEvent.HLT_IsoMu20 == true && dataset.find("SM")!=std::string::npos)
				if(abs(myEvent.lep1_pdgid) == 13) passMuTrigger = true;
			if(myEvent.HLT_IsoMu22  == true && dataset.find("SM")!=std::string::npos)
				if(abs(myEvent.lep1_pdgid) == 13) passMuTrigger = true; //string in recompute that I have the right dataset name
		 //}

		 if(passMETMHTTrigger == true || passElTrigger == true || passMuTrigger == true) trigged = true;
         }
         else
         {
             trigged = true;
         }
}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && DPhi_jmyEvent.pfmet>=0.8 && trigged();}
bool SR1l() { return ( baseline() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto ); }
bool CR1l() { return ( baseline() && myEvent.ngoodbtags==0 && myEvent.ngoodleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto ); }
bool CR2l() { return ( baseline() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==2 && myEvent.PassTrackVeto && myEvent.PassTauVeto ); }
bool SR1l23jLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l23jLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l23jLowMlb_MET450to550() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l23jLowMlb_MET550toInf() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=550);}
bool SR1l23jHighMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l23jHighMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l23jHighMlb_MET450to550() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l23jHighMlb_MET550toInf() { return (SR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=550);}
bool SR1l4jLowMT2WLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jLowMT2WLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jLowMT2WLowMlb_MET450to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l4jLowMT2WLowMlb_MET550to650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l4jLowMT2WLowMlb_MET650toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=650);}
bool SR1l4jLowMT2WHighMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jLowMT2WHighMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jLowMT2WHighMlb_MET450to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l4jLowMT2WHighMlb_MET550toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=550);}
bool SR1l4jMidMT2WLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jMidMT2WLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jMidMT2WLowMlb_MET450toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450);}
bool SR1l4jMidMT2WHighMlb_MET250to400() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<400);}
bool SR1l4jMidMT2WHighMlb_MET400toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=400);}
bool SR1l4jHighMT2WLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jHighMT2WLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jHighMT2WLowMlb_MET450to600() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l4jHighMT2WLowMlb_MET600toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=600);}
bool SR1l4jHighMT2WHighMlb_MET250to400() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<400);}
bool SR1l4jHighMT2WHighMlb_MET400to650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=400 && myEvent.pfmet<650);}
bool SR1l4jHighMT2WHighMlb_MET650toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=650);}
bool CR1l23jLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l23jLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l23jLowMlb_MET450to550() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l23jLowMlb_MET550toInf() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=550);}
bool CR1l23jHighMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l23jHighMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l23jHighMlb_MET450to550() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l23jHighMlb_MET550toInf() { return (CR1l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=550);}
bool CR1l4jLowMT2WLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jLowMT2WLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jLowMT2WLowMlb_MET450to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l4jLowMT2WLowMlb_MET550to650() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR1l4jLowMT2WLowMlb_MET650toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=650);}
bool CR1l4jLowMT2WHighMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jLowMT2WHighMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jLowMT2WHighMlb_MET450to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l4jLowMT2WHighMlb_MET550toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=550);}
bool CR1l4jMidMT2WLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jMidMT2WLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jMidMT2WLowMlb_MET450toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450);}
bool CR1l4jMidMT2WHighMlb_MET250to400() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<400);}
bool CR1l4jMidMT2WHighMlb_MET400toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=400);}
bool CR1l4jHighMT2WLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jHighMT2WLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jHighMT2WLowMlb_MET450to600() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR1l4jHighMT2WLowMlb_MET600toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=600);}
bool CR1l4jHighMT2WHighMlb_MET250to400() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<400);}
bool CR1l4jHighMT2WHighMlb_MET400to650() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=400 && myEvent.pfmet<650);}
bool CR1l4jHighMT2WHighMlb_MET650toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=650);}
bool CR2l23jLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l23jLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l23jLowMlb_MET450to550() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l23jLowMlb_MET550toInf() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb<175 && myEvent.pfmet>=550);}
bool CR2l23jHighMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l23jHighMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l23jHighMlb_MET450to550() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l23jHighMlb_MET550toInf() { return (CR2l() && myEvent.ngoodjets<3 && myEvent.topness>7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=550);}
bool CR2l4jLowMT2WLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jLowMT2WLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jLowMT2WLowMlb_MET450to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l4jLowMT2WLowMlb_MET550to650() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR2l4jLowMT2WLowMlb_MET650toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb<175 && myEvent.pfmet>=650);}
bool CR2l4jLowMT2WHighMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jLowMT2WHighMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jLowMT2WHighMlb_MET450to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l4jLowMT2WHighMlb_MET550toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<0 && myEvent.Mlb>=175 && myEvent.pfmet>=550);}
bool CR2l4jMidMT2WLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jMidMT2WLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jMidMT2WLowMlb_MET450toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450);}
bool CR2l4jMidMT2WHighMlb_MET250to400() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<400);}
bool CR2l4jMidMT2WHighMlb_MET400toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=400);}
bool CR2l4jHighMT2WLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jHighMT2WLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jHighMT2WLowMlb_MET450to600() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR2l4jHighMT2WLowMlb_MET600toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb<175 && myEvent.pfmet>=600);}
bool CR2l4jHighMT2WHighMlb_MET250to400() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=250 && myEvent.pfmet<400);}
bool CR2l4jHighMT2WHighMlb_MET400to650() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=400 && myEvent.pfmet<650);}
bool CR2l4jHighMT2WHighMlb_MET650toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=7.5 && myEvent.Mlb>=175 && myEvent.pfmet>=650);}
