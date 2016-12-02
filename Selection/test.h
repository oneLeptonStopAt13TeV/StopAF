#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;
bool trigged;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

void CheckTrigger(bool data, string dataset){
        /*if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");

	bool passMETMHTTrigger = false;
        bool passElTrigger = false;
        bool passMuTrigger = false;
	
	trigged = false;

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
         }*/
}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.8;}
bool SR1l() { return ( baseline() && myEvent.ngoodbtags>=1 && OneLep() && ( (myEvent.lep1_passMediumID && abs(myEvent.lep1_pdgid)==11) || (myEvent.lep1_passTightID && abs(myEvent.lep1_pdgid)==13))  ); }
bool CR1l() { return ( baseline() && myEvent.ngoodbtags==0 && OneLep() ); }
bool CR2l() { return ( baseline() && myEvent.ngoodbtags>=1 && TwoLep() ); }
bool SR1l2j_MET250to350() { return (SR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l2j_MET350to450() { return (SR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l2j_MET450toInf() { return (SR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=450);}
bool SR1l3j_MET250to350() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l3j_MET350to450() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l3j_MET450to550() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l3j_MET550toInf() { return (SR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=550);}
bool SR1l4jLow_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jLow_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jLow_MET450toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=450);}
bool SR1l4jHighMT2W_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jHighMT2W_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jHighMT2W_MET450to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l4jHighMT2W_MET550to650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l4jHighMT2W_MET650toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=650);}
bool CR1l2j_MET250to350() { return (CR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l2j_MET350to450() { return (CR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l2j_MET450toInf() { return (CR1l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=450);}
bool CR1l3j_MET250to350() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l3j_MET350to450() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l3j_MET450to550() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l3j_MET550toInf() { return (CR1l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=550);}
bool CR1l4jLow_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jLow_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jLow_MET450toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=450);}
bool CR1l4jHighMT2W_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jHighMT2W_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jHighMT2W_MET450to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l4jHighMT2W_MET550to650() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR1l4jHighMT2W_MET650toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=650);}
bool CR2l2j_MET250to350() { return (CR2l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l2j_MET350to450() { return (CR2l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l2j_MET450toInf() { return (CR2l() && myEvent.ngoodjets==2 && myEvent.topnessMod>6.4 && myEvent.pfmet>=450);}
bool CR2l3j_MET250to350() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l3j_MET350to450() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l3j_MET450to550() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l3j_MET550toInf() { return (CR2l() && myEvent.ngoodjets==3 && myEvent.MT2W>=200 && myEvent.pfmet>=550);}
bool CR2l4jLow_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jLow_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jLow_MET450toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W<200 && myEvent.pfmet>=450);}
bool CR2l4jHighMT2W_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jHighMT2W_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jHighMT2W_MET450to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l4jHighMT2W_MET550to650() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR2l4jHighMT2W_MET650toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.MT2W>=200 && myEvent.pfmet>=650);}
