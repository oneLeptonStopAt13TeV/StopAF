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

