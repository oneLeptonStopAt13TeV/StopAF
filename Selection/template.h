#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

bool CheckTrigger(bool data, string dataset){
       if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
        
        bool trigged = false;

        if(data)
        {
          if( dataset.find("data_single_electron")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl)
          {
            //cout << "electron trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << endl; 
            trigged = true;
          }
          else if( dataset.find("data_single_muon")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
          {
            //cout << "muon trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << endl; 
            trigged = true;
          }
          else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET100_MHT100 || myEvent.HLT_MET) && ( (fabs(myEvent.lep1_pdgid)==11 && !myEvent.HLT_SingleEl ) ||  (fabs(myEvent.lep1_pdgid)==13 && !myEvent.HLT_SingleMu) ))
          {
            //cout << "HT trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << " MET trigger " << myEvent.HLT_MET << endl; 
            trigged = true;
          }
          else
          {
            //cout << "no trigger; dataset " << dataset << " lep " << myEvent.lep1_pdgid << " el trigger " << myEvent.HLT_SingleEl << " mu trigger " << myEvent.HLT_SingleMu << " HT trigger " << myEvent.HLT_MET100_MHT100 << " MET trigger " << myEvent.HLT_MET << endl; 
     
          }

         }
         else
         {
             trigged = true;
         }
         return trigged;

}

