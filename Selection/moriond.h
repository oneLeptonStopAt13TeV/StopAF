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
          if( dataset.find("data_single_electron")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl && !myEvent.HLT_MET100_MHT100 )
            trigged = true;
          if( dataset.find("data_single_muon")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu && !myEvent.HLT_MET100_MHT100 )
            trigged = true;
          if( dataset.find("data_met")!=std::string::npos &&  myEvent.HLT_MET100_MHT100 && ( (fabs(myEvent.lep1_pdgid)==11 && !myEvent.HLT_SingleEl ) ||  (fabs(myEvent.lep1_pdgid)==13 && !myEvent.HLT_SingleMu) ))
            trigged = true;

         }
         else
         {
             trigged = true;
         }
         return trigged;

}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.trigger;}
bool SR1l() { return ( baseline() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && ( (myEvent.lep1_passMediumID && abs(myEvent.lep1_pdgid)==11) || (myEvent.lep1_passTightID && abs(myEvent.lep1_pdgid)==13) ) ); }
bool CR1l() { return ( baseline() && myEvent.ngoodbtags==0 && myEvent.ngoodleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto ); }
bool CR2l() { return ( baseline() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==2 && myEvent.PassTrackVeto && myEvent.PassTauVeto ); }
bool CR3l2b() { return ( baseline() && myEvent.ngoodbtags>=2 && myEvent.ngoodleps==3 && myEvent.PassTrackVeto && myEvent.PassTauVeto ); }
bool SR1l23jLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l23jLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l23jLowMlb_MET450to600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l23jLowMlb_MET600toInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool SR1l23jHighMlb_MET250to450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l23jHighMlb_MET450to600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l23jHighMlb_MET600toInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=600);}
bool SR1l4jLowMT2WLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jLowMT2WLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jLowMT2WLowMlb_MET450to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l4jLowMT2WLowMlb_MET550to650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l4jLowMT2WLowMlb_MET650toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=650);}
bool SR1l4jLowMT2WHighMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jLowMT2WHighMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jLowMT2WHighMlb_MET450to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l4jLowMT2WHighMlb_MET550toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=550);}
bool SR1l4jMidMT2WLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jMidMT2WLowMlb_MET350to550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l4jMidMT2WLowMlb_MET550toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=550);}
bool SR1l4jMidMT2WHighMlb_MET250to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l4jMidMT2WHighMlb_MET450toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool SR1l4jHighMT2WLowMlb_MET250to350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l4jHighMT2WLowMlb_MET350to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l4jHighMT2WLowMlb_MET450to600() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l4jHighMT2WLowMlb_MET600toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool SR1l4jHighMT2WHighMlb_MET250to450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l4jHighMT2WHighMlb_MET450toInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool CR1l23jLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l23jLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l23jLowMlb_MET450to600() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR1l23jLowMlb_MET600toInf() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool CR1l23jHighMlb_MET250to450() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR1l23jHighMlb_MET450to600() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR1l23jHighMlb_MET600toInf() { return (CR1l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=600);}
bool CR1l4jLowMT2WLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jLowMT2WLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jLowMT2WLowMlb_MET450to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l4jLowMT2WLowMlb_MET550to650() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR1l4jLowMT2WLowMlb_MET650toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=650);}
bool CR1l4jLowMT2WHighMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jLowMT2WHighMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jLowMT2WHighMlb_MET450to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR1l4jLowMT2WHighMlb_MET550toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=550);}
bool CR1l4jMidMT2WLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jMidMT2WLowMlb_MET350to550() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool CR1l4jMidMT2WLowMlb_MET550toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=550);}
bool CR1l4jMidMT2WHighMlb_MET250to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR1l4jMidMT2WHighMlb_MET450toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool CR1l4jHighMT2WLowMlb_MET250to350() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR1l4jHighMT2WLowMlb_MET350to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR1l4jHighMT2WLowMlb_MET450to600() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR1l4jHighMT2WLowMlb_MET600toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool CR1l4jHighMT2WHighMlb_MET250to450() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR1l4jHighMT2WHighMlb_MET450toInf() { return (CR1l() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool CR2l23jLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l23jLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l23jLowMlb_MET450to600() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR2l23jLowMlb_MET600toInf() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool CR2l23jHighMlb_MET250to450() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR2l23jHighMlb_MET450to600() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR2l23jHighMlb_MET600toInf() { return (CR2l() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=600);}
bool CR2l4jLowMT2WLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jLowMT2WLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jLowMT2WLowMlb_MET450to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l4jLowMT2WLowMlb_MET550to650() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR2l4jLowMT2WLowMlb_MET650toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=650);}
bool CR2l4jLowMT2WHighMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jLowMT2WHighMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jLowMT2WHighMlb_MET450to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR2l4jLowMT2WHighMlb_MET550toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=550);}
bool CR2l4jMidMT2WLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jMidMT2WLowMlb_MET350to550() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool CR2l4jMidMT2WLowMlb_MET550toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=550);}
bool CR2l4jMidMT2WHighMlb_MET250to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR2l4jMidMT2WHighMlb_MET450toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool CR2l4jHighMT2WLowMlb_MET250to350() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR2l4jHighMT2WLowMlb_MET350to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR2l4jHighMT2WLowMlb_MET450to600() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR2l4jHighMT2WLowMlb_MET600toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool CR2l4jHighMT2WHighMlb_MET250to450() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR2l4jHighMT2WHighMlb_MET450toInf() { return (CR2l() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool CR3l2b23jLowMlb_MET250to350() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR3l2b23jLowMlb_MET350to450() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR3l2b23jLowMlb_MET450to600() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR3l2b23jLowMlb_MET600toInf() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool CR3l2b23jHighMlb_MET250to450() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR3l2b23jHighMlb_MET450to600() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR3l2b23jHighMlb_MET600toInf() { return (CR3l2b() && myEvent.ngoodjets<=3 && myEvent.topness>10 && myEvent.Mlb>175 && myEvent.pfmet>=600);}
bool CR3l2b4jLowMT2WLowMlb_MET250to350() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR3l2b4jLowMT2WLowMlb_MET350to450() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR3l2b4jLowMT2WLowMlb_MET450to550() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR3l2b4jLowMT2WLowMlb_MET550to650() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool CR3l2b4jLowMT2WLowMlb_MET650toInf() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=650);}
bool CR3l2b4jLowMT2WHighMlb_MET250to350() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR3l2b4jLowMT2WHighMlb_MET350to450() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR3l2b4jLowMT2WHighMlb_MET450to550() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool CR3l2b4jLowMT2WHighMlb_MET550toInf() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness<=0 && myEvent.Mlb>175 && myEvent.pfmet>=550);}
bool CR3l2b4jMidMT2WLowMlb_MET250to350() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR3l2b4jMidMT2WLowMlb_MET350to550() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool CR3l2b4jMidMT2WLowMlb_MET550toInf() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>0  && myEvent.topness<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=550);}
bool CR3l2b4jMidMT2WHighMlb_MET250to450() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR3l2b4jMidMT2WHighMlb_MET450toInf() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>=0  && myEvent.topness<=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
bool CR3l2b4jHighMT2WLowMlb_MET250to350() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool CR3l2b4jHighMT2WLowMlb_MET350to450() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool CR3l2b4jHighMT2WLowMlb_MET450to600() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool CR3l2b4jHighMT2WLowMlb_MET600toInf() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool CR3l2b4jHighMT2WHighMlb_MET250to450() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool CR3l2b4jHighMT2WHighMlb_MET450toInf() { return (CR3l2b() && myEvent.ngoodjets>=4 && myEvent.topness>=10 && myEvent.Mlb>175 && myEvent.pfmet>=450);}
