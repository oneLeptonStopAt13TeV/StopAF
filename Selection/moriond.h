#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
//bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

bool CheckTrigger(bool data, string dataset){
       if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
        
        bool trigged = false;

        if(data)
        {
            if( dataset.find("data_single_electron")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl)
                trigged = true;
            else if( dataset.find("data_single_muon")!=std::string::npos &&  fabs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
                trigged = true;
            else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET100_MHT100 || myEvent.HLT_MET) && ( (fabs(myEvent.lep1_pdgid)==11 && !myEvent.HLT_SingleEl ) ||  (fabs(myEvent.lep1_pdgid)==13 && !myEvent.HLT_SingleMu) ))
                trigged = true;
            else
            {
            }

        }
        else
        {
            trigged = true;
        }
        return trigged;

}

bool tightCSVV2()
{
    for (uint32_t j =0 ; j<myEvent.ak4pfjets_CSV->size(); j++)
    {
        if( myEvent.ak4pfjets_CSV->at(j) >0.935)
            return true;
    }
    return false;
}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0;}
bool SR1l() { return ( baseline() && myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && ( (myEvent.lep1_passMediumID && abs(myEvent.lep1_pdgid)==11) || (myEvent.lep1_passTightID && abs(myEvent.lep1_pdgid)==13) ) ); }
bool SR1l_A_250lessMETless350() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && tightCSVV2() && myEvent.pfmet>=450);}
