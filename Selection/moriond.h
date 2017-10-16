#include "../common/Reader_CommonFormat_CommonBabies.h"

//global variable
babyEvent myEvent;



bool OneLep(){ return (myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto); }
//bool TwoLep(){ return ( (myEvent.ngoodleps+myEvent.nvetoleps)==2  && myEvent.PassTauVeto); }

bool CheckTrigger(bool data, string dataset){
       if(dataset == "")
            throw std::runtime_error("dataset for using the trigger info not specified");
        
        bool trigged = false;

        if(data && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 )
        {
            cout << "1 lep data event " << endl;
            if( dataset.find("data_single_electron")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl)
                trigged = true;
            else if( dataset.find("data_single_muon")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
                trigged = true;
            else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET110_MHT110 || myEvent.HLT_MET120_MHT120 || myEvent.HLT_MET) ) 
            {
                if(abs(myEvent.lep1_pdgid)==11 && myEvent.HLT_SingleEl ) //they are not doing that!!!
                    trigged = false;
                else if(abs(myEvent.lep1_pdgid)==13 && myEvent.HLT_SingleMu)
                    trigged = false;
                else
                    trigged = true;
            }
            else
            {
                cout << "for data no trigger was found 1l " << endl;
                cout << "dataset " << dataset << " myEvent.HLT_SingleEl "<< myEvent.HLT_SingleEl << " myEvent.HLT_SingleMu " <<  myEvent.HLT_SingleMu << endl;
            }
            cout << "datset " << dataset << " triggered by 1lep trigger " << trigged << endl;

        }
        else if(data && (myEvent.ngoodleps + myEvent.nvetoleps)>2 )
        {
            if( dataset.find("data_double_eg")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==11 && abs(myEvent.lep2_pdgid)==11  && myEvent.HLT_DiEl)
                trigged = true;
            else if( dataset.find("data_double_mu")!=std::string::npos &&  abs(myEvent.lep1_pdgid)==13 &&  abs(myEvent.lep2_pdgid)==13 && myEvent.HLT_DiMu)
                trigged = true;
            else if( dataset.find("data_muon_eg")!=std::string::npos &&  abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==24 && myEvent.HLT_MuE)
                trigged = true;
            else if( dataset.find("data_met")!=std::string::npos &&  (myEvent.HLT_MET110_MHT110 || myEvent.HLT_MET120_MHT120 || myEvent.HLT_MET))
            {
                if(  abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==22 && !myEvent.HLT_DiEl )
                    trigged = false;
                else if( abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==26 && !myEvent.HLT_DiMu)
                    trigged = false;
                else if(abs(myEvent.lep1_pdgid)+abs(myEvent.lep2_pdgid)==24 && !myEvent.HLT_MuE)
                    trigged = false;
                else
                    trigged = true;
            }
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

bool passFilters()
{
    bool filtRes = false;
    if(!myEvent.is_data)
        filtRes = true;
    else if(myEvent.filt_met && 
            myEvent.filt_badChargedCandidateFilter && 
            myEvent.filt_jetWithBadMuon && 
            myEvent.filt_pfovercalomet && 
            myEvent.filt_badMuonFilter && 
            !myEvent.filt_duplicatemuons && 
            !myEvent.filt_badmuons && 
            myEvent.filt_nobadmuons )
        filtRes = true;
    else
        filtRes = false ;  

    cout << "result of filters " << filtRes << endl;
    return filtRes;
        
}

bool passGoodVtx()
{
    bool selRes = false;
    if(myEvent.nvertex>=1)
        selRes= true;

    cout << "good vertex returned " << selRes << endl;
    return selRes;
}

bool tightCSVV2()
{
    if(myEvent.ngoodbtags ==0 )
        return true;
    /*for (uint32_t j =0 ; j<myEvent.ak4pfjets_CSV->size(); j++)
    {
        if( myEvent.ak4pfjets_CSV->at(j) >0.935)
            return true;
    }
    return false;*/
    if(myEvent.ntightbtags> 0)
        return true;

    return false;
}


float Mlb()
{

    if(myEvent.ngoodbtags==0)
    {
        cout << "csv size " << myEvent.ak4pfjets_CSV->size()  << endl;
        if(myEvent.ak4pfjets_CSV->size()==0)
            return -13;
        if(myEvent.ngoodleps==0)
            return -13;

 
         ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >   lb =  myEvent.ak4pfjets_leadbtag_p4 + myEvent.lep1_p4;
         float Mlb2 = lb.M2();
         if(Mlb2<0)
             throw std::runtime_error("Square root of Mlb is smaller than zero");
         return sqrt(Mlb2);
    }
    else
        return myEvent.Mlb;
}

bool dilepSel()
{

    if ( ( (myEvent.ngoodleps + myEvent.nvetoleps)>2) &&
         ( (abs(myEvent.lep1_pdgid)==13 && myEvent.lep1_passMediumID) || 
	   (abs(myEvent.lep1_pdgid)==11 && fabs(myEvent.lep1_eta)<1.4442 && myEvent.lep1_passMediumID ) ) &&
         ( (abs(myEvent.lep2_pdgid)==13 && myEvent.lep2_passMediumID) || 
	   (abs(myEvent.lep2_pdgid)==11 && fabs(myEvent.lep2_eta)<1.4442 && myEvent.lep2_passMediumID ) ) )
    {
        myEvent.pfmet = myEvent.pfmet_rl;
        myEvent.pfmet_phi = myEvent.pfmet_phi_rl;
        myEvent.MT2W = myEvent.MT2W_rl;
        myEvent.dphi_ak4pfjets_met = myEvent.dphi_ak4pfjets_met_rl;
        myEvent.mt_met_lep = myEvent.mt_met_lep_rl;
        myEvent.topnessMod = myEvent.topnessMod_rl;
        myEvent.lep1_dphiMET = myEvent.lep1_dphiMET_rl;
        myEvent.lep2_dphiMET = myEvent.lep2_dphiMET_rl;
        return true;
    }

    else
        return false;   
}

bool baseline(){ return myEvent.pfmet>=250 && myEvent.mt_met_lep>=150 && myEvent.ngoodjets>=2 && myEvent.dphi_ak4pfjets_met>=0.5 && myEvent.trigger && myEvent.topnessMod>-1000 && myEvent.Mlb>=0 && passGoodVtx() && passFilters();}
bool SR1l() { return (myEvent.ngoodbtags>=1 && myEvent.ngoodleps==1 && myEvent.nvetoleps==1 && myEvent.PassTrackVeto && myEvent.PassTauVeto && baseline()  ); }
bool SR1l_A_250lessMETless350() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_A_350lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_A_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_A_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_B_250lessMETless450() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_B_450lessMETless600() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_B_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets<=3 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=600);}
bool SR1l_C_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_C_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_C_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_C_550lessMETless650() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550 && myEvent.pfmet<650);}
bool SR1l_C_650lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=650);}
bool SR1l_D_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_D_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_D_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_D_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod<=0 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=550);}
bool SR1l_E_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_E_350lessMETless550() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<550);}
bool SR1l_E_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=550);}
bool SR1l_F_250lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_F_450lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>0  && myEvent.topnessMod<=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_G_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_G_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_G_450lessMETless600() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=450 && myEvent.pfmet<600);}
bool SR1l_G_600lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb<=175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.pfmet>=600);}
bool SR1l_H_250lessMETless450() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=250 && myEvent.pfmet<450);}
bool SR1l_H_450lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=4 && myEvent.topnessMod>=10 && myEvent.Mlb>175 && myEvent.dphi_ak4pfjets_met>=0.8 && myEvent.ntightbtags>=1 && myEvent.pfmet>=450);}
bool SR1l_I_250lessMETless350() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=250 && myEvent.pfmet<350);}
bool SR1l_I_350lessMETless450() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=350 && myEvent.pfmet<450);}
bool SR1l_I_450lessMETless550() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=450 && myEvent.pfmet<550);}
bool SR1l_I_550lessMETlessInf() { return (SR1l() && myEvent.ngoodjets>=5 && myEvent.lep1_pt<150 && myEvent.lep1_dphiMET<2 && myEvent.ak4pfjets_passMEDbtag->at(0) == false && myEvent.pfmet>=550);}
