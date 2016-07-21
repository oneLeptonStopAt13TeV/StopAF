

//------------------------------------
// Struct use to categorize the event
//------------------------------------


struct EventCategory{
	// One lepton events (CR or CR1l) - includes the vetoes
	bool OneLep;
	bool OneLep_e;
	bool OneLep_mu;

	// Two lepton events (CR-2l)
	bool TwoLep;
	bool TwoLep_ee;
	bool TwoLep_emu;
	bool TwoLep_mumu;

	// At least one b-tagged jets
	bool OneBtag;

	//Baseline: >=2 jets, MET, MT, DPhi
	bool PassBaselineCuts;
	
	// 2 jets bins (include mod_top>6.4
	bool TwoJets_MET_250;
	bool TwoJets_MET_350;
	bool TwoJets_MET_450;

	// 3 jets bins
	bool FourJets_hightMT2W_MET_250;
	bool FourJets_hightMT2W_MET_350;
	bool FourJets_hightMT2W_MET_450;
	bool FourJets_hightMT2W_MET_550;

	// >=4 jets bins
	
	// low MT2W
	bool FourJets_lowMT2W_MET_250;
	bool FourJets_lowMT2W_MET_350;
	bool FourJets_lowMT2W_MET_450;

	//- high MT2W
	bool FourJets_highMT2W_MET_250;
	bool FourJets_highMT2W_MET_350;
	bool FourJets_highMT2W_MET_450;
	bool FourJets_highMT2W_MET_550;
	bool FourJets_highMT2W_MET_650;


	void Reset();
	void Update();
};

EventCategory eventCat;

bool EventCategory::Reset(){
	// One lepton events (CR or CR1l) - includes the vetoes
	 OneLep = false;
	 OneLep_e = false;
	 OneLep_mu = false;
	// Two lepton events (CR-2l)
	 TwoLep = false;
	 TwoLep_ee = false;
	 TwoLep_emu = false;
	 TwoLep_mumu = false;
	// At least one b-tagged jets
	 OneBtag = false;
	//Baseline: >=2 jets, MET, MT, DPhi
	 PassBaselineCuts = false;
	// 2 jets bins (include mod_top>6.4
	 TwoJets_MET_250 = false;
	 TwoJets_MET_350 = false;
	 TwoJets_MET_450 = false;
	// 3 jets bins
	 FourJets_hightMT2W_MET_250 = false;
	 FourJets_hightMT2W_MET_350 = false;
	 FourJets_hightMT2W_MET_450 = false;
	 FourJets_hightMT2W_MET_550 = false;
	// >=4 jets bins
	// low MT2W
	 FourJets_lowMT2W_MET_250 = false;
	 FourJets_lowMT2W_MET_350 = false;
	 FourJets_lowMT2W_MET_450 = false;
	//- high MT2W
	 FourJets_highMT2W_MET_250 = false;
	 FourJets_highMT2W_MET_350 = false;
	 FourJets_highMT2W_MET_450 = false;
	 FourJets_highMT2W_MET_550 = false;
	 FourJets_highMT2W_MET_650 = false;
}

bool EventCategory::Update(){
	// One lepton events (CR or CR1l) - includes the vetoes
	if(myEvent.ngoodleps == 1 && myEvent.PassTrackVeto && myEvent.PassTauVeto){
		OneLep = true;
		if(abs(myEvent.lep1_pdgid) == 11) OneLep_e = true;
		if(abs(myEvent.lep1_pdgid) == 13) OneLep_mu = true;
	}

	if( (myEvent.ngoodleps+nvetoleps)==2 ){
		TwoLep = true;
		if(abs(myEvent.lep1_pdgid) == 11 && abs(myEvent.lep2_pdgid) == 11) TwoLep_ee = true;
		if(abs(myEvent.lep1_pdgid) == 11 && abs(myEvent.lep2_pdgid) == 13) TwoLep_emu = true;
		if(abs(myEvent.lep1_pdgid) == 13 && abs(myEvent.lep2_pdgid) == 11) TwoLep_emu = true;
		if(abs(myEvent.lep1_pdgid) == 13 && abs(myEvent.lep2_pdgid) == 13) TwoLep_mumu = true;
	}

	// At least one b-tagged jets
	if( myEvent.ngoodbtags>=1 ) OneBtag = true; 

	//Baseline: >=2 jets, MET, MT, DPhi
    	if (myEvent.pfmet >= 250 && myEvent.mt_met_lep > 150 && myEvent.ngoodjets >=2 myEvent.dphi_ak4pfjets_met >= 0.8 )   PassBaselineCuts = true;
}
	
	// 2 jets bins (include mod_top>6.4
	if(myEvent.ngoodjets == 2 && myEvent.topness>6.4){
		if(myEvent.pfmet >= 250 && myEvent.pfmet <350) TwoJets_MET_250 = true;
		if(myEvent.pfmet >= 350 && myEvent.pfmet <450) TwoJets_MET_350 = true;
		if(myEvent.pfmet >= 450) TwoJets_MET_450 = true;
	}

	// 3 jets bins
	if(myEvent.ngoodjets == 3){
		if(myEvent.pfmet >= 250 && myEvent.pfmet <350) FourJets_hightMT2W_MET_250 = true;
		if(myEvent.pfmet >= 350 && myEvent.pfmet <450) FourJets_hightMT2W_MET_350 = true;
		if(myEvent.pfmet >= 450 && myEvent.pfmet <550) FourJets_hightMT2W_MET_450 = true;
		if(myEvent.pfmet >= 550) FourJets_hightMT2W_MET_550 = true;
	}

	// >=4 jets bins
	if(myEvent.ngoodjets == 4){
		if(myEvent.MT2W<200){
			if(myEvent.pfmet >= 250 && myEvent.pfmet <350) FourJets_lowMT2W_MET_250 = true;
			if(myEvent.pfmet >= 350 && myEvent.pfmet <450) FourJets_lowMT2W_MET_350 = true;
			if(myEvent.pfmet >= 450) FourJets_lowMT2W_MET_450 = true;
		}
		else{
			if(myEvent.pfmet >= 250 && myEvent.pfmet <350) FourJets_highMT2W_MET_250 = true;
			if(myEvent.pfmet >= 350 && myEvent.pfmet <450) FourJets_highMT2W_MET_350 = true;
			if(myEvent.pfmet >= 450 && myEvent.pfmet <550) FourJets_highMT2W_MET_450 = true;
			if(myEvent.pfmet >= 550 && myEvent.pfmet <550) FourJets_highMT2W_MET_450 = true;
			if(myEvent.pfmet >= 650) FourJets_highMT2W_MET_450 = true;
		}
	}
	
}


bool 

// -- Function used for CR2l events where 2nd lepton is added to the MET
void RecomputeMETAndDerivedQuantities(){
}

//--------------
//   Channels 
//--------------

//-- Signal regions
bool SR1l() {return eventCat.OneLep && eventCat.PassBaselineCuts && eventCat.OneBtag;}
bool SR1l_e() {return eventCat.OneLep_e && eventCat.PassBaselineCuts && eventCat.OneBtag;}
bool SR1l_mu() {return eventCat.OneLep_mu && eventCat.PassBaselineCuts && eventCat.OneBtag;}

//-- CR1l regions
bool CR1l() {return eventCat.OneLep && eventCat.PassBaselineCuts && !eventCat.OneBtag;}
bool CR1l_e() {return eventCat.OneLep_e && eventCat.PassBaselineCuts && !eventCat.OneBtag;}
bool CR1l_mu() {return eventCat.OneLep_mu && eventCat.PassBaselineCuts && !eventCat.OneBtag;}


// Should be add Z-veto and b-tagging requirement ??
bool CR2l() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};
bool CR2l_ee() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};
bool CR2l_emu() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};
bool CR2l_mumu() {return eventCat.TwoLep && eventCat.PassBaselineCuts;};

//-----------------
// Selection  Bins
//-----------------


// 2 jets bins (include mod_top>6.4
bool TwoJets_MET_250(){ return eventCat.TwoJets_MET_350;}
bool TwoJets_MET_350(){ return eventCat.TwoJets_MET_350;}
bool TwoJets_MET_450(){ return eventCat.TwoJets_MET_450;}

// 3 jets bins
bool FourJets_hightMT2W_MET_250(){ return eventCat.FourJets_hightMT2W_MET_250;}
bool FourJets_hightMT2W_MET_350(){ return eventCat.FourJets_hightMT2W_MET_350;}
bool FourJets_hightMT2W_MET_450(){ return eventCat.FourJets_hightMT2W_MET_450;}
bool FourJets_hightMT2W_MET_550(){ return eventCat.FourJets_hightMT2W_MET_550;}

// >=4 jets bins

// low MT2W
bool FourJets_lowMT2W_MET_250(){ return eventCat.FourJets_lowMT2W_MET_250;}
bool FourJets_lowMT2W_MET_350(){ return eventCat.FourJets_lowMT2W_MET_350;}
bool FourJets_lowMT2W_MET_450(){ return eventCat.FourJets_lowMT2W_MET_450;}

//- high MT2W
bool FourJets_highMT2W_MET_250(){ return eventCat.FourJets_highMT2W_MET_250;}
bool FourJets_highMT2W_MET_350(){ return eventCat.FourJets_highMT2W_MET_350;}
bool FourJets_highMT2W_MET_450(){ return eventCat.FourJets_highMT2W_MET_450;}
bool FourJets_highMT2W_MET_550(){ return eventCat.FourJets_highMT2W_MET_550;}
bool FourJets_highMT2W_MET_650(){ return eventCat.FourJets_highMT2W_MET_650;}


