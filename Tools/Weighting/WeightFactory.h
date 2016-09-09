#ifndef _WEIGHT_FACTORY_
#define _WEIGHT_FACTORY_

#include "../../sonicScrewdriver/test/common.h"
#include "BTagCalibrationStandalone.h"

#include "TFile.h"
#include "TH2D.h"

#include <vector>
#include <string>
#include <iostream>

using namespace std;

#define BTAG_DISCRI_CUT 0.80

//---------------------------
// To Be DONE
// -btagging weights: done as US group
//   does not treat SF for C ...


class WeightFactory{


	public:
	  WeightFactory();
	  ~WeightFactory();

	//-- Initalialization

	void Initialize();
	void InitializeBtagSFTool(string CSVfileFullSim, string CSVfileFastSim, string BtagEffFullSimFilename, string BtagEffFastSimFilename);
	//many things are hard-coded in that function such as the WPoint etc ...

	//-- Computation
	float GetBtagEff(const float& pt, const float& eta, const int& mcFlavour); // called by BtagWeighComputor
	void BtagWeighComputor (const vector<float>& jets_pt, const vector<float>& jets_eta, const vector<int>& jets_hadronFlavour, const vector<float>& jets_CSV);// should be called once per event

	//-- Accessors
	double GetGlobalW() {return Wglobal;}
	double GetMCW(){return Wmc;}
	double GetBtagW(){return Wbtag;}
	double GetLepW(){return Wlep;}
	
	//-- Mutators
	//-- should be called for each dataset
	void SetIsFastSim(bool isFS) {isFastSim = isFS;}
	void SetIsData(bool isD) {isData = isD;}

	private:
	  
	  //---- Global evnt info ----//
	  bool isFastSim;
	  bool isData;

	  //---- Global Weights ----//
	  double Wglobal;
	  double Wmc;
	  double Wbtag;
	  double Wlep;

	  //---- break down of the weights ----//


	  //---- Btag weighting tools 	----//

		//-- for fullsim --//
	  BTagCalibration* calib;
	  BTagCalibrationReader* reader_heavy;
	  BTagCalibrationReader* reader_heavy_UP;
	  BTagCalibrationReader* reader_heavy_DN;
	  BTagCalibrationReader* reader_light;
	  BTagCalibrationReader* reader_light_UP;
	  BTagCalibrationReader* reader_light_DN;
		
		//-- fullsim plots
	  TH2D* h_btag_eff_b_fullSim;
	  TH2D* h_btag_eff_c_fullSim;
	  TH2D* h_btag_eff_udsg_fullSim;
	  
		//-- for fastsim --//
	  
	  BTagCalibration* calib_fastsim;
	  BTagCalibrationReader* reader_fastsim;
	  BTagCalibrationReader* reader_fastsim_UP;
	  BTagCalibrationReader* reader_fastsim_DN;

		//-- fastsim plots
	  TH2D* h_btag_eff_b_fastsim;
	  TH2D* h_btag_eff_c_fastsim;
	  TH2D* h_btag_eff_udsg_fastsim;;


};

#endif
