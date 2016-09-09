#include "WeightFactory.h"

WeightFactory::WeightFactory ()
{
}

WeightFactory::~WeightFactory ()
{
}

/*
void WeightFactory::Initialize ()
{
}
*/

//do we also need a file for Btag eff FastStim sample ?
void WeightFactory::InitializeBtagSFTool (string CSVfileFullSim, string CSVfileFastSim, string BtagEffFullSimFilename, string BtagEffFastSimFilename)
{
  //CSVfileFullSim = "btagsf/CSVv2_ichep_slimmed.csv";
  //CSVfileFastSim = "btagsf/CSV_13TEV_Combined_14_7_2016.csv"
  printBoxedMessage ("WeightFactory: initialization of BtagSFTool");
  printBoxedMessage (" > FullSim SF file: " + CSVfileFullSim);
  printBoxedMessage (" > FastSim SF file: " + CSVfileFastSim);
  printBoxedMessage (" > FullSim BtagEff file: " + BtagEffFullSimFilename);
  printBoxedMessage (" > FastSim BtagEff file: " + BtagEffFastSimFilename);

  //-- Loading SF --//
  //
  // fullSim
  cout<<"toto"<<endl;
  calib = new BTagCalibration ("csvv2", CSVfileFullSim.c_str ());
  cout<<"toto"<<endl;
  //reader_heavy = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "central");	// central
  reader_heavy = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  cout<<"toto"<<endl;
  reader_heavy->load(*calib,BTagEntry::FLAV_B,"comb");
  cout<<"toto"<<endl;
  //reader_heavy_UP = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "up");	// sys up
  reader_heavy_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_heavy_UP->load(*calib,BTagEntry::FLAV_B,"comb");
  cout<<"toto"<<endl;
  //reader_heavy_DN = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "comb", "down");	// sys down
  reader_heavy_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_heavy_DN->load(*calib,BTagEntry::FLAV_B,"comb");
  cout<<"toto"<<endl;
  //reader_light = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "incl", "central");	// central
  reader_light = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  reader_light->load(*calib,BTagEntry::FLAV_UDSG,"incl");
  cout<<"toto"<<endl;
  //reader_light_UP = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "incl", "up");	// sys up
  reader_light_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_light_UP->load(*calib,BTagEntry::FLAV_UDSG,"incl");
  cout<<"toto"<<endl;
  //reader_light_DN = new BTagCalibrationReader (calib, BTagEntry::OP_MEDIUM, "incl", "down");	// sys down
  reader_light_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_light_DN->load(*calib,BTagEntry::FLAV_UDSG,"incl");
  cout<<"toto"<<endl;
  
  // fastSim
  calib_fastsim = new BTagCalibration("CSV", CSVfileFastSim.c_str()); 
  //reader_fastsim = new BTagCalibrationReader (calib_fastsim, BTagEntry::OP_MEDIUM, "fastsim", "central");	// central
  reader_fastsim = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "central");	// central
  cout<<"toto"<<endl;
  reader_fastsim->load(*calib_fastsim,BTagEntry::FLAV_UDSG,"fastsim");
  cout<<"toto"<<endl;
  //reader_fastsim_UP = new BTagCalibrationReader (calib_fastsim, BTagEntry::OP_MEDIUM, "fastsim", "up");	// sys up
  reader_fastsim_UP = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "up");	// sys up
  reader_fastsim_UP->load(*calib_fastsim,BTagEntry::FLAV_UDSG,"fastsim");
  cout<<"toto"<<endl;
  //reader_fastsim_DN = new BTagCalibrationReader (calib_fastsim, BTagEntry::OP_MEDIUM, "fastsim", "down");	// sys down
  reader_fastsim_DN = new BTagCalibrationReader ( BTagEntry::OP_MEDIUM, "down");	// sys down
  reader_fastsim_DN->load(*calib_fastsim,BTagEntry::FLAV_UDSG,"fastsim");

  //-- Loading efficiency --//
  //-- fullSim
  TFile *file_eff = TFile::Open (BtagEffFullSimFilename.c_str ());
  h_btag_eff_b_fullSim = (TH2D *) file_eff->Get ("h2_BTaggingEff_csv_med_Eff_b")->Clone ();
  h_btag_eff_c_fullSim = (TH2D *) file_eff->Get ("h2_BTaggingEff_csv_med_Eff_c")->Clone ();
  h_btag_eff_udsg_fullSim = (TH2D *) file_eff->Get ("h2_BTaggingEff_csv_med_Eff_udsg")->Clone ();
  //-- fastSim
  TFile *file_eff_FastSim = TFile::Open (BtagEffFastSimFilename.c_str ());
  h_btag_eff_b_fastsim = (TH2D *) file_eff_FastSim->Get ("h2_BTaggingEff_csv_med_Eff_b")->Clone ();
  h_btag_eff_c_fastsim = (TH2D *) file_eff_FastSim->Get ("h2_BTaggingEff_csv_med_Eff_c")->Clone ();
  h_btag_eff_udsg_fastsim = (TH2D *) file_eff_FastSim->Get ("h2_BTaggingEff_csv_med_Eff_udsg")->Clone ();

  if (!h_btag_eff_b_fullSim || !h_btag_eff_c_fullSim || !h_btag_eff_udsg_fullSim)
    std::cerr << " WeightFactory::InitializeBtagSFTool:\t fullSim efficiency plots were not retrieved .." << endl;
  if (!h_btag_eff_b_fastsim || !h_btag_eff_c_fastsim || !h_btag_eff_udsg_fastsim)
    std::cerr << " WeightFactory::InitializeBtagSFTool:\t fastSim efficiency plots were not retrieved .." << endl;

}

float WeightFactory::GetBtagEff (const float& pt, const float& eta, const int& mcFlavour)
{
  // only use pt bins up to 400 GeV for charm and udsg
  float pt_cutoff = std::max (20., std::min (399., double (pt)));
  TH2D *h (0);
  if (abs (mcFlavour) == 5) {
    if (isFastSim)
      h = h_btag_eff_b_fastsim;
    else
      h = h_btag_eff_b_fullSim;
    // use pt bins up to 600 GeV for b
    pt_cutoff = std::max (20., std::min (599., double (pt)));
  }
  else if (abs (mcFlavour) == 4) {
    if (isFastSim)
      h = h_btag_eff_c_fastsim;
    else
      h = h_btag_eff_c_fullSim;
  }
  else {
    if (isFastSim)
      h = h_btag_eff_udsg_fastsim;
    else
      h = h_btag_eff_udsg_fullSim;
  }

  int binx = h->GetXaxis ()->FindBin (pt_cutoff);
  int biny = h->GetYaxis ()->FindBin (fabs (eta));
  return h->GetBinContent (binx, biny);
}
void WeightFactory::BtagWeighComputor (const vector<float>& jets_pt, const vector<float>& jets_eta, const vector<int>& jets_hadronFlavour, const vector<float>& jets_CSV)
{
  double btagprob_data (1.), btagprob_mc (1.), btagprob_heavy_UP (1.), btagprob_heavy_DN (1.), btagprob_light_UP (1.), btagprob_light_DN (1.), btagprob_FS_UP (1.), btagprob_FS_DN (1.);

  for (size_t i = 0; i < jets_pt.size (); ++i) {
    // Check if the jet pass the medium btag cut
    if (jets_CSV[i] > BTAG_DISCRI_CUT) {
      if (!isData) {
	float eff = GetBtagEff (jets_pt[i], jets_eta[i], jets_hadronFlavour[i]);
	BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;

	if (abs (jets_hadronFlavour[i]) == 5)
	  flavor = BTagEntry::FLAV_B;
	else if (abs (jets_hadronFlavour[i]) == 4)
	  flavor = BTagEntry::FLAV_C;

	float pt_cutoff = std::max (30., std::min (669., double (jets_pt[i])));
	float eta_cutoff = std::min (2.39, fabs (double (jets_eta[i])));
	float weight_cent (1.), weight_UP (1.), weight_DN (1.), weight_FS_UP (1.), weight_FS_DN (1.);

	if (flavor == BTagEntry::FLAV_UDSG) {
	  weight_cent = reader_light->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP = reader_light_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_light_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	else {
	  weight_cent = reader_heavy->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP = reader_heavy_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_heavy_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	if (isFastSim) {
	  weight_FS_UP = weight_cent * reader_fastsim_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_FS_DN = weight_cent * reader_fastsim_DN->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_cent *= reader_fastsim->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP *= reader_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	  weight_DN *= reader_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	}
	btagprob_data *= weight_cent * eff;
	btagprob_mc *= eff;
	if (flavor == BTagEntry::FLAV_UDSG) {
	  btagprob_light_UP *= weight_UP * eff;
	  btagprob_light_DN *= weight_DN * eff;
	  btagprob_heavy_UP *= weight_cent * eff;
	  btagprob_heavy_DN *= weight_cent * eff;
	}
	else {
	  btagprob_light_UP *= weight_cent * eff;
	  btagprob_light_DN *= weight_cent * eff;
	  btagprob_heavy_UP *= weight_UP * eff;
	  btagprob_heavy_DN *= weight_DN * eff;
	}
	if (isFastSim) {
	  btagprob_FS_UP *= weight_FS_UP * eff;
	  btagprob_FS_DN *= weight_FS_DN * eff;
	}
      }
    }
    else {			// event if it fails med btag -- it's  needed for SF event weights
      if (!isData) {
	float eff = GetBtagEff (jets_eta[i], jets_eta[i], jets_hadronFlavour[i]);
	BTagEntry::JetFlavor flavor = BTagEntry::FLAV_UDSG;
	if (abs (jets_hadronFlavour[i]) == 5)
	  flavor = BTagEntry::FLAV_B;
	else if (abs (jets_hadronFlavour[i]) == 4)
	  flavor = BTagEntry::FLAV_C;
	float pt_cutoff = std::max (30., std::min (669., double (jets_eta[i])));
	float eta_cutoff = std::min (2.39, fabs (double (jets_eta[i])));
	float weight_cent (1.), weight_UP (1.), weight_DN (1.), weight_FS_UP (1.), weight_FS_DN (1.);
	if (flavor == BTagEntry::FLAV_UDSG) {
	  weight_cent = reader_light->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP = reader_light_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_light_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	else {
	  weight_cent = reader_heavy->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP = reader_heavy_UP->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_DN = reader_heavy_DN->eval (flavor, eta_cutoff, pt_cutoff);
	}
	if (isFastSim) {
	  weight_FS_UP = weight_cent * reader_fastsim_UP->eval (flavor, eta_cutoff, pt_cutoff);	//this is pure fastsimSF
	  weight_FS_DN = weight_cent * reader_fastsim_DN->eval (flavor, eta_cutoff, pt_cutoff);	//this is pure fastsimSF
	  weight_cent *= reader_fastsim->eval (flavor, eta_cutoff, pt_cutoff);
	  weight_UP *= reader_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	  weight_DN *= reader_fastsim->eval (flavor, eta_cutoff, pt_cutoff);	//this is still just btagSF
	}

	btagprob_data *= (1. - weight_cent * eff);
	btagprob_mc *= (1. - eff);
	if (flavor == BTagEntry::FLAV_UDSG) {
	  btagprob_light_UP *= (1. - weight_UP * eff);
	  btagprob_light_DN *= (1. - weight_DN * eff);
	  btagprob_heavy_UP *= (1. - weight_cent * eff);
	  btagprob_heavy_DN *= (1. - weight_cent * eff);
	}
	else {
	  btagprob_light_UP *= (1. - weight_cent * eff);
	  btagprob_light_DN *= (1. - weight_cent * eff);
	  btagprob_heavy_UP *= (1. - weight_UP * eff);
	  btagprob_heavy_DN *= (1. - weight_DN * eff);
	}
	if (isFastSim) {
	  btagprob_FS_UP *= (1. - weight_FS_UP * eff);
	  btagprob_FS_DN *= (1. - weight_FS_DN * eff);
	}
      }
    }
  }
  //Compute a weight as the ratio of data over mc
  Wbtag = btagprob_data / btagprob_mc;
}
