#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THashList.h"

#include <iostream>
#include <vector>

void AddDivideMultiply()
{
    //steering file with information
    TString input = gApplication->Argv(5);
    //information if addition, division or multiplication required
    TString process = gApplication->Argv(6);
    
    vector<TString> in_stream;

    ifstream is;
    is.open(input);
    if( is.is_open() )
    {
	char * sline;
	TString tsline("");
	Int_t line = 0;
	string strline;
	while(getline(is, strline)) {
		tsline = strline;
		if(tsline.BeginsWith(" ")) continue;
		if(tsline.Length()==0)     continue;
		in_stream.push_back(tsline);
		line++;
	}
	is.close();
    }

    uint32_t nrOfVars = in_stream.at(0).Atof();

    TCanvas* can = new TCanvas("can", "can");
    can->Divide(1,1);

    //this will be read from steering
    vector<TString> inputRootFile;
    //this will be read from steering
    vector<TString> variable;

    for(uint32_t f = 0; f < nrOfVars; f++)
    {
        variable.push_back(in_stream.at(1+f));
        inputRootFile.push_back(in_stream.at(1+nrOfVars+f));
    }

    vector<TFile> inputTFiles;
    vector<TH1D*> inputHist;

    for(uint32_t files; files < inputRootFile.size(); files++)
    {
        inputTFiles.push_back(new TFile(inputRootFile.at(file)));
        inputHist.push_back((TH1D*)data->Get(variable.at(file))->Clone());      
        
    }


}
