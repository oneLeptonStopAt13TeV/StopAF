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
#include "TApplication.h"
#include "TGraph.h"

#include <fstream>
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
    TString dir = in_stream.at(1+2*nrOfVars);    


    //this will be read from steering
    vector<TString> inputRootFile;
    //this will be read from steering
    vector<TString> variable;

    for(uint32_t f = 0; f < nrOfVars; f++)
    {
        variable.push_back(in_stream.at(1+f));
        inputRootFile.push_back(in_stream.at(1+nrOfVars+f));
    }

    vector<TFile*> inputTFiles;
    vector<TH1D*> inputHist;

    for(uint32_t files = 0; files < inputRootFile.size(); files++)
    {
        inputTFiles.push_back(new TFile(inputRootFile.at(files)));
        cout << "here1" << endl;
        TH1D* some = reinterpret_cast<TH1D*>(inputTFiles.at(files)->Get("muon/preselection/" + variable.at(files))->Clone()); //@MJ@ TODO this is not good
        inputHist.push_back(some);
        //TH1D* hist;
        //inputTFiles.at(files)->GetObject("muon/preselection/" + variable.at(files), hist); 
        //inputHist.push_back(hist);
    
    }

    cout << "nr of input files; " << inputRootFile.size() << endl;
    
    /*TCanvas* can = new TCanvas("can", "bla", 600, 600);
    can->Divide(1,2);
    can->cd(1);
    */

    if(process == "naddition" && inputRootFile.size() > 1)
    {        
        for(uint32_t i = 0; i < inputHist.size(); i++)
        {  
            Double_t norm = inputHist.at(i)->Integral(); 
            cout << "bin content" << inputHist.at(i)->GetBinContent(5) << inputHist.at(i)->GetBinContent(6) << inputHist.at(i)->GetBinContent(7) << endl; 
            cout << "norm " << norm << endl;
            inputHist.at(i)->Scale(9000, "width");//(1/norm);
            cout << inputHist.at(i)->Integral() << endl;
        }
        
        for(uint32_t h = 0; h < inputHist.size(); h++)
        {  //protect if only one hist
            cout << "adding histos: " << inputHist.size() << endl;
           //inputHist.at(0)->Add(inputHist.at(1+h), 1);
        }

        //inputHist.at(0)->GetYaxis()->SetTitle("whatever");
    }

    if(process == "merge")
    {   
        for(uint32_t h = 0; h <inputRootFile.size(); h++ )
        {   
            if(h == 0)
                inputHist.at(h)->Draw();
            else
                inputHist.at(h)->Draw("same");
        }
    }

    if(process == "bsefficiency")
    {
        if(inputRootFile.size() == 2)
        {
            Int_t nrOfBins = inputHist.at(0)->GetNbinsX();
            TGraph* output = new TGraph(nrOfBins);
            cout << "nr of bins: " << nrOfBins << endl;  
            
            for(uint32_t bin = 1; bin<(nrOfBins+1); bin++)
            {
                output->SetPoint(bin, inputHist.at(0)->Integral(bin, nrOfBins), inputHist.at(1)->Integral(bin, nrOfBins));
            }
            
           output->Draw("AC*"); 
        }
    }
    
    /*TCanvas* can = new TCanvas("can", "bla", 600, 600);
    can->Divide(1,2);
    can->cd(1);
    */
    /*
    cout << "here2" << endl;
    cout << "integral: " << inputHist.at(0)->Integral() << endl;
    inputHist.at(0)->Draw();
    //inputHist.at(1)->Draw("same hist e");
    //can->SaveAs(process + variable.at(0) + ".eps");
    cout << "drawing" << endl;
    */
    
    TString rootOut = process + variable.at(0) + ".root";
    TFile fi(rootOut,"RECREATE");
    cout << "saving" << endl;
	
    inputHist.at(0)->Write();
    fi.Close();
}
