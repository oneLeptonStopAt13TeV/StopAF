#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//./PUYields signalRegPU.tab datacards/signalReg.txt

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 3)
            throw std::runtime_error("Bad number of arguments!");

        string inputTab = argv[1];
        TString inputFile = argv[2];

        vector<string> regions;
        vector<string> datasets = {"500.000000_400.000000", "1000.000000_1.000000", "1000.000000_650.000000"}; //@MJ@ TODO change the datasets to meaningful ones

	string line;
        ifstream regfile(inputFile);
        if (regfile.is_open())
        {
            while ( getline (regfile,line) )
            {
                regions.push_back(line);
            }
            regfile.close();
        }

       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

       vector<TH1F*> histo;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), regions.size(), 0, regions.size()));
       }
       for(uint32_t r=0; r<regions.size();r++)
       {
           for(uint32_t d=0; d<datasets.size();d++)
           {
               theDoctor::Figure resLPU = tab.Get(regions.at(r)+"LowPU", datasets.at(d));
               theDoctor::Figure resHPU = tab.Get(regions.at(r)+"HighPU", datasets.at(d));
               Figure res = resLPU/resHPU;
               histo.at(d)->SetBinContent(r+1,res.value());
               histo.at(d)->SetBinError(r+1,res.error());
 
           }
       }

       TCanvas *can = new TCanvas("can","can");
       can->cd();
       TFile fi2("PUYields.root","RECREATE");
       for(uint32_t c=0; c<histo.size(); c++)
       {
           histo.at(c)->SetTitle("LowPU/HighPU");
           for(uint32_t b=0; b<regions.size(); b++)
           {
               histo.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
           }
           if(c==0)
           {
              histo.at(c)->Draw("hist e");
                  
           }
           else
           {               
              histo.at(c)->SetLineColor(c+4);
              histo.at(c)->Draw("same hist e");
           }
           histo.at(c)->Write();
       }
       fi2.Close();
       can->SaveAs("plotPUyield.eps"); //@MJ@ TODO svae to root file
}
