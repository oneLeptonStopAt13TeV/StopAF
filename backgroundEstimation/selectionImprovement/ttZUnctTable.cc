#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//./ttZUnctTable yield.tab signalReg.txt nUnc+1(=nparameters per one signal region)

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 4)
            throw std::runtime_error("Bad number of arguments!");
         
         

        string inputTab = argv[1];
        string inputFile = argv[2];
        TString nUncs = argv[3];
        int nUnc = nUncs.Atoi();

        vector<string> regions;
        vector<string> datasets = {"ttZ"}; //@MJ@ TODO change the datasets to meaningful ones

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
      
       if(((regions.size())%9) != 0)
           throw std::runtime_error("wrong number of regions");

       TH1::SetDefaultSumw2();
       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

  
        vector<string> colI = {"1", "2", "3", "4", "5", "6", "7", "8", "9"}; //@MJ@ TODO real names //@MJ@ TODO do some checks taht I am not owerflowing the table
        uint32_t rowIdI = 0;
        vector<Double_t> error;
        Double_t result;
        vector<string> realReg;
        for(uint32_t e=0; e< regions.size(); e++)
        {
            if(e==0 || (e%nUnc)==0) //nUnct for 8 syst should be 9!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                realReg.push_back(regions.at(e));
        }
       vector<TH1F*> histo;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), realReg.size(), 0, realReg.size()));
       }
 
       Table tI(colI, realReg);
       uint32_t uncLine = 0;
       float uncTot = 0;
       float uncBef = 0;      
       float binValue = 0;

       for(uint32_t r=0; r<regions.size();r++)
       {
           theDoctor::Figure resall = tab.Get(regions.at(r), datasets.at(0));

           if(r == 0 || r%nUnc == 0)
           {
               uncLine++;
               histo.at(0)->SetBinContent(uncLine,resall.value());
               if(r != 0)
                   histo.at(0)->SetBinError(uncLine-1,uncTot); //fill preceeding error
           }
           else
           {
               if(uncBef == 0)
               {
                   uncBef = resall.value();
               }
               else
               {
                   uncTot += (abs(resall.value() - uncBef)* abs(resall.value() - uncBef)) ;//@MJ@TODO now average unc, change to maximal?!
                   uncBef = 0;
               }
           }
           cout << "colI " << r-((uncLine-1)*(nUnc)) << endl;
           tI.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), resall);
           cout << realReg.at(uncLine-1) << " value " <<  colI.at(r-((uncLine-1)*(nUnc)))<< endl;
       }
        tI.Print(static_cast<string>("tableUnc.tab"));
        tI.PrintLatex(static_cast<string>("tableUnc.tex"));

       //TCanvas *can = new TCanvas("can","can");
       //can->cd();
       TFile fi2("ttZuncertainties.root","RECREATE");
       histo.at(0)->Write();
       /*for(uint32_t c=0; c<histo.size(); c++)
       {
           histo.at(c)->SetTitle("LowPU/HighPU");
           histoH.at(c)->SetTitle("HighPU/all");
           histoL.at(c)->SetTitle("LowPU/all");
           for(uint32_t b=0; b<regions.size(); b++)
           {
               histo.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
               histoH.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
               histoL.at(c)->GetXaxis()->SetBinLabel(b+1,regions.at(b).c_str());
           }
           //if(c==0)
           //{
           //   histo.at(c)->Draw("hist e");
           //}
           //else
           //{               
           //histo.at(c)->SetLineColor(c+4);
           //histoL.at(c)->SetLineColor(c+4);
           //histoH.at(c)->SetLineColor(c+4);
              //histo.at(c)->Draw("same hist e");
           //}
           histo.at(c)->Write();
           histoH.at(c)->Write();
           histoL.at(c)->Write();
       }*/
       fi2.Close();
       //can->SaveAs("plotPUyield.eps"); //@MJ@ TODO svae to root file
}
