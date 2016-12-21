#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//./ZnunuUnctTablePU yieldMorPU.tab signalRegMorPU.txt statNamesPU.txt 3 yieldJECDownMor.tab yieldJECUpMor.tab(=nparameters per one signal region, nUnc+1 = 16 currently)

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 7)
            throw std::runtime_error("Bad number of arguments!");
         
         

        string inputTab = argv[1];
        string inputFile = argv[2];
        string uncNames = argv[3];
        TString nUncs = argv[4];
        string JECDownTab = argv[5];
        string JECUpTab = argv[6];
        int nUnc = nUncs.Atoi();

        vector<string> regions;
        vector<string> datasets = {"Znunu"}; //@MJ@ TODO change the datasets to meaningful ones

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
      
       if(((regions.size())%nUnc) != 0)
           throw std::runtime_error("wrong number of regions");

       TH1::SetDefaultSumw2();
       Table tab(inputTab);

       Table tJECdown(JECDownTab);
       Table tJECup(JECUpTab);
       

  
        vector<string> colI;
        
	string line2;
        ifstream systfile(uncNames);
        if (systfile.is_open())
        {
            colI.push_back("value");
            while ( getline (systfile,line2) )
            { 
                if(colI.size() == nUnc)
                    break;
                colI.push_back(line2);
            }
            systfile.close();
        }
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
       Table trel(colI, realReg);
       uint32_t uncLine = 0;
       float uncTot = 0;
       float uncBef = 0;      
       float uncHigher = 0;      
       float binValue = 0;
       theDoctor::Figure yieldVal;
       for(uint32_t r=0; r<regions.size();r++)
       {
           theDoctor::Figure resall = tab.Get(regions.at(r), datasets.at(0));

           if(r == 0 || r%nUnc == 0)
           {
               uncLine++;
               yieldVal = resall;
               histo.at(0)->SetBinContent(uncLine,resall.value());
               if(r != 0)
                   histo.at(0)->SetBinError(uncLine-1,uncTot);
               uncTot = 0;
               uncBef = 0;
 
           }
           else
           {
               if(uncBef == 0)
               {
                   uncBef = resall.value();
               }
               else
               { 
                   binValue = histo.at(0)->GetBinContent(uncLine);
                   uncHigher = abs(uncBef-binValue) > abs(resall.value() - binValue) ? abs(uncBef-binValue): abs(resall.value() - binValue);   
                   uncTot += uncHigher*uncHigher ;//@MJ@TODO now average unc, change to maximal?!
                   uncBef = 0;
               }
           }
           cout << "colI " << colI.at(r-((uncLine-1)*(nUnc))) << endl;
	   tI.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), resall);
           Figure g = (resall-yieldVal)/yieldVal;
           if(r == 0 || r%nUnc == 0)
           {
           trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(abs(yieldVal.error()/yieldVal.value())*100,abs(0)));
           }
           else
           {
           trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(abs(g.value())*100,abs(0)));
           }
           cout << realReg.at(uncLine-1) << " value " <<  colI.at(r-((uncLine-1)*(nUnc)))<< endl;
       }
        tI.Print(static_cast<string>("tableUncZnunuPU.tab"),2);
        tI.PrintLatex(static_cast<string>("tableUncZnunuPU.tex"),2);
        trel.Print(static_cast<string>("tableUncZnunuRelPU.tab"),2);
        trel.PrintLatex(static_cast<string>("tableUncZunuRelPU.tex"),2);

       TFile fi2("ttZuncertaintiesPU.root","RECREATE");
       histo.at(0)->Write();
       fi2.Close();
}
