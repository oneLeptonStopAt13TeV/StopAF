#include <exception>
#include <iostream>

#include "TH1F.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"


//usage
//./ZnunuUnctTablePerProcess yieldMorttZ.tab signalRegMor.txt statNames.txt 17 yieldMorZNuNuJECDown.tab yieldMorZNuNuJECUp.tab ttZ notgroup (=nparameters per one signal region, nUnc+1 = 15 currently)
//./ZnunuUnctTablePerProcess yieldMorWZ.tab signalRegMor.txt statNames.txt 17 yieldMorZNuNuJECDown.tab yieldMorZNuNuJECUp.tab WZ notgroup (=nparameters per one signal region, nUnc+1 = 15 currently)
//./ZnunuUnctTablePerProcess yieldMorZZ.tab signalRegMor.txt statNames.txt 17 yieldMorZNuNuJECDown.tab yieldMorZNuNuJECUp.tab ZZ notgroup (=nparameters per one signal region, nUnc+1 = 15 currently)

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 9)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTab = argv[1];
        string inputFile = argv[2];
        string uncNames = argv[3];
        TString nUncs = argv[4];
        string JECDownTab = argv[5];
        string JECUpTab = argv[6];
        string currentProcess = argv[7];
        string group = argv[8];
        int nUnc = nUncs.Atoi();
        if(nUnc%2 !=1)
            cout << "wrong number of systematics!! tables will not be filled correctly" << endl;

        //@MJ@ TODO apply SF
        Figure SF = Figure(1,0);
        if( currentProcess == "ttZ")
           SF = Figure(1.37,0.16);
        else if(currentProcess == "WZ")
           SF = Figure(0.93,0.14);

        vector<string> regions;
        //specify dataset to be read
        vector<string> datasets = {currentProcess}; //@MJ@ TODO change the datasets to meaningful ones

        //read signal regions 
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
       //theDoctor::SonicScrewdriver sonic;
       Table tab(inputTab);

       Table tJECdown(JECDownTab);
       Table tJECup(JECUpTab);
       

        //read names of sysematics 
        vector<string> colI;
	string line2;
        ifstream systfile(uncNames);
        ofstream myfile2;
        myfile2.open ("uncertainties.txt");
        if (systfile.is_open())
        {
            colI.push_back("yield");
            while ( getline (systfile,line2) )
            { 
                if(colI.size() == nUnc)
                    break;
                colI.push_back(line2);
                myfile2 << line2 << endl;
            }
            colI.push_back("jesDN");
            myfile2 << "jesDN" << endl;
            colI.push_back("jesUP");
            myfile2 << "jesUP" << endl;
            colI.push_back("normalizationDN");
            myfile2 << "normalizationDN" << endl;
            colI.push_back("normalizationUP");
            myfile2 << "normalizationUP" << endl;
            systfile.close();
        }
        //get nemaes of signal regions without regions for systematics
        uint32_t rowIdI = 0;
        vector<Double_t> error;
        Double_t result;
        vector<string> realReg;
        ofstream myfile;
        myfile.open ("realregions.txt");
        for(uint32_t e=0; e< regions.size(); e++)
        {
            if(e==0 || (e%nUnc)==0) //nUnct for 8 syst should be 9!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            {
                realReg.push_back(regions.at(e));
                myfile<< regions.at(e) << endl;
            }
        }
        myfile.close();
       //define output histogram
       vector<TH1F*> histo;
       for(uint32_t h =0; h<datasets.size(); h++)
       {
           histo.push_back(new TH1F(datasets.at(h).c_str(), datasets.at(h).c_str(), realReg.size(), 0, realReg.size()));
       }

       //output syst names
       vector<string> systOutNames;
       for(uint32_t s=0; s<colI.size(); s++)
       {
               //systOutNames.push_back("yield");
           if(s==0 || s%2 == 1)
               systOutNames.push_back(colI.at(s));
       }
       //tables for higher uncertainty and summary table
       Table thigher(systOutNames, realReg);
       vector<string> uncRange = {"syst low (%)", "syst high (%)"};
       Table tsummary( uncRange, systOutNames);

       
       //tables for all realtive and absolute uncertainties
       Table tI(colI, realReg);
       Table trel(colI, realReg);
       uint32_t uncLine = 0;
       float uncTot = 0;
       float uncBef = 0;      
       float uncHigher = 0;      
       float binValue = 0;
       uint32_t l = 0;
       theDoctor::Figure yieldVal;

       //loop over all regions
       for(uint32_t r=0; r<regions.size();r++)
       {
           theDoctor::Figure resall = tab.Get(regions.at(r), datasets.at(0)) * SF;
           
           //first there is real region value
           if(r == 0 || r%nUnc == 0)
           {
               l=0;
               uncLine++;
               yieldVal = resall;
               histo.at(0)->SetBinContent(uncLine,resall.value());
               if(r != 0)
                   histo.at(0)->SetBinError(uncLine-1,uncTot);
               uncTot = 0;
               uncBef = 0;

               if(resall.value() !=0) 
                   thigher.Set(systOutNames.at(0), realReg.at(uncLine-1), Figure((resall.error()/resall.value())*100, 0));
               else
                   thigher.Set(systOutNames.at(0), realReg.at(uncLine-1), Figure(-1, 0));
      
               l++;
           }
           //then all systematics
           else
           {
               if(uncBef == 0)
               {
                   uncBef = resall.value();
               }
               else
               { 
                   uncHigher = abs(uncBef-yieldVal.value()) > abs(resall.value() - yieldVal.value()) ? abs(uncBef-yieldVal.value()): abs(resall.value() - yieldVal.value());   
                   uncTot += uncHigher*uncHigher ;//@MJ@TODO now average unc, change to maximal?!

                
                   if(systOutNames.at(l) == "lepSFDN" || systOutNames.at(l) == "lepSFUP")
                   {
                       if(yieldVal.value() != 0)
                           thigher.Set(systOutNames.at(l), realReg.at(uncLine-1), Figure((uncHigher/yieldVal.value())*100, 0));
                       else
                           thigher.Set(systOutNames.at(l), realReg.at(uncLine-1), Figure(0, 0));

                   }
        	   else
                   {
                   
                       float center;
                       if(uncBef+resall.value() != 0)
                           center = (abs(uncBef-resall.value()))/(uncBef+resall.value());
                       else
                           center = 0;

                       thigher.Set(systOutNames.at(l), realReg.at(uncLine-1), Figure((center)*100, 0));

                   }
                   l++;
                   uncBef = 0;


               }
           }
           //fill the tables
           //cout << "colI " << colI.at(r-((uncLine-1)*(nUnc))) << endl;
	   tI.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), resall);
           Figure g = (resall-yieldVal)/yieldVal;
           if(r == 0 || r%nUnc == 0)
           {
               if(yieldVal.value() != 0)
                   trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure((yieldVal.error()/yieldVal.value())*100, 0));
               else
                   trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(-1,0));
               //thigher.Set(systOutNames.at(0), realReg.at(uncLine-1), Figure((abs(yieldVal.error()/yieldVal.value())*100, 0)));
           }
           else
           {
               if(yieldVal.value() != 0)
                   trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(abs(g.value())*100,abs(0)));
               else
                   trel.Set(colI.at(r-((uncLine-1)*(nUnc))), realReg.at(uncLine-1), Figure(-1,0));
           }
           //cout << realReg.at(uncLine-1) << " value " <<  colI.at(r-((uncLine-1)*(nUnc)))<< endl;
       }
       //fill JEC + norm
       for(uint32_t j = 0; j<realReg.size(); j++)
       {     
           theDoctor::Figure JECd = tJECdown.Get(realReg.at(j), currentProcess ) * SF;
           theDoctor::Figure JECu = tJECup.Get(realReg.at(j), currentProcess) * SF;
   
           Figure SFDN = Figure(SF.value()-SF.error(),0 );
           Figure SFUP = Figure(SF.value()+SF.error(),0 );
           Figure nDN = tab.Get(realReg.at(j), currentProcess ) * SFDN;
           Figure nUP = tab.Get(realReg.at(j), currentProcess ) * SFUP;
           Figure relZZNorm(0.06,0);
           Figure nyield = tI.Get("yield", realReg.at(j));

           if(currentProcess == "ZZ")
           {
               nDN =  nyield - (nyield*relZZNorm);
               nUP = nyield + (nyield*relZZNorm);
           }
           tI.Set("normalizationDN", realReg.at(j), nDN);
           tI.Set("normalizationUP", realReg.at(j), nUP);

           
           //JES yields
           tI.Set("jesDN", realReg.at(j), JECd);
           tI.Set("jesUP", realReg.at(j), JECu);

           //relative unc
           Figure a;
           Figure c;
           Figure an;
           Figure cn;
           if(nyield.value() > 0)
           {
               a = (JECd - nyield)/(nyield);
               c = (JECu - nyield)/(nyield);
               an = (nDN - nyield)/(nyield);
               cn = (nUP - nyield)/(nyield);
           }
           else
           {
               a = -0.01;
               c = -0.01;
               an = -0.01;
               cn = -0.01;
           }

           Figure h(-0.01,0);
           if((JECd.value() + JECu.value()) !=0 )
               h = Figure((abs(JECu.value()-JECd.value()))/(JECu.value()+JECd.value()), 0);

          
           Figure hn(-0.01,0);
           hn = an.value() > cn.value() ? an :cn;
             
           //@MJ@ TODO finish it here!!!
           trel.Set("jesDN", realReg.at(j), Figure(abs(a.value())*100,abs(0)));
           trel.Set("jesUP", realReg.at(j), Figure(abs(c.value())*100,abs(0)));
           trel.Set("normalizationDN", realReg.at(j), Figure(abs(an.value())*100,abs(0)));
           trel.Set("normalizationUP", realReg.at(j), Figure(abs(cn.value())*100,abs(0)));

           thigher.Set("jesDN", realReg.at(j), Figure(h.value()*100,0));
           thigher.Set("normalizationDN", realReg.at(j), Figure(hn.value()*100,0));
           
       }

        //recompute PU and set the normalization unc
        Figure PUDNtot(0.0);
        Figure PUUPtot(0.0);
        //Figure PUyieldtot(0.0);
        for(uint32_t p=0; p<realReg.size(); p++)
        {
            PUDNtot += tI.Get("PUdown", realReg.at(p));
            PUUPtot += tI.Get("PUup", realReg.at(p));
            //PUyieldtot += tI.Get("yield", realReg.at(p));
        }
        Figure two(2,0);
        Figure perDiff = (PUUPtot-PUDNtot)/(PUUPtot+PUDNtot);
        Figure fPerDiff = Figure(abs(perDiff.value()), perDiff.error());
        Figure PUY(0.0);
        for(uint32_t p=0; p<realReg.size(); p++)
        {
            
            PUY = tI.Get("yield", realReg.at(p));
            tI.Set("PUdown", realReg.at(p), PUY - (PUY*fPerDiff));
            tI.Set("PUup", realReg.at(p), PUY + (PUY*fPerDiff));
            trel.Set("PUdown", realReg.at(p), Figure(fPerDiff.value()*100,0));
            trel.Set("PUup", realReg.at(p),  Figure(fPerDiff.value()*100,0));
            thigher.Set("PUdown", realReg.at(p),  Figure(fPerDiff.value()*100,0));

        }

        //group values to compute btetter yields
        //first compute the grouped yields and put them to vector

        if(group == "group")
        {
        //vector of groups
        vector<string> groups = {"_A_", "_B_", "_C_", "_D_", "_E_", "_F_", "_G_", "_H_" };
        vector<string> systDN = { "pdfDN", "alphaSDN", "Q2DN" ,"jesDN"};
        vector<string> systUP = { "pdfUP", "alphaSUP", "Q2UP" ,"jesUP"};
        uint32_t begin = 0;
        uint32_t end = 0;
        uint32_t nGroup = 0;
        Figure gYield(0,0);
        Figure gDN(0,0);
        Figure gUP(0,0);

        //vector<float> groupedYields;
        for(uint32_t u=0; u<systDN.size(); u++)
        {
            string substr = groups.at(nGroup);
            for(uint32_t r2 = 0; r2<realReg.size(); r2++)
            {
                if (realReg.at(u).find(substr) == std::string::npos) 
                {
                    nGroup++;
                    substr = groups.at(nGroup);
                    end = r2;
                    for(uint32_t g = begin; g< end; g++)
                    {
                        //set
                    }
                    begin = r2;
                    gYield = Figure(0,0);
                    gDN = Figure(0,0);
                    gUP = Figure(0,0);
                
                }
                gYield += tI.Get("yield", realReg.at(r2));
                gDN += tI.Get(systDN.at(u), realReg.at(r2));
                gUP += tI.Get(systUP.at(u), realReg.at(r2));

            }
        }
        }
        //now find the extreme values
        for(uint32_t s=0; s<systOutNames.size(); s++)
        {
            vector<float> oneSysts;
            for(uint32_t t=0; t<realReg.size(); t++)
            {
                Figure oneSyst = thigher.Get(systOutNames.at(s), realReg.at(t));
                oneSysts.push_back(oneSyst.value());
            }
            sort(oneSysts.begin(), oneSysts.end());
            for(uint32_t z=0; z<oneSysts.size(); z++)
            {
                if(oneSysts.at(z) != 0 )
                {
                    tsummary.Set(uncRange.at(0), systOutNames.at(s), oneSysts.at(z));
                    break;
                }
            }
            tsummary.Set(uncRange.at(1), systOutNames.at(s), oneSysts.at(realReg.size()-1));
        }

        //write output tables
        if(group == "group")
            currentProcess+=group;
        tI.Print(static_cast<string>(currentProcess+"tableUnc.tab"),4);
        tI.PrintLatex(static_cast<string>(currentProcess+"tableUnc.tex"),4);
        trel.Print(static_cast<string>(currentProcess+"tableUncRel.tab"),4, "noError");
        trel.PrintLatex(static_cast<string>(currentProcess+"tableUncRel.tex"),4, "noError");
        thigher.Print(static_cast<string>(currentProcess+"tableUncRelHiger.tab"),4, "noError");
        thigher.PrintLatex(static_cast<string>(currentProcess+"tableUncRelHigher.tex"),4, "noError");
        tsummary.Print(static_cast<string>(currentProcess+"tableUncRelSummaryTab.tab"),4, "noError");
        tsummary.PrintLatex(static_cast<string>(currentProcess+"tableUncRelSummaryTab.tex"),4, "noError");


       //TFile fi2("Znunuuncertainties.root","RECREATE");
       //histo.at(0)->Write();
       //fi2.Close();
}
