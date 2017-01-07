#include <exception>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "../sonicScrewdriver/interface/SonicScrewdriver.h"
#include "../sonicScrewdriver/interface/Table.h"

//usage
//vim uncertainties.txt
//./ZnunuUnctSummedTables ttZtableUnc.tab WZtableUnc.tab ZZtableUnc.tab realregions.txt uncertainties.txt

using namespace std;
using namespace theDoctor;

int main(int argc, char *argv[]){

        if(argc != 6)
            throw std::runtime_error("Bad number of arguments!");
         
        //get info from arguments
        string inputTab1 = argv[1];
        string inputTab2 = argv[2];
        string inputTab3 = argv[3];
        string inputFile = argv[4];
        string uncNames = argv[5];

        vector<string> regions;
        float lumi = 0.062;
        //specify dataset to be read

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
     
        vector<string> regLabels = regions;
        const std::string ss = "less";
        const std::string t = "<";
        
        for(uint32_t reg =0; reg<regLabels.size(); reg++)
        {
            std::string::size_type n = 0;
            while ( ( n = regLabels.at(reg).find( ss, n ) ) != std::string::npos )
            {
                regLabels.at(reg).replace( n, ss.size(), t );
                n += t.size();
            } 
        }              

       Table tab1(inputTab1);
       Table tab2(inputTab2);
       Table tab3(inputTab3);



        //read names of sysematics 
        vector<string> colI;
	string line2;
        ifstream systfile(uncNames);
        if (systfile.is_open())
        {
            colI.push_back("yield");
            //colI.push_back("mcStatsDN");
            //colI.push_back("mcStatsUP");
            //colI.push_back("lumiDN");
            //colI.push_back("lumiUP");
            //if(currentProcess == "ZZ")
            //{
            //    colI.push_back("crossSectionDN");
            //    colI.push_back("crossSectionUP");
            //}
       
            while ( getline (systfile,line2) )
            {
              /*  if(line2 == "topPtmodeling" && currentProcess != "ttZ") 
                    continue;
                if(line2 == "topPtmodeling2") 
                    continue;*/
                colI.push_back(line2);
            }
            systfile.close();
        }
      
       vector<string> systOutNames;
       for(uint32_t s=0; s<colI.size(); s++)
       {
               if(s==0 || s%2 == 1)
                   systOutNames.push_back(colI.at(s));
       }
 
        vector<string> systOutNamesAll = systOutNames;
        systOutNamesAll.push_back("luminosity");
        //systOutNamesAll.push_back("pile-up");
        //systOutNamesAll.push_back("ZZ normalization");
 
        Table tabSum( colI, regions );
        Table thigher(systOutNames, regions);
        vector<string> uncRange = {"syst low (%)", "syst high (%)"};
        Table tsummary( uncRange, systOutNamesAll);
        vector<string> processes = {"ttZ", "WZ", "ZZ", "ZNuNu"};
        Table totYields(processes, regions );
        Table totYieldsStatSyst(processes, regions );

        for(uint32_t r = 0; r<regions.size(); r++)
        {
            Figure y1;
            Figure y2;
            Figure y3;
            float unc1;
            float unc2;
            float unc3;
            float bef1;
            float bef2;
            float bef3;
            float uncHigher1;
            float uncHigher2;
            float uncHigher3;
            float unc1Stat;
            float unc2Stat;
            float unc3Stat;
            float unc1Syst;
            float unc2Syst;
            float unc3Syst;
            for(uint32_t c = 0; c< colI.size(); c++)
            {
                Figure f1 = tab1.Get(colI.at(c), regions.at(r));
                Figure f2 = tab2.Get(colI.at(c), regions.at(r));
                Figure f3 = tab3.Get(colI.at(c), regions.at(r));

                Figure f = f1+f2+f3;
                tabSum.Set(colI.at(c), regions.at(r),f );

                if(c==0)
                {
                    y1 = f1;
                    y2 = f2;
                    y3 = f3;
                    unc1 = 0;
                    unc2 = 0;
                    unc3 = 0;
                    unc1Stat = 0;
                    unc2Stat = 0;
                    unc3Stat = 0;
                    unc1Syst = 0;
                    unc2Syst = 0;
                    unc3Syst = 0;

                    unc1 += y1.error()*y1.error();
                    unc2 += y2.error()*y2.error();
                    unc3 += y3.error()*y3.error();
                    unc1Stat = y1.error()*y1.error();
                    unc2Stat = y2.error()*y2.error();
                    unc3Stat = y3.error()*y3.error();

                }
                else if(c%2 != 0)
                {
                    bef1 = f1.value();
                    bef2 = f2.value();
                    bef3 = f3.value();
                }
                else if(c%2 ==0)
                {
                   uncHigher1 = abs(bef1-y1.value()) > abs(f1.value() - y1.value()) ? abs(bef1-y1.value()): abs(f1.value() - y1.value());
                   uncHigher2 = abs(bef2-y2.value()) > abs(f2.value() - y2.value()) ? abs(bef2-y2.value()): abs(f2.value() - y2.value());
                   uncHigher3 = abs(bef3-y3.value()) > abs(f3.value() - y3.value()) ? abs(bef3-y3.value()): abs(f3.value() - y3.value());

                   if(colI.at(c) == "lepSFDN" || colI.at(c) == "lepSFUP")
                   {
                       unc1 += uncHigher1*uncHigher1;
                       unc2 += uncHigher2*uncHigher2;
                       unc3 += uncHigher3*uncHigher3;
                       unc1Syst += uncHigher1*uncHigher1;
                       unc2Syst += uncHigher2*uncHigher2;
                       unc3Syst += uncHigher3*uncHigher3;
                   }
                   else
                   {
                       unc1 += (abs(bef1-f1.value())/2) * (abs(bef1-f1.value())/2) ;
                       unc2 += (abs(bef2-f2.value())/2) * (abs(bef2-f2.value())/2) ;
                       unc3 += (abs(bef3-f3.value())/2) * (abs(bef3-f3.value())/2) ;
                       unc1Syst += (abs(bef1-f1.value())/2) * (abs(bef1-f1.value())/2) ;
                       unc2Syst += (abs(bef2-f2.value())/2) * (abs(bef2-f2.value())/2) ;
                       unc3Syst += (abs(bef3-f3.value())/2) * (abs(bef3-f3.value())/2) ;

                   }

                }
                else
                {
                    cout << "something went terribly wrong" << endl;
                }

                if(c == colI.size() -1)
                {
                    unc1 += (lumi*y1.value())*(lumi*y1.value());
                    unc2 += (lumi*y2.value())*(lumi*y2.value());
                    unc3 += (lumi*y3.value())*(lumi*y3.value());
                    unc1Syst += (lumi*y1.value())*(lumi*y1.value());
                    unc2Syst += (lumi*y2.value())*(lumi*y2.value());
                    unc3Syst += (lumi*y3.value())*(lumi*y3.value());

                    totYields.Set(processes.at(0), regions.at(r), Figure(y1.value(),sqrt(unc1)));
                    totYields.Set(processes.at(1), regions.at(r), Figure(y2.value(),sqrt(unc2)));
                    totYields.Set(processes.at(2), regions.at(r), Figure(y3.value(),sqrt(unc3)));
                    cout << "stat " << sqrt(unc1Stat) << " syst " << sqrt(unc1Syst) << endl;
                    totYieldsStatSyst.Set(processes.at(0), regions.at(r), Figure(y1.value(),sqrt(unc1Stat),sqrt(unc1Syst)));
                    totYieldsStatSyst.Set(processes.at(1), regions.at(r), Figure(y2.value(),sqrt(unc2Stat),sqrt(unc2Syst)));
                    totYieldsStatSyst.Set(processes.at(2), regions.at(r), Figure(y3.value(),sqrt(unc3Stat),sqrt(unc3Syst)));
                }
            }
        }

        float uncBef = 0;
        float uncHigher = 0;
        for(uint32_t r=0; r<regions.size();r++)
        {
            Figure yieldVal;
            float uncTot = 0;
            float uncTotStat = 0;
            float uncTotSyst = 0;
            for(uint32_t c = 0; c< colI.size(); c++)
            {
                Figure resall = tabSum.Get(colI.at(c), regions.at(r));
                
                if(c==0)
                {
                    yieldVal = resall;
                    uncTot += yieldVal.error()*yieldVal.error() ;
                    uncTotStat = yieldVal.error()*yieldVal.error() ;
                    if(resall.value() !=0)
                        thigher.Set(systOutNames.at(0), regions.at(r), Figure((resall.error()/resall.value())*100, 0));
                    else
                        thigher.Set(systOutNames.at(0), regions.at(r), Figure(-1, 0));

                }
                else if(c%2 != 0)
                {
                    uncBef = resall.value();
                }
                else if(c%2 == 0)
                {
                   uncHigher = abs(uncBef-yieldVal.value()) > abs(resall.value() - yieldVal.value()) ? abs(uncBef-yieldVal.value()): abs(resall.value() - yieldVal.value());


                   if(colI.at(c) == "lepSFDN" || colI.at(c) == "lepSFUP")
                   {
                       uncTot += uncHigher*uncHigher ;
                       uncTotSyst += uncHigher*uncHigher ;
                       if(yieldVal.value() != 0)
                           thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((uncHigher/yieldVal.value())*100, 0));
                       else
                           thigher.Set(systOutNames.at(c/2), regions.at(r), Figure(0, 0));

                   }
                   else
                   {

                       float center;
                       if(uncBef+resall.value() != 0)
                           center = (abs(uncBef-resall.value()))/(uncBef+resall.value());
                       else
                           center = 0;
                       
                       uncTot += (abs(uncBef-resall.value())/2) * (abs(uncBef-resall.value())/2) ;
                       uncTotSyst += (abs(uncBef-resall.value())/2) * (abs(uncBef-resall.value())/2) ;

                       thigher.Set(systOutNames.at(c/2), regions.at(r), Figure((center)*100, 0));

                   }

                   if(c==colI.size() -1)
                       cout << "yield" << yieldVal.value() << "uncertainty is" << sqrt(uncTot) << endl;

                }
                else
                {
                    cout << "there is a mistake somewhere " << endl;
                }
                if(c == colI.size() -1)
                {
                    uncTot += (lumi*yieldVal.value())*(lumi*yieldVal.value());
                    uncTotSyst += (lumi*yieldVal.value())*(lumi*yieldVal.value());
                    totYields.Set(processes.at(3), regions.at(r), Figure(yieldVal.value(),sqrt(uncTot)));
                    totYieldsStatSyst.Set(processes.at(3), regions.at(r), Figure(yieldVal.value(),sqrt(uncTotStat),sqrt(uncTotSyst)));
                }

            }

        }

       //summary plots
        for(uint32_t s=0; s<systOutNames.size(); s++)
        {
            vector<float> oneSysts;
            for(uint32_t x=0; x<regions.size(); x++)
            {
                Figure oneSyst = thigher.Get(systOutNames.at(s), regions.at(x));
                oneSysts.push_back(oneSyst.value());
            }
            sort(oneSysts.begin(), oneSysts.end());
            for(uint32_t z=0; z<oneSysts.size(); z++)
            {
                if(oneSysts.at(z) != 0 )
                {
                    tsummary.Set(uncRange.at(0), systOutNamesAll.at(s), oneSysts.at(z));
                    break;
                }
            }
            tsummary.Set(uncRange.at(1), systOutNamesAll.at(s), oneSysts.at(regions.size()-1));
        }
        tsummary.Set(uncRange.at(1), systOutNamesAll.at(systOutNames.size()), Figure(6.2,0));

        //@MJ@ TODO add lumi, CS

 
        tabSum.Print(static_cast<string>("summedprocessesZnunu.tab"),4);
        tabSum.PrintLatex(static_cast<string>("summedprocessesZnunu.tex"),4); 
        thigher.Print(static_cast<string>("summedprocessesRelHigherZnunu.tab"),4, "noError");
        thigher.PrintLatex(static_cast<string>("summedprocessesRelHigherZnunu.tex"),4, "noError"); 
        tsummary.Print(static_cast<string>("summedprocessesRelHigherZnunuSummary.tab"),4, "noError");
        tsummary.PrintLatex(static_cast<string>("summedprocessesRelHigherZnunuSummary.tex"),4, "noError"); 
        totYields.Print(static_cast<string>("summedprocessesZnunuYieldsAndErrors.tab"),4);
        totYields.PrintLatex(static_cast<string>("summedprocessesZnunuYieldsAndErrors.tex"),4); 
        totYieldsStatSyst.Print(static_cast<string>("summedprocessesZnunuYieldsAndErrorsStatSystBreakdown.tab"),4, "systError");
        totYieldsStatSyst.PrintLatex(static_cast<string>("summedprocessesZnunuYieldsAndErrorsStatSystBreakdown.tex"),4, "systError"); 


}
