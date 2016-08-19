#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

#include <iostream>

using namespace std;

int main(){

        vector<string> processes;
	theDoctor::SonicScrewdriver sonic;
	sonic.LoadXMLConfig("config.xml");
        //sonic.GetProcessClassTagList(&label);
        //for(uint32_t i = 0; i< label.size(); i++)
        //{cout << label.at(i) << endl;}
	//sonic.SchedulePlots("1DSuperimposed");
        //sonic.WritePlots("./plotsTest2/");
        sonic.GetProcessClassTagList(&processes);
        vector<TH1D*> h;
        h =  sonic.Get1DHistoCloneFromFile("plotsTest", "1DSuperimposed","MET",processes,"CR1l","lepChannel");
	//sonic.ImportHistosFromFile("plotsTest/1DSuperimposed.root");		
	//sonic.MakePlots();
	///TH1D* h = sonic.Get1DHistoClone();
	cout<<"size of histos: " <<h.size()<<endl;
	//cout<< h->GetBinContent(1)<<endl;
	std::cout<<h.at(0)->GetEntries()<<std::endl;
}
