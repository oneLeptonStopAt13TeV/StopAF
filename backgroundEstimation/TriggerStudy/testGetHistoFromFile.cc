//to get histogram(s) from file and canvas
#include "../../sonicScrewdriver/interface/SonicScrewdriver.h"

#include <iostream>

using namespace std;

int main(){

        vector<string> processes;
	theDoctor::SonicScrewdriver sonic;
        sonic.LoadXMLConfig("config.xml");
        sonic.GetProcessClassTagList(&processes);
        vector<TH1D*> h;
        h = sonic.Get1DHistoCloneFromFile("plotsTest", "1DSuperimposed","MET",processes,"CR1l","lepChannel");
        cout<<"size of histos: " <<h.size()<<endl;
       	std::cout<<h.at(0)->GetEntries()<<std::endl;
}
