#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include <sstream>
#include <thread>
#include <iostream>

using namespace std;

#define DPM_DIR "FlatTree/triggerV1/SingleMuon/Run2016B_PromptReco_v2_MINIAOD/160622_212210/0000/"
#define NUMBER_OF_FILES 499
#define DATASET "SingleElectron"

//Skimming the files for selection
//met_pt>150 [line 28]


void SkimOneFile(int number){
	cout<<"@@ START @@ " <<endl;
	cout<<"taking care of file "<<number<<endl;
	string dpm_path = "root://sbgse1.in2p3.fr//cms/phedex/store/user/echabert/";
	stringstream fullfilename;
	fullfilename<<dpm_path<<DPM_DIR<<"output_"<<number<<".root";
	cout<<"filename: "<<fullfilename.str()<<endl;
	TFile* fin = TFile::Open(fullfilename.str().c_str());
	TTree* tree = (TTree*) fin->Get("FlatTree/tree");
	cout<<"Oringal tree have<<"<<tree->GetEntries()<<endl;
	TTree* selTree = tree->CopyTree("met_pt>150");	
	cout<<"Selected tree have "<<selTree->GetEntries()<<endl;
	stringstream ofilename;
	ofilename<<"DATASET"<<"_"<<number<<".root";
	TFile* output = new TFile(ofilename.str().c_str(),"RECREATE");
	selTree->Write();
	output->Write();
	output->Close();
	cout<<"@@ DONE @@"<<endl;
}

int main(int argc, char** argv){
	if(argc!=2){
		cerr<<"Need to give a number"<<endl;
		return -1;
	}
	int number = atoi(argv[1]);
	if(number>NUMBER_OF_FILES){
		cerr<<"This number exceed the maximum allowed !"<<endl;
		return -1;
	}
	cout<<"Run on file "<<number<<endl;
	
	gSystem->Load("/usr/lib64/libdpm.so");
	std::thread th(SkimOneFile,number);
	th.join();
	/*
	int nofloops = NUMBER_OF_FILES/nWorkers;
	std::vector<thread> workers;
	for(int i=0;i<NUMBER_OF_FILES;i++){
		workers.push_back(std::thread(SkimOneFile,i));
	}
	std::for_each(workers.begin(),workers.end(),[] (std::thread &t){t.join();});
	*/
	/*
	cout<<"@@@ LAST FILES"<<endl;
	// run on the last files [not //]
	// to be improved
	for(int i=nofloops*nWorkers;i<NUMBER_OF_FILES+1;i++){
		std::thread th(SkimOneFile,i);
		
	}
	*/
}
        /*
	int count150 = 0;
	for(unsigned int i=0;i<tree->GetEntries();i++){
		metBranch->GetEntry(i);
		//cout<<met<<endl;
		if(met>150){
			count50++;
			tree->GetEntry();
			
		}
		if(i%1000 == 0) cout<<count50*1./i<<endl;
	}
	*/



