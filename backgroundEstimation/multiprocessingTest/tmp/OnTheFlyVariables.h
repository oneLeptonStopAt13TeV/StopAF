#ifndef ON_THE_FLY_VARIABLES
#define ON_THE_FLY_VARIABLES

struct OnTheFlyVariables
{
public:
   OnTheFlyVariables()
    :m_mumIdx(-1),
     m_nrOfGenW(0),
     m_nrOfGenTop(0)
    {
    }

public:
    float someCoolVariable;

    int SRbins;

    float m_genWPt;

    float m_genTopPt;

    float m_genWdR;

    float m_genTopdR;

    float m_recWPt;

    float m_recTopPt;

    int m_mumIdx;

    int m_nrOfGenW;

    int m_nrOfGenTop;
};

OnTheFlyVariables onTheFlyVariables;

void ComputeOnTheFlyVariables()
{
    onTheFlyVariables.someCoolVariable = 3.14;


    //-- Compute SRbin variables (MT2W/MET binning)
    onTheFlyVariables.SRbins = 0;

    cout << "log1" << myEvent.MT2W << endl;
    if(myEvent.MT2W<200 && myEvent.MT2W>=0){
    	if(myEvent.pfmet>=250 && myEvent.pfmet<300) onTheFlyVariables.SRbins = -1;
	if(myEvent.pfmet>=300 && myEvent.pfmet<350) onTheFlyVariables.SRbins = -2;
    	if(myEvent.pfmet>=350 && myEvent.pfmet<400) onTheFlyVariables.SRbins = -3;
    	if(myEvent.pfmet>=400) onTheFlyVariables.SRbins = -4;
    }
    if(myEvent.MT2W>=200){
    	if(myEvent.pfmet>=250 && myEvent.pfmet<300) onTheFlyVariables.SRbins = 1;
	if(myEvent.pfmet>=300 && myEvent.pfmet<350) onTheFlyVariables.SRbins = 2;
    	if(myEvent.pfmet>=350 && myEvent.pfmet<400) onTheFlyVariables.SRbins = 3;
    	if(myEvent.pfmet>=400 && myEvent.pfmet<450) onTheFlyVariables.SRbins = 4;
    	if(myEvent.pfmet>=450 ) onTheFlyVariables.SRbins = 5;
    }
    //-- End SRbin computation

}

void findGenParticleProps(int idx, float* genPt, float* gendR)
{
    cout << "log2" << myEvent.MT2W << endl;
    
    //cout << "moth " << partId << "daugh " << daughterId1 << endl;
    vector<int> partId = myEvent.gen_id;
    vector<int> quarksIdx;
    vector<int> wIdx;
    cout << "particle id size: " << partId.size() << " from tree " << myEvent.gen_id.size() << endl; 
    
    for(uint32_t i = 0; i< partId.size(); i++)
    {
       cout << "mother id: " << partId.at(i) << endl;

       if((onTheFlyVariables.m_mumIdx == -1) && (abs(partId.at(i)) == idx)) // find W or t = mother particles
       {
         onTheFlyVariables.m_mumIdx = myEvent.gen_index.at(i);
         cout << "mother index: " << onTheFlyVariables.m_mumIdx << endl;
         //@MJ@ TODO there should be more of them - status?
       }

       // find daughter particles
       if(abs(partId.at(i)) <= 6) 
         quarksIdx.push_back(myEvent.gen_index.at(i));
       if(abs(partId.at(i)) == 24) 
         wIdx.push_back(myEvent.gen_index.at(i));
    }
    cout << "mother id reading finished: " << endl;
    
    vector<int> daughtersIdx;

    //match the generated daughters with mother 
    if(idx == 24 && onTheFlyVariables.m_mumIdx != -1)
    {
       for(uint32_t j = 0; j<quarksIdx.size(); j++)
       {
           if(myEvent.gen_mother_index.at(quarksIdx.at(j)) == onTheFlyVariables.m_mumIdx)
               daughtersIdx.push_back(quarksIdx.at(j));
       }
       cout << "size of daughters is: " << daughtersIdx.size() << endl;
    } 
    else if (idx == 6 && onTheFlyVariables.m_mumIdx != -1) //what or bar t???
    {
       for(uint32_t k = 0; k<quarksIdx.size(); k++)
       {
           if(myEvent.gen_mother_index.at(quarksIdx.at(k)) == onTheFlyVariables.m_mumIdx && abs(myEvent.gen_id.at(quarksIdx.at(k))) == 5)
               daughtersIdx.push_back(quarksIdx.at(k));
       }
       for(uint32_t l = 0; l<wIdx.size(); l++)
       {
           if(myEvent.gen_mother_index.at(wIdx.at(l)) == onTheFlyVariables.m_mumIdx)
               daughtersIdx.push_back(wIdx.at(w));
       }
       cout << "size of daughters is: " << daughtersIdx.size() << endl;
    }
    else
    {
        cout << "no mother particle was found! " << "mum Id is: " << onTheFlyVariables.m_mumIdx << endl;
    }

    if(daughtersIdx.size() == 2)
    {
        //@MJ@ TODO do I need this?
        //if(abs(myEvent.gen_id.at(onTheFlyVariables.m_mumIdx)) == 24)
         //   onTheFlyVariables.m_nrOfGenW++;
          //  onTheFlyVariables.m_nrOfGenTop++;
    *genPt = myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx);
    
    float deta = myEvent.gen_eta.at(daughtersIdx.at(1)) - myEvent.gen_eta.at(daughtersIdx.at(2));
    float dphi = myEvent.gen_phi.at(daughtersIdx.at(1)) - myEvent.gen_phi.at(daughtersIdx.at(2));;
    *gendR =  TMath::Sqrt( deta*deta+dphi*dphi );
    }
    else
    {
        cout << "number of generated particles from W/top is not two, it is: " << daughtersIdx.size();
    }
}

void findRecoParticleProps()
{
    //then obatin pt of reconstructed pt/W
    onTheFlyVariables.m_recWPt = 300;
    onTheFlyVariables.m_recTopPt = 500;    
}

#endif
