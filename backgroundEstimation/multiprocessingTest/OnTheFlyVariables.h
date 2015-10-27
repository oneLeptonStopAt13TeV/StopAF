#ifndef ON_THE_FLY_VARIABLES
#define ON_THE_FLY_VARIABLES

struct OnTheFlyVariables
{
public:
   OnTheFlyVariables()
    : m_genWPt(0),
      m_genTopPt(0),
      m_genWdR(0),
      m_genTopdR(0),
      m_mumIdx(-1)
     {
     }

public:
    float someCoolVariable;

    int SRbins;

public:
    float m_genWPt;

    float m_genTopPt;

    float m_genWdR;

    float m_genTopdR;

    float m_recWPt;

    float m_recTopPt;

    int m_mumIdx;
    
    vector<int> m_daughtersIdx;

public:
    float m_ak8WMass;
    float m_ak8WPrunedMass;
    float m_ak8WTrimmedMass;
    float m_ak8WSoftDropMass;
    float m_ak8recoWPt;

    float m_ak8WFakeMass;
    float m_ak8WFakePrunedMass;
    float m_ak8WFakeTrimmedMass;
    float m_ak8WFakeSoftDropMass;
    float m_ak8recoWFakePt;

};

OnTheFlyVariables onTheFlyVariables;

void ComputeOnTheFlyVariables()
{
    onTheFlyVariables.someCoolVariable = 3.14;


    //-- Compute SRbin variables (MT2W/MET binning)
    onTheFlyVariables.SRbins = 0;

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

//find mother particle with given index
//looking for particle which is decaying to two particles, which are not lepton
void findMotherAndDaughters(int idx)
{
    //initial clean-up
    onTheFlyVariables.m_mumIdx = -1;
    
    //loop over particles and find mum candidates
    for(uint32_t i = 0; i < myEvent.gen_id.size(); i++)
    {

       if((onTheFlyVariables.m_mumIdx == -1) && (abs(myEvent.gen_id.at(i)) == idx) && myEvent.gen_daughter_index.at(i).size() == 2) // find W or t = mother particles
       {
           int dp1  = myEvent.gen_daughter_index.at(i).at(0);
           int dp2  = myEvent.gen_daughter_index.at(i).at(1);
           if( !( (11 <= myEvent.gen_id.at(dp1) && myEvent.gen_id.at(dp1) <= 18) || (11 <= myEvent.gen_id.at(dp2) && myEvent.gen_id.at(dp2)<= 18) ))
           {	
                onTheFlyVariables.m_mumIdx = myEvent.gen_index.at(i);
                onTheFlyVariables.m_daughtersIdx.push_back(dp1);
                onTheFlyVariables.m_daughtersIdx.push_back(dp2);  
                //@MJ@ TODO we need t and bar t!!!
           }
       }
       
       //W decaying to two quarks found
       if(idx == 24 && onTheFlyVariables.m_mumIdx != -1 && abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0))) <= 6 && abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1))) <= 6)
       {
               break;
       }
       
      //t decaying to W and sth found 
       else if (idx == 6 && onTheFlyVariables.m_mumIdx != -1 && (abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0))) == 24 || abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1))) == 24))
       {
               break;
       }
       //this was not correct mother particle
       else
       {
	       onTheFlyVariables.m_daughtersIdx.clear();
	       onTheFlyVariables.m_mumIdx = -1;
       }

    }
}

void findGenParticleProps(int idx, float* genPt, float* gendR)
{ 
    //do initial clean up
    onTheFlyVariables.m_mumIdx = -1;
  
    onTheFlyVariables.m_daughtersIdx.clear();
    
    //loop over the generated paricles, find mother and its daughters
    findMotherAndDaughters(idx);

    //no mother found
    if(onTheFlyVariables.m_mumIdx == -1)
    {
        //cout << "no mother particle was found! " << "mum Id is: " << onTheFlyVariables.m_mumIdx << endl;
    }

    //fill the kinematic variables
    if(onTheFlyVariables.m_mumIdx != -1)
    {
        *genPt = myEvent.gen_pt.at(onTheFlyVariables.m_mumIdx);
    
        float deta = myEvent.gen_eta.at(onTheFlyVariables.m_daughtersIdx.at(0)) - myEvent.gen_eta.at(onTheFlyVariables.m_daughtersIdx.at(1));
        float dphi = myEvent.gen_phi.at(onTheFlyVariables.m_daughtersIdx.at(0)) - myEvent.gen_phi.at(onTheFlyVariables.m_daughtersIdx.at(1));
        *gendR =  TMath::Sqrt( deta*deta+dphi*dphi );
    }
    else
    {
        *genPt = -1;
        *gendR = -1;
    }
    
    //final cleanup
    onTheFlyVariables.m_daughtersIdx.clear();
    onTheFlyVariables.m_mumIdx = -1; 
}

//count efficiency of particle tagging
void countEfficiency(int idx)
{
  //initial clean up
  onTheFlyVariables.m_daughtersIdx.clear();
  onTheFlyVariables.m_ak8WMass = -1;
  onTheFlyVariables.m_ak8WPrunedMass = -1;
  onTheFlyVariables.m_ak8WTrimmedMass = -1;
  onTheFlyVariables.m_ak8WSoftDropMass = -1;
  onTheFlyVariables.m_ak8recoWPt = -1;


  onTheFlyVariables.m_ak8WFakeMass = -1;
  onTheFlyVariables.m_ak8WFakePrunedMass = -1;
  onTheFlyVariables.m_ak8WFakeTrimmedMass = -1;
  onTheFlyVariables.m_ak8WFakeSoftDropMass = -1;
  onTheFlyVariables.m_ak8recoWFakePt = -1;
  
  uint32_t nrOfJets = myEvent.ak8pfjets_mass.size();
  
  for(uint32_t i = 0; i < nrOfJets; i++ )
  {
       //@MJ@ TODO add another info about W tagged jet in order to discriminate fakes rate
       // e.g. sufficient mass of jet, number of subjets
       // find mother particle
       findMotherAndDaughters(idx);
       float deta = -13;
       float dphi = -13;
       float dR = -13;
       
       //count dR
       if(onTheFlyVariables.m_mumIdx != -1)
       {
           deta = myEvent.ak8pfjets_eta.at(i) - myEvent.gen_eta.at(onTheFlyVariables.m_mumIdx);
           dphi = myEvent.ak8pfjets_phi.at(i) - myEvent.gen_phi.at(onTheFlyVariables.m_mumIdx);
           dR = TMath::Sqrt( deta*deta+dphi*dphi );
       }
       //Ws
       if(abs(dR) < 0.1 && myEvent.ak8pfjets_pt.at(i) > 100) //@MJ@ TODO put the correct number
       {
           //@MJ@ TODO count subjets?
           onTheFlyVariables.m_ak8WMass = myEvent.ak8pfjets_mass.at(i);
           onTheFlyVariables.m_ak8WPrunedMass = myEvent.ak8pfjets_pruned_mass.at(i);
           onTheFlyVariables.m_ak8WTrimmedMass = myEvent.ak8pfjets_trimmed_mass.at(i);
           onTheFlyVariables.m_ak8WSoftDropMass = myEvent.ak8pfjets_softdrop_mass.at(i);
           onTheFlyVariables.m_ak8recoWPt = myEvent.ak8pfjets_pt.at(i);
           
           onTheFlyVariables.m_ak8WFakeMass = -1;
           onTheFlyVariables.m_ak8WFakePrunedMass = -1;
           onTheFlyVariables.m_ak8WFakeTrimmedMass = -1;
           onTheFlyVariables.m_ak8WFakeSoftDropMass = -1;
           onTheFlyVariables.m_ak8recoWFakePt = -1;
           break;
       }
       //fakes
       else if(myEvent.ak8pfjets_pt.at(i) > 100)
       {
           onTheFlyVariables.m_ak8WFakeMass = myEvent.ak8pfjets_mass.at(i);
           onTheFlyVariables.m_ak8WFakePrunedMass = myEvent.ak8pfjets_pruned_mass.at(i);
           onTheFlyVariables.m_ak8WFakeTrimmedMass = myEvent.ak8pfjets_trimmed_mass.at(i);
           onTheFlyVariables.m_ak8WFakeSoftDropMass = myEvent.ak8pfjets_softdrop_mass.at(i);
           onTheFlyVariables.m_ak8recoWFakePt = myEvent.ak8pfjets_pt.at(i);
           
           onTheFlyVariables.m_ak8WMass = -1;
           onTheFlyVariables.m_ak8WPrunedMass = -1;
           onTheFlyVariables.m_ak8WTrimmedMass = -1;
           onTheFlyVariables.m_ak8WSoftDropMass = -1;
           onTheFlyVariables.m_ak8recoWPt = -1;
       }
       //nothing found
       else
       {
           cout << "W boosted jet not found" << endl;
       }
  }  
}

#endif
