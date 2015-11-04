#ifndef ON_THE_FLY_VARIABLES
#define ON_THE_FLY_VARIABLES

struct OnTheFlyVariables
{
public:
   OnTheFlyVariables()
    : m_genWPt(0),
      m_genTopPt(0),
      m_genBarTopPt(0),
      m_genWdR(0),
      m_genTopdR(0),
      m_genBarTopdR(0),
      m_mumIdx(-1)
     {
     }

public:
    float someCoolVariable;

    int SRbins;

public:
    float m_genWPt;
    float m_genTopPt;
    float m_genBarTopPt;
    
    float m_genWdR;
    float m_genTopdR;
    float m_genBarTopdR;

    float m_recWPt;

    //float m_recTopPt; //@MJ@ TODO not used

    int m_mumIdx;
    
    vector<int> m_daughtersIdx;

public:
    float m_ak8WMass;
    float m_ak8WPrunedMass;
    float m_ak8WTrimmedMass;
    float m_ak8WSoftDropMass;
    float m_ak8recoWPt;
    float m_ak8WSubJets;
    float m_ak8MET;

    float m_ak8WFakeMass;
    float m_ak8WFakePrunedMass;
    float m_ak8WFakeTrimmedMass;
    float m_ak8WFakeSoftDropMass;
    float m_ak8recoWFakePt;
    float m_ak8WFakeSubJets;
    float m_ak8FakeMET;

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

       if((onTheFlyVariables.m_mumIdx == -1) && (abs(myEvent.gen_id.at(i)) == abs(idx)) && myEvent.gen_daughter_index.at(i).size() == 2) // find W or t = mother particles
       {
           int dp1  = myEvent.gen_daughter_index.at(i).at(0);
           int dp2  = myEvent.gen_daughter_index.at(i).at(1);
           if( !( (11 <= myEvent.gen_id.at(dp1) && myEvent.gen_id.at(dp1) <= 18) || (11 <= myEvent.gen_id.at(dp2) && myEvent.gen_id.at(dp2)<= 18) ))
           {	
                onTheFlyVariables.m_mumIdx = myEvent.gen_index.at(i);
                onTheFlyVariables.m_daughtersIdx.push_back(dp1);
                onTheFlyVariables.m_daughtersIdx.push_back(dp2);  
           }
       }

       if(onTheFlyVariables.m_mumIdx == -1) 
           continue;
       if(abs(myEvent.gen_id.at(onTheFlyVariables.m_mumIdx)) != abs(idx))
       {
	   onTheFlyVariables.m_daughtersIdx.clear();
	   onTheFlyVariables.m_mumIdx = -1;
           continue;
       }

       //W decaying to two quarks found
       if(idx == 24  && abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0))) <= 6 && abs(myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1))) <= 6)
       {
               //cout << "W found" << endl;
               break;
       }
       
      //t decaying to W+ and sth found 
       else if (idx == 6 && (myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0)) == 24 || myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1)) == 24))
       {
               //cout << "top found" << endl;
               break;
       }
       
      //bar t decaying to W- and sth found 
       else if (idx == -6 && (myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(0)) == -24 || myEvent.gen_id.at(onTheFlyVariables.m_daughtersIdx.at(1)) == -24))
       {
               //cout << "bar top found" << endl;
               break;
       }
       //this was not correct mother particle
       else
       {
               //cout << "mother condition not fulfilled" << endl;
	       onTheFlyVariables.m_daughtersIdx.clear();
	       onTheFlyVariables.m_mumIdx = -1;
       }

    }
}

      //@MJ@ TODO what about ak10 jets?
      void fillRealAndFakes(float deltaR, uint32_t i)
      {
           if(abs(deltaR) < 0.1)
           {
                onTheFlyVariables.m_recWPt = myEvent.ak8pfjets_pt.at(i);
           }
    
           //Ws
           if(abs(deltaR) < 0.1 && myEvent.ak8pfjets_pt.at(i) > 100) //@MJ@ TODO put the correct number
           {
               //@MJ@ TODO count subjets?
               onTheFlyVariables.m_ak8WMass = myEvent.ak8pfjets_mass.at(i);
               onTheFlyVariables.m_ak8WPrunedMass = myEvent.ak8pfjets_pruned_mass.at(i);
               onTheFlyVariables.m_ak8WTrimmedMass = myEvent.ak8pfjets_trimmed_mass.at(i);
               onTheFlyVariables.m_ak8WSoftDropMass = myEvent.ak8pfjets_softdrop_mass.at(i);
               onTheFlyVariables.m_ak8recoWPt = myEvent.ak8pfjets_pt.at(i);
               onTheFlyVariables.m_ak8WSubJets = myEvent.ak8pfjets_nSubJets.at(i); //for the event
               onTheFlyVariables.m_ak8MET = myEvent.pfmet; //for the event
           
               onTheFlyVariables.m_ak8WFakeMass = -1;
               onTheFlyVariables.m_ak8WFakePrunedMass = -1;
               onTheFlyVariables.m_ak8WFakeTrimmedMass = -1;
               onTheFlyVariables.m_ak8WFakeSoftDropMass = -1;
               onTheFlyVariables.m_ak8recoWFakePt = -1;
               onTheFlyVariables.m_ak8WFakeSubJets = -1; //for the event
               onTheFlyVariables.m_ak8FakeMET = -1; //for the event
           }
           //fakes
           else if(abs(deltaR) > 0.2 && myEvent.ak8pfjets_pt.at(i) > 100)
           {
               onTheFlyVariables.m_ak8WFakeMass = myEvent.ak8pfjets_mass.at(i);
               onTheFlyVariables.m_ak8WFakePrunedMass = myEvent.ak8pfjets_pruned_mass.at(i);
               onTheFlyVariables.m_ak8WFakeTrimmedMass = myEvent.ak8pfjets_trimmed_mass.at(i);
               onTheFlyVariables.m_ak8WFakeSoftDropMass = myEvent.ak8pfjets_softdrop_mass.at(i);
               onTheFlyVariables.m_ak8recoWFakePt = myEvent.ak8pfjets_pt.at(i);
               onTheFlyVariables.m_ak8WFakeSubJets = myEvent.ak8pfjets_nSubJets.at(i); //for the event
               onTheFlyVariables.m_ak8FakeMET = myEvent.pfmet; //for the event
           
               onTheFlyVariables.m_ak8WMass = -1;
               onTheFlyVariables.m_ak8WPrunedMass = -1;
               onTheFlyVariables.m_ak8WTrimmedMass = -1;
               onTheFlyVariables.m_ak8WSoftDropMass = -1;
               onTheFlyVariables.m_ak8recoWPt = -1;
               onTheFlyVariables.m_ak8WSubJets = -1; //for the event
               onTheFlyVariables.m_ak8MET = -1; //for the event
           }
           //nothing found
           else
           {
               //cout << "W boosted jet not found" << endl;
           }
       }

void findGenParticleProps(int idx, float* genPt, float* gendR)
{ 
    //do initial clean up
    onTheFlyVariables.m_mumIdx = -1;
    *genPt = -1;
    *gendR = -1;
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
void countEfficiency(int idx, string cutName = "", float lowCut = -1, float upCut = -1)
{
  //initial clean up
  onTheFlyVariables.m_daughtersIdx.clear();
  onTheFlyVariables.m_ak8WMass = -1;
  onTheFlyVariables.m_ak8WPrunedMass = -1;
  onTheFlyVariables.m_ak8WTrimmedMass = -1;
  onTheFlyVariables.m_ak8WSoftDropMass = -1;
  onTheFlyVariables.m_ak8recoWPt = -1;
  onTheFlyVariables.m_ak8WSubJets = -1; //for the event
  onTheFlyVariables.m_ak8MET = -1;
  onTheFlyVariables.m_recWPt = -1;
  
  onTheFlyVariables.m_ak8WFakeMass = -1;
  onTheFlyVariables.m_ak8WFakePrunedMass = -1;
  onTheFlyVariables.m_ak8WFakeTrimmedMass = -1;
  onTheFlyVariables.m_ak8WFakeSoftDropMass = -1;
  onTheFlyVariables.m_ak8recoWFakePt = -1;
  onTheFlyVariables.m_ak8WFakeSubJets = -1; //for the event
  onTheFlyVariables.m_ak8FakeMET = -1;
  
   //choose on which variable cut
   vector<float> masses(myEvent.ak8pfjets_mass.size(), -1.0);
   if(cutName == "raw")
       masses = myEvent.ak8pfjets_mass;
   else if(cutName == "pruned")
       masses = myEvent.ak8pfjets_pruned_mass;
   else if(cutName == "trimmed")
       masses = myEvent.ak8pfjets_trimmed_mass;
   else if(cutName == "softdrop")
       masses = myEvent.ak8pfjets_softdrop_mass;
   else
   {
       //cout << "no cut" << endl;
   }

  uint32_t nrOfJets = myEvent.ak8pfjets_mass.size();
  
  for(uint32_t i = 0; i < nrOfJets; i++ )
  {

       //@MJ@ TODO add another info about W tagged jet in order to discriminate fakes rate
       // e.g. sufficient mass of jet, number of subjets, N-subjettiness, PUPPI
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
       else
          continue;

       //perform cuts and filling fakes/reals
       if(cutName == "")
       {
           fillRealAndFakes(dR, i);
           if(onTheFlyVariables.m_ak8WMass == myEvent.ak8pfjets_mass.at(i)) //end the loop when W jet found //@MJ@ TODO find better solution
               break;
       }
       else if(cutName != "" && lowCut != -1 && upCut != -1)
       {
           if(masses.at(i) > lowCut && masses.at(i) < lowCut)
               fillRealAndFakes(dR, i);
           if(onTheFlyVariables.m_ak8WMass == myEvent.ak8pfjets_mass.at(i)) //end the loop when W jet found
               break;
       }
       else if(cutName != "" && (lowCut != -1 || upCut != -1))
       {
           if(lowCut != -1 && masses.at(i) > lowCut)
               fillRealAndFakes(dR, i);
           if(upCut != -1 && masses.at(i) < upCut)
               fillRealAndFakes(dR, i);
           if(onTheFlyVariables.m_ak8WMass == myEvent.ak8pfjets_mass.at(i))  //end the loop when W jet found
               break;
       }
      else
      {
          cout << "cuts were defined in wrog way!" << endl;
      }
  }  
}


#endif
