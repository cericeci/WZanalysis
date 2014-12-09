
#include "WZEventMCOld.h"

#include "TLorentzVector.h"


#include <iostream>

//
//  Exteran function declarations
//

bool Z_muons(WZBASECLASS *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, float* pt, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, float* phi, float *eta);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);


WZEventMCOld::WZEventMCOld(TTree * tree) :
  WZBASECLASS(tree),
  //
  numZ(0),
  numW(0),
  numMET(0),
  num3e(0),
  num2e1mu(0),
  num1e2mu(0),
  num3mu(0),
  numMET3e(0),
  numMET2e1mu(0),
  numMET1e2mu(0),
  numMET3mu(0)
{
  pt.push_back(&pt1);
  pt.push_back(&pt2);
  pt.push_back(&pt3);
  pt.push_back(&pt4);
  //
  bdt.push_back(&bdt1);
  bdt.push_back(&bdt2);
  bdt.push_back(&bdt3);
  bdt.push_back(&bdt4);
  // 
  ch.push_back(&ch1);
  ch.push_back(&ch2);
  ch.push_back(&ch3);
  ch.push_back(&ch4);
  //
  eta.push_back(&eta1);
  eta.push_back(&eta2);
  eta.push_back(&eta3);
  eta.push_back(&eta4);
  // 
  phi.push_back(&phi1);
  phi.push_back(&phi2);
  phi.push_back(&phi3);
  phi.push_back(&phi4);
  // 
  /*  pdgid.push_back(&pdgid1);
  pdgid.push_back(&pdgid2);
  pdgid.push_back(&pdgid3);
  pdgid.push_back(&pdgid4);
  */
  // 
  pass2012ICHEP.push_back(&pass2012ICHEP1);
  pass2012ICHEP.push_back(&pass2012ICHEP2);
  pass2012ICHEP.push_back(&pass2012ICHEP3);
  pass2012ICHEP.push_back(&pass2012ICHEP4);
  //
  iso.push_back(&iso1);
  iso.push_back(&iso2);
  iso.push_back(&iso3);
  iso.push_back(&iso4);
  // 
  isomva.push_back(&isomva1);
  isomva.push_back(&isomva2);
  isomva.push_back(&isomva3);
  isomva.push_back(&isomva4);
}


void WZEventMCOld::ReadEvent()
{

  //  final_state     = undefined;
  //  selection_level = failsSelection;

  wLeptonIndex     = -1;
  zLeptonsIndex[0] = -1;
  zLeptonsIndex[1] = -1;


  // Red Reco leptons

  leptons.clear();

  float lepton_pt[4]  = {pt1,pt2,pt3,pt4};
  float lepton_eta[4] = {eta1,eta2,eta3,eta4};
  float lepton_phi[4] = {phi1,phi2,phi3,phi4};
  float lepton_ch[4]  = {ch1,ch2,ch3,ch4};
  //  float lepton_id[4]  = {pdgid1,pdgid2,pdgid3,pdgid4};
  float bdt[4] = {bdt1, bdt2, bdt3, bdt4};
  float lepton_id[4];
  
  for (int il=0; il<4; il++) {
    if (lepton_pt[il]>1.) {
      if (bdt[il]<100) lepton_id[il]=11.;
      else lepton_id[il]=13.;
      leptons.push_back(RecoLepton(lepton_pt[il],
				   lepton_eta[il],
				   lepton_phi[il],
				   lepton_ch[il],
				   lepton_id[il]));
    }

  }
  
  // Read recoJets  
  recoJets.clear();

  float recoJets_pt[6]  = {jetpt1,jetpt2,jetpt3,jetpt4,jetpt5,jetpt6};
  float recoJets_eta[6] = {jeteta1,jeteta2,jeteta3,jeteta4,jeteta5,jeteta6};
  float recoJets_phi[6] = {jetphi1,jetphi2,jetphi3,jetphi4,jetphi5,jetphi6};
  float recoJets_m[6]   = {jetmass1,jetmass2,jetmass3,jetmass4,jetmass5,jetmass6};

  TLorentzVector rj;  

  for (int i=0; i<6; i++) {
    if (recoJets_pt[i]>0.) {
      rj.SetPtEtaPhiM(recoJets_pt[i],recoJets_eta[i],recoJets_phi[i],recoJets_m[i]);
      recoJets.push_back(rj);
    }
  }
  //  std::cout<<"SIZE INSIDE:" <<recoJets.size()<<std::endl;


}
