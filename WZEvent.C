#include "WZEvent.h"

#include "TLorentzVector.h"


#include <iostream>

//
//  Exteran function declarations
//

bool Z_muons(WZBASECLASS *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, float* pt, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, float* phi, float *eta);


float RecoLepton::GetScaleFactor() {

  float factor = 1.;

  return factor;

}

WZEvent::WZEvent(TTree * tree) :
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
  pdgid.push_back(&pdgid1);
  pdgid.push_back(&pdgid2);
  pdgid.push_back(&pdgid3);
  pdgid.push_back(&pdgid4);
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

void WZEvent::ReadEvent()
{

  final_state     = undefined;
  selection_level = failsSelection;

  // Red Reco leptons

  leptons.clear();

  float lepton_pt[4]  = {pt1,pt2,pt3,pt4};
  float lepton_eta[4] = {eta1,eta2,eta3,eta4};
  float lepton_phi[4] = {phi1,phi2,phi3,phi4};
  float lepton_ch[4]  = {ch1,ch2,ch3,ch4};
  float lepton_id[4]  = {pdgid1,pdgid2,pdgid3,pdgid4};

  for (int il=0; il<4; il++) {
    if (lepton_pt[il]>1.) {
      leptons.push_back(RecoLepton(lepton_pt[il],
				   lepton_eta[il],
				   lepton_phi[il],
				   lepton_ch[il],
				   lepton_id[il]));
    }
  }


  // Read generated leptons

  genLeptons.clear();
  TLorentzVector gl1,gl2,gl3; // ,gl;
  //  gl1.SetPtEtaPhiM(leptonGenpt1,leptonGeneta1,leptonGenphi1,0.);
  //  gl2.SetPtEtaPhiM(leptonGenpt2,leptonGeneta2,leptonGenphi2,0.);
  //  gl3.SetPtEtaPhiM(leptonGenpt3,leptonGeneta3,leptonGenphi3,0.);

  float genLeptons_phi[3] = {genVV_S1lepton1_phi,genVV_S1lepton2_phi,genVV_S1lepton3_phi};
  float genLeptons_eta[3] = {genVV_S1lepton1_eta,genVV_S1lepton2_eta,genVV_S1lepton3_eta};
  float genLeptons_pt[3]  = {genVV_S1lepton1_pt,genVV_S1lepton2_pt,genVV_S1lepton3_pt};
  float genLeptons_pid[3]  = {genVV_S1lepton1_pid,genVV_S1lepton2_pid,genVV_S1lepton3_pid};
  float genLeptons_oVpid[3]  = {genVV_S1lepton1_oVpid,genVV_S1lepton2_oVpid,genVV_S1lepton3_oVpid};
  float genLeptons_imTau[3]  = {genVV_S1lepton1_imTau,genVV_S1lepton2_imTau,genVV_S1lepton3_imTau};

  //  gl1.SetPtEtaPhiM(genVV_S1lepton1_pt,genVV_S1lepton1_eta,genVV_S1lepton1_phi,0);
  //  gl2.SetPtEtaPhiM(genVV_S1lepton2_pt,genVV_S1lepton2_eta,genVV_S1lepton2_phi,0);
  //  gl3.SetPtEtaPhiM(genVV_S1lepton3_pt,genVV_S1lepton3_eta,genVV_S1lepton3_phi,0);
  //  genLeptons.push_back(gl1);
  //  genLeptons.push_back(gl2);
  //  genLeptons.push_back(gl3);

  for (int i=0; i<3; i++) {
    if (genLeptons_pt[i] > -1.) {
      GenS1Lepton gl(genLeptons_pt[i],genLeptons_eta[i],genLeptons_phi[i],
		     genLeptons_pid[i],genLeptons_oVpid[i],genLeptons_imTau[i]);

      //      gl.SetPtEtaPhiM(genLeptons_pt[i],genLeptons_phi[i],genLeptons_eta[i],0);
      genLeptons.push_back(gl);
    } 
  }


  for (int i=0; i<genLeptons.size(); i++) {
    float ptl = genLeptons[i].Pt();
    if (ptl>900.) {
      std::cout << "Gen lepton " << i << " \t Pt = " << ptl
		<< "\t channel: " << WZchan << std::endl;
    }
  }

  float genJets_phi[5] = {genVV_jet1_phi, genVV_jet2_phi, genVV_jet3_phi, genVV_jet4_phi, genVV_jet5_phi};
  float genJets_eta[5] = {genVV_jet1_eta, genVV_jet2_eta, genVV_jet3_eta, genVV_jet4_eta, genVV_jet5_eta};
  float genJets_pt[5]  = {genVV_jet1_pt, genVV_jet2_pt, genVV_jet3_pt, genVV_jet4_pt, genVV_jet5_pt};


  // Read generated jets
  TLorentzVector gj;  
  genJets.clear();
  for (int i=0; i<3; i++) {
    if (genJets_pt[i]>0.) {
      gj.SetPtEtaPhiM(genJets_pt[i],genJets_phi[i],genJets_eta[i],0);
      genJets.push_back(gj);
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
      rj.SetPtEtaPhiM(recoJets_pt[i],recoJets_phi[i],recoJets_eta[i],recoJets_m[i]);
      recoJets.push_back(rj);
    }
    // else {
    //      std::cout << "RECO JET WITH NEGATIVE PT: " << i << "\t"
    //		<< recoJets_pt[i] << " PHI = " 
    //		<< recoJets_phi[i] << std::endl;
    //    }
  }
  //  for (int i=0; i<genJets.size(); i++) {
  //    std::cout << "ZERO GEN JET: " << i << " " ;
  //    genJets[i].Print();
  //  }

}


void WZEvent::PrintSummary()
{

  std::cout<<" SUMMARY " << std::endl << std::endl;
  std::cout<<" ==================== \n \n";

}

bool WZEvent::passesSelection(){

  bool passed = false;


  selectedZPt=-88888.;
  //rejecting run 201191

  if (run==201191) return false;
  if (!trigger)  return false;

  const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);

  float pileUpWeight=this->puW;
    
  //rejecting run 201191
  if (this->run==201191) return passed;
  
  //  if (!(this->trigger)) return passed;
    
  float pts[leptonNumber]={pt1, pt2, pt3, pt4};
  float charges[leptonNumber]={ch1, ch2, ch3, ch4};
  float phis[leptonNumber]={phi1, phi2, phi3, phi4};
  float etas[leptonNumber]={eta1, eta2, eta3, eta4};
  //find Z boson, save the index of it
  std::vector<int> good_muons;
  std::vector<int> good_electrons;
  TLorentzVector v_nizEl[9];
  TLorentzVector v_nizMu[9];
  TLorentzVector v_3Lepton(0.,0.,0.,0.);
  
  int WZcandidates[3]; 
    
  //here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
  int lepNum(0);
  for (int i=0; i<leptonNumber; i++){
    if ((*pt[i]>10) && (*pt[i]!=-9999))
      lepNum++;
    //      if ((*bdt[i]>-10) && (*pass2012ICHEP[i]) && (*pt[i]>10)){
    if (( *pdgid[i]==-11 || *pdgid[i]==11)   && (*pt[i]>10)){
      good_electrons.push_back(i);
      v_nizEl[i].SetPtEtaPhiM(*pt[i],*eta[i], *phi[i], electronMass);
      v_3Lepton=v_3Lepton+v_nizEl[i];
    }
    //      if ((*bdt[i]<-10) && (*pass2012ICHEP[i]) && (*pt[i]>10)){
    if (( *pdgid[i]==-13 || *pdgid[i]==13)   && (*pt[i]>10)){
      good_muons.push_back(i);
      v_nizMu[i].SetPtEtaPhiM(*pt[i],*eta[i], *phi[i], muonMass);
      v_3Lepton=v_3Lepton+v_nizMu[i];
    }
  }
  
  if (lepNum!=3) return false;
  
  if (v_3Lepton.M()<100) return false;
  
  bool foundZel(false), foundZmu(false);
  double massMu(-999), massEl(0), Zpt(0);
  
  foundZmu=  Z_muons(this, &good_muons, WZcandidates, v_nizMu, pts, charges, massMu, Zpt);
  
  foundZel= Z_muons(this, &good_electrons, WZcandidates, v_nizEl, pts, charges, massEl, Zpt);
  
  // reject all events without Z boson--and candidates with double Z boson
  
  ///////////////////////////rejecting Zel&Zmu
  if ((foundZel) && (foundZmu)){
    return false;
  }
  
  /////////////////////////////rejecting two independent Zels
  if (Z_independent(charges, &good_muons, WZcandidates, v_nizMu)) 
    return false;
  
  
  ////////////////////////////rejecting two independent Zmus
  if (Z_independent(charges, &good_electrons, WZcandidates, v_nizEl)) 
    return false;
  
  ///////////////////////////rejecting when found no Zel& no Zmu
  if ((!foundZel) && (!foundZmu))
    return false;
  
  selection_level = passesZSelection;
  
  ///////////////////////FILLING HISTOGRAMS/////////////////////////////////////
  //later...
  
  //jet number...
  
  
  //     std::cout << "Run: " << run << "\t Event: " << event 
  // 	      << " Found Z: " << foundZel << foundZmu 
  // 	      << "\t Good lepton numbers: " 
  // 	      << good_electrons.size() 
  // 	      << "\t " << good_muons.size() << std::endl;
  
  
  //////////////////////FIND W BOSON///////////////////////////////////////////
  int ZmuWelCounter(0),ZmuWmuCounter(0), ZelWelCounter(0), ZelWmuCounter(0); 
  int ZmuWel_index(0), ZmuWmu_index(0), ZelWel_index(0), ZelWmu_index(0);
  bool testEl1(false),testMu1(false), testEl2(false), testMu2(false);
  bool ZmuWel(false), ZmuWmu(false), ZelWel(false),ZelWmu(false);
  bool ev3mu(false), ev1e2mu(false), ev2e1mu(false), ev3e(false);
  
  //*****Z->mumu*************
  if (foundZmu){
    //      std::cout << "Found Zmumu : " << ZmuWel << ZmuWmu << std::endl;
    for (int iel1=0; iel1< (good_electrons.size()); iel1++){
      int elIndex1= good_electrons[iel1];
      if (*pt[elIndex1]>20){
	ZmuWelCounter++;
	ZmuWel_index=elIndex1; 
	testEl1=true;  
      }
    }
    if (ZmuWelCounter==1) ZmuWel=true;
    int ZmuWmuCounter(0);
    
    for (int imu1=0; imu1< (good_muons.size());imu1++){
      int muIndex1=good_muons[imu1];
      if ((muIndex1==WZcandidates[0]) || (muIndex1==WZcandidates[1])) continue;
      
      if (*pt[muIndex1]>20){
	ZmuWmuCounter++;
	ZmuWmu_index=muIndex1;
	testMu1=true;
      }
    }
    if (ZmuWmuCounter==1) ZmuWmu=true;
    
    if ((ZmuWmu) && (ZmuWel))   return false;
    
    if ((!ZmuWmu) && (!ZmuWel)) return false;
    
    if ((ZmuWmu) && (!ZmuWel))
      {
	WZcandidates[2]=ZmuWmu_index;
	ev3mu=true;
      }
    
    if ((ZmuWel) && (!ZmuWmu))
      {
	WZcandidates[2]=ZmuWel_index;
	ev1e2mu=true;
      }
  }
  
  //**********Z->ee***************
  if (foundZel) {
    int ZelWelCouter(0);
    for (int iel2=0; iel2< (good_electrons.size()); iel2++){
      int elIndex2=good_electrons[iel2];
      if ((elIndex2==WZcandidates[0]) || (elIndex2==WZcandidates[1])) continue;
      if (*pt[elIndex2]>20)
	{
	  ZelWelCounter++;
	  ZelWel_index=elIndex2;
	  testEl2=true;
	}
    }
    if (ZelWelCounter==1)
      { 
	ZelWel=true;
      }
    int ZelWmuCounter(0);
    for (int imu2=0; imu2<(good_muons.size()); imu2++){
      int muIndex2=good_muons[imu2];
      if (*pt[muIndex2]>20)
	{
	  ZelWmuCounter++;
	  ZelWmu_index=muIndex2;
	  testMu2=true;
	}
    }
    if (ZelWmuCounter==1) ZelWmu=true;
    
    if ((ZelWel) && (ZelWmu)){
      return false;
    }
    if ((!ZelWel) && (!ZelWmu)) return false;
    
    if ((ZelWmu) && (!ZelWel))
      {
	WZcandidates[2]=ZelWmu_index;
	ev2e1mu=true;
      }
    if ((ZelWel) && (!ZelWmu))
      {
	WZcandidates[2]=ZelWel_index;
	ev3e=true;
      }
  }
  
  //    std::cout << ev3e << ev3mu << ev1e2mu << ev2e1mu << std::endl;    
  
  if ((!ev3e) && (!ev3mu) && (!ev1e2mu) && (!ev2e1mu)) return false;
  
  //deltaR condition
  if (!passDeltaRWleptZlept(WZcandidates, phis, etas)) return false;
  numW+=pileUpWeight;  
  
  selection_level = passesWSelection;
  
  if (ev3e){
    final_state = eee;
  }
  if (ev2e1mu){
    final_state = eem;
  }
  if (ev1e2mu){
    final_state = mme;
  }
  if (ev3mu){
    final_state = mmm;
  }
  //////////////////////////////////////MET CUT//////////////////////////
  
  if ((this->pfmet)<30)  return false;
  //    if ((this->pfmetTypeI)<30) continue;   ///CHANGE THIS
  
  selection_level = passesFullSelection;
  
  //    std::cout << "PASSED FINAL SELECTION: " << state << std::endl;


  TLorentzVector zl1,zl2,zp4;;
  

  zl1.SetPtEtaPhiM(*pt[WZcandidates[0]],*eta[WZcandidates[0]],
		   *phi[WZcandidates[0]], muonMass);

  zl2.SetPtEtaPhiM(*pt[WZcandidates[1]],*eta[WZcandidates[1]],
		   *phi[WZcandidates[1]], muonMass);

  zp4 = zl1+zl2;

  selectedZPt = zp4.Pt();
  
  return true;
  
}

bool WZEvent::PassesGenCuts(){

  return (MZ>71. && MZ<111.);

}
