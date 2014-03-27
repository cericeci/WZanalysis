#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

// Replace this with the new tree
//#include "WZ.h"
#include "WZ2012Data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>


//function declaration goes here:
//
//
//
bool Z_muons(WZ2012Data *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, float* pt, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, float* phi, float *eta);

int determineLabel(int * pass2012ICHEP, int * WZcandidates);

double weight(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);

float GetError(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

float MMerror(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label, float weight);

//float MMerror(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label);

int main()
{
  using namespace std;
  
  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;

  myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3e_Lucija.txt");
  myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts2e1mu_Lucija.txt");
  myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts1e2mu_Lucija.txt");
  myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3mu_Lucija.txt");
  myfileAll.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/all_Lucija.txt");


  TFile * fout= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data.root", "RECREATE");

  TH1F * hZmassMu1         = new TH1F ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
  TH1F * hZmassEl1         = new TH1F ("hZmassEl1", "hZmassEl1", 100, 60, 120);  
  TH2F* MuonFR;
  TH2F* MuonPR;
  TH2F* ElectronFR;
  TH2F* ElectronPR;
  
  MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet20_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  MuonPR=LoadHistogram("auxiliaryFiles/MuPR_Moriond13_2012.root", "h2inverted", "MuonPR");
  ElectronFR=LoadHistogram("auxiliaryFiles/EleFR_Moriond13_jet35_EWKcorr.root", "fakeElH2", "ElectronFR");
  ElectronPR=LoadHistogram("auxiliaryFiles/ElePR_Moriond13_2012.root", "h2inverted", "ElectronPR");
    
  const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);
  int numZ(0), numW(0), numMET(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0), numMET2e1mu(0), numMET1e2mu(0), numMET3mu(0);
  float N_good_3e(0), N_good_2e1mu(0), N_good_1e2mu(0), N_good_3mu(0);
  float N_fake_3e(0), N_fake_2e1mu(0), N_fake_1e2mu(0), N_fake_3mu(0);
  float ferror_3e(0),ferror_2e1mu(0), ferror_1e2mu(0), ferror_3mu(0);
  float error_3e(0),error_2e1mu(0), error_1e2mu(0), error_3mu(0);

  //test variables
  int numZ1(0), numZ2(0), numZ3(0), numZ4(0), numZ5(0), numZ6(0);
  int electronCount(0);
  int justCount(0);
  
  std::vector<TString>inputName;
  TChain wz("latino");
  
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_103_DoubleMuon2012A.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_104_MuEG2012A.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_112_DoubleElectron2012A.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_102_DoubleElectron2012A.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_113_DoubleMuon2012A.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_114_MuEG2012A.root");

    
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_202_DoubleElectron2012B.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_203_DoubleMuon2012B.root ");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_204_MuEG2012B.root");
  
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_308_DoubleMuon2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_309_MuEG2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_312_DoubleElectron2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_313_DoubleMuon2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_314_MuEG2012C.root");
  
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_322_DoubleElectron2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_323_DoubleMuon2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_324_MuEG2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_332_DoubleElectron2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_333_DoubleMuon2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_334_MuEG2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_342_DoubleElectron2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_343_DoubleMuon2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_344_MuEG2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_352_DoubleElectron2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_353_DoubleMuon2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_354_MuEG2012C.root");
  

  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_402_DoubleElectron2012D.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_307_DoubleElectron2012C.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_403_DoubleMuon2012D.root");
  inputName.push_back("/STORE/lucija/latinosTrees/Data_LooseLooseTypeI/latino_404_MuEG2012D.root");
  

  
  for (int input=0; input< inputName.size(); input++){
    wz.Add(inputName[input]);
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZ2012Data *cWZ= new WZ2012Data(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;


  for  (Int_t k = 0; k<events /*&& k<10000*/;k++) {
    wz_tTree->GetEntry(k);

    //rejecting run 201191
    if (cWZ->run==201191) continue;
    if (!(cWZ->trigger)) continue;

    
    float pt[leptonNumber]={cWZ->pt1, cWZ->pt2, cWZ->pt3, cWZ->pt4};
    float bdt[leptonNumber]={cWZ->bdt1, cWZ->bdt2, cWZ->bdt3, cWZ->bdt4};
    float ch[leptonNumber]={cWZ->ch1, cWZ->ch2, cWZ->ch3, cWZ->ch4};
    float eta[leptonNumber]={cWZ->eta1, cWZ->eta2, cWZ->eta3, cWZ->eta4};
    float phi[leptonNumber]={cWZ->phi1, cWZ->phi2, cWZ->phi3, cWZ->phi4};
    int pass2012ICHEP[leptonNumber]={cWZ->pass2012ICHEP1, cWZ->pass2012ICHEP2, cWZ->pass2012ICHEP3, cWZ->pass2012ICHEP4};
    float iso[leptonNumber]={cWZ->iso1, cWZ->iso2, cWZ->iso3, cWZ->iso4};
    float isomva[leptonNumber]={cWZ->isomva1, cWZ->isomva2, cWZ->isomva3, cWZ->isomva4};
    //find Z boson, save the index of it
    std::vector<int> good_muons;
    std::vector<int> good_electrons;
    TLorentzVector v_nizEl[9];
    TLorentzVector v_nizMu[9];
    TLorentzVector v_3Lepton(0.,0.,0.,0.);
    
    int WZcandidates[3]; //->here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
    int lepNum(0);
    for (int i=0; i<leptonNumber; i++){
      if ((pt[i]>10) && (pt[i]!=-9999))
	lepNum++;
      if ((bdt[i]<100) && (pt[i]>10)){
	good_electrons.push_back(i);
	v_nizEl[i].SetPtEtaPhiM(pt[i],eta[i], phi[i], electronMass);
	v_3Lepton=v_3Lepton+v_nizEl[i];
      }
      if ((bdt[i]>100) && (pt[i]>10)){
	good_muons.push_back(i);
	v_nizMu[i].SetPtEtaPhiM(pt[i],eta[i], phi[i], muonMass);
	v_3Lepton=v_3Lepton+v_nizMu[i];
      }
    }
    if (lepNum!=3) continue; ///???

    if (v_3Lepton.M()<100) continue;
    
    justCount++;
    bool foundZel(false), foundZmu(false);
    double massMu(-999), massEl(0), Zpt(0);

    foundZmu=  Z_muons(cWZ, &good_muons, WZcandidates, v_nizMu, pt, ch, massMu, Zpt);
    if (foundZmu){
      hZmassMu1->Fill(massMu);
    }

    foundZel= Z_muons(cWZ, &good_electrons, WZcandidates, v_nizEl, pt, ch, massEl, Zpt);

    // reject all events without Z boson--and candidates with double Z boson
        
    ///////////////////////////rejecting Zel&Zmu
    if ((foundZel) && (foundZmu)){
      continue;
    }
       
    /////////////////////////////rejecting two independent Zels
    if (Z_independent(ch, &good_muons, WZcandidates, v_nizMu)) 
      continue;

        
    ////////////////////////////rejecting two independent Zmus
    if (Z_independent(ch, &good_electrons, WZcandidates, v_nizEl)) 
      continue;
    
        
    ///////////////////////////rejecting when found no Zel& no Zmu
    if ((!foundZel) && (!foundZmu))
      continue;

    numZ++;

    ///////////////////////FILLING HISTOGRAMS/////////////////////////////////////
    //later...
 
    //jet number...

    
    //////////////////////FIND W BOSON///////////////////////////////////////////
    int ZmuWelCounter(0),ZmuWmuCounter(0), ZelWelCounter(0), ZelWmuCounter(0); 
    int ZmuWel_index(0), ZmuWmu_index(0), ZelWel_index(0), ZelWmu_index(0);
    bool testEl1(false),testMu1(false), testEl2(false), testMu2(false);
    bool ZmuWel(false), ZmuWmu(false), ZelWel(false),ZelWmu(false);
    bool ev3mu(false), ev1e2mu(false), ev2e1mu(false), ev3e(false);
      
    //*****Z->mumu*************
    if (foundZmu){
      for (int iel1=0; iel1< (good_electrons.size()); iel1++){
	int elIndex1= good_electrons[iel1];
	if (pt[elIndex1]>20){
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
	if (pt[muIndex1]>20){
	  ZmuWmuCounter++;
	  ZmuWmu_index=muIndex1;
	  testMu1=true;
	}
      }
      if (ZmuWmuCounter==1) ZmuWmu=true;

      if ((ZmuWmu) && (ZmuWel)) continue;  
      
      if ((!ZmuWmu) && (!ZmuWel)) continue;
      
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
    if (foundZel){
      int ZelWelCouter(0);
      for (int iel2=0; iel2< (good_electrons.size()); iel2++){
	int elIndex2=good_electrons[iel2];
	if ((elIndex2==WZcandidates[0]) || (elIndex2==WZcandidates[1])) continue;
	if (pt[elIndex2]>20)
	  {
	    ZelWelCounter++;
	    ZelWel_index=elIndex2;
	    testEl2=true;
	  }
      }
      if (ZelWelCounter==1)
	{ 
	  ZelWel=true;
	  electronCount++;
	}
      int ZelWmuCounter(0);
      for (int imu2=0; imu2<(good_muons.size()); imu2++){
	int muIndex2=good_muons[imu2];
	if (pt[muIndex2]>20)
	  {
	    ZelWmuCounter++;
	    ZelWmu_index=muIndex2;
	    testMu2=true;
	  }
      }
      if (ZelWmuCounter==1) ZelWmu=true;
      
      if ((ZelWel) && (ZelWmu)){
	continue;
      }
      if ((!ZelWel) && (!ZelWmu)) continue;
      
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
    

    if ((!ev3e) && (!ev3mu) && (!ev1e2mu) && (!ev2e1mu)) continue;

    //deltaR condition
    if (!passDeltaRWleptZlept(WZcandidates, phi, eta)) continue;
    numW++;  
    
    int type;
    //0-EEE, 1-EEM, 2-EMM, 3-MMM
    
    if (ev3e){
      num3e++;
      type=0;
    }
    if (ev2e1mu){
      num2e1mu++;
      type=1;
    }
    if (ev1e2mu){
      num1e2mu++;
      type=2;
    }

    if (ev3mu){
      num3mu++;
      type=3;
    }
    //////////////////////////////////////MET CUT//////////////////////////
    
    if ((cWZ->pfmet)<30) continue;
    //    if ((cWZ->pfmetTypeI)<30) continue;   ///CHANGE THIS

    ////////////////////MATRIX METHOD PART//////////////////////////

    int label= determineLabel(pass2012ICHEP, WZcandidates);
    float Nttt(0);
    if (label==1) Nttt=1; 
    
    if (ev3e){
      double w_event_3e=weight(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label);
      N_good_3e+=w_event_3e;
      //      error_3e+=(w_event_3e*w_event_3e);
      error_3e+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, w_event_3e);
      N_fake_3e+=(Nttt-w_event_3e);
      //      ferror_3e+=(Nttt-w_event_3e)*(Nttt-w_event_3e);
      ferror_3e+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, (Nttt-w_event_3e));
    }

 
    if (ev2e1mu){
      double w_event_2e1mu=weight(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label);
      N_good_2e1mu+=w_event_2e1mu;
      //      double tpe= thirdPartError(w_event_2e1mu, epsilon1, epsilon2, epsilon3, epsilonError1, epsilonError2, epsilonError3, p1, p2, p3);
      // sigma2e1mu=fpe+spe+tpe;
      //      error_2e1mu+=(w_event_2e1mu*w_event_2e1mu);
      error_2e1mu+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, w_event_2e1mu);
      N_fake_2e1mu+=(Nttt-w_event_2e1mu);
      //      ferror_2e1mu+=(Nttt-w_event_2e1mu)*(Nttt-w_event_2e1mu);
      ferror_2e1mu+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, (Nttt-w_event_2e1mu));
    }
  
    if (ev1e2mu){
      double w_event_1e2mu=weight(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label);
      N_good_1e2mu+=w_event_1e2mu;
      //      double tpe= thirdPartError(w_event_1e2mu, epsilon1, epsilon2, epsilon3, epsilonError1, epsilonError2, epsilonError3, p1, p2, p3);
      //sigma1e2mu=fpe+spe+tpe;
      //      error_1e2mu+=(w_event_1e2mu*w_event_1e2mu);
      error_1e2mu+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, w_event_1e2mu);
      N_fake_1e2mu+=(Nttt-w_event_1e2mu);
      //      ferror_1e2mu+=(Nttt-w_event_1e2mu)*(Nttt-w_event_1e2mu);
      ferror_1e2mu+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, (Nttt-w_event_1e2mu));
    }
    
    if (ev3mu){
      double w_event_3mu=weight(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label);
      N_good_3mu+=w_event_3mu;
      //      double tpe= thirdPartError(w_event_3mu, epsilon1, epsilon2, epsilon3, epsilonError1, epsilonError2, epsilonError3, p1, p2, p3);
      //sigma3mu=fpe+spe+tpe;
      //      error_3mu+=(w_event_3mu*w_event_3mu);
      error_3mu+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, w_event_3mu);
      N_fake_3mu+=(Nttt-w_event_3mu);
      //      ferror_3mu+=(Nttt-w_event_3mu)*(Nttt-w_event_3mu);
      ferror_3mu+=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label, (Nttt-w_event_3mu));
    }
    //    std::cout<<MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, type, pt, eta, label)<<std::endl;
     
    if (ev3e){
      numMET3e++;
    }
    if (ev2e1mu){
      numMET2e1mu++;
    }
    if (ev1e2mu){
      numMET1e2mu++;
    }

    if (ev3mu){
      numMET3mu++;
    }
        
    
   
  }
  std::cout<<"Zyield: "<<numZ<<std::endl;
  std::cout<<"Wyield: "<<numW<<std::endl;
  std::cout<<"W3e:     "<<num3e<<std::endl;
  std::cout<<"W2e1mu:  "<<num2e1mu<<std::endl;
  std::cout<<"W1e2mu:  "<<num1e2mu<<std::endl;
  std::cout<<"W3mu:    "<<num3mu<<std::endl;
  std::cout<<"****************************"<<std::endl;
  std::cout<<"3e:     "<<numMET3e<<std::endl;
  std::cout<<"2e1mu:  "<<numMET2e1mu<<std::endl;
  std::cout<<"1e2mu:  "<<numMET1e2mu<<std::endl;
  std::cout<<"3mu:    "<<numMET3mu<<std::endl;

  std::cout<<justCount<<std::endl;
  std::cout<<"/////////////////MATRIX METHOD RESULT//////////////////"<<std::endl;
  std::cout<<"N_good_3e:     "<<N_good_3e<<std::endl;
  std::cout<<"error_3e:     "<<sqrt(error_3e)<<std::endl;
  std::cout<<"N_good_2e1mu:  "<<N_good_2e1mu<<std::endl;
  std::cout<<"error_2e1mu:  "<<sqrt(error_2e1mu)<<std::endl;
  std::cout<<"N_good_1e2mu:  "<<N_good_1e2mu<<std::endl;
  std::cout<<"error_1e2mu:  "<<sqrt(error_1e2mu)<<std::endl;
  std::cout<<"N_good_3mu:    "<<N_good_3mu<<std::endl;
  std::cout<<"error_3mu:    "<<sqrt(error_3mu)<<std::endl;
  std::cout<<"*****************************************************"<<std::endl;
  std::cout<<"N_fake_3e:     "<<N_fake_3e<<std::endl;
  std::cout<<"ferror_3e:     "<<sqrt(ferror_3e)<<std::endl;
  std::cout<<"N_fake_2e1mu:  "<<N_fake_2e1mu<<std::endl;
  std::cout<<"ferror_2e1mu:  "<<sqrt(ferror_2e1mu)<<std::endl;
  std::cout<<"N_fake_1e2mu:  "<<N_fake_1e2mu<<std::endl;
  std::cout<<"ferror_1e2mu:  "<<sqrt(ferror_1e2mu)<<std::endl;
  std::cout<<"N_fake_3mu:    "<<N_fake_3mu<<std::endl;
  std::cout<<"ferror_3mu:    "<<sqrt(ferror_3mu)<<std::endl;  
  fout->cd();
  fout->Write();
  fout->Close();
}
