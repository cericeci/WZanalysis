#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

// Replace this with the new tree
#include "WZGenEvent.h"
//#include "WZ.h"
//#include "WZ2012Data.h"

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
bool Z_muons(WZGenEvent *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, float* pt, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, float* phi, float *eta);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax=-999.);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

float ScaleFactors(TH2F* MuonSF, TH2F* ElecSF,int type, int * WZCandidates, float *pt, float * eta);

float trigger3sameLeptons(int* eL, int* eT);

float trigger2sameLeptons(int* eL, int* eT);

float TriggerWeight(int* WZcandidates, TH2F* DoubleElLead, TH2F* DoubleMuLead, TH2F* DoubleElTrail, TH2F* DoubleMuTrail, int type,float* pt, float* eta);

void readFileFromList(TString fileList, std::vector<TString>* inputFile);

float AxeError(float weight, TH2F * MuonSF, TH2F * ElecSF, float pileUpWeight, int * WZCandidates,int type, float* pt, float* eta, float TriggerEff, float triggerError=0);

//int genInformation(WZ *cWZ, int* Wdecay, int* Zdecay);

int main()
{
  using namespace std;

  //write output numbers
  bool writeOutputNumbers(false);
  //  bool writeOutputNumbers(true);

  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;
  ofstream fileNumGEN; 
  if (writeOutputNumbers){
    fileNumGEN.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/numGEN.h");
    fileNumGEN<<"#ifndef numGEN_h"<<std::endl;
    fileNumGEN<<"#define numGEN_h"<<std::endl;
  }
  myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3eMC_Lucija.txt");
  myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts2e1muMC_Lucija.txt");
  myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts1e2muMC_Lucija.txt");
  myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3muMC_Lucija.txt");
  myfileAll.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_all.txt");


  /*
  myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3e_Lucija.txt");
  myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts2e1mu_Lucija.txt");
  myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts1e2mu_Lucija.txt");
  myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3mu_Lucija.txt");
  myfileAll.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/all_Lucija.txt");
  */

  //TFile * fout= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data.root", "RECREATE");

  TH1F * hZmassMu1         = new TH1F ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
  TH1F * hZmassEl1         = new TH1F ("hZmassEl1", "hZmassEl1", 100, 60, 120);  


  TH2F* MuonSF;
  TH2F* ElecSF;
  TH2F* DoubleElLead;
  TH2F* DoubleElTrail;
  TH2F* DoubleMuLead;
  TH2F* DoubleMuTrail;


  MuonSF=LoadHistogram("auxiliaryFiles/MuSF_2012.root", "h2inverted", "MuonSF");
  ElecSF=LoadHistogram("auxiliaryFiles/EleSF_2012.root", "h2inverted", "ElecSF");
  DoubleElLead=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleElLead", "DoubleElLead");
  DoubleMuLead=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleMuLead", "DoubleMuLead");
  DoubleElTrail=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleElTrail", "DoubleElTrail");
  DoubleMuTrail=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleMuTrail", "DoubleMuTrail");

  const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);
  const float luminosity(19602);
  float numZ(0), numW(0), numMET(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0), numMET2e1mu(0), numMET1e2mu(0), numMET3mu(0),numMET3e_1(0), numMET2e1mu_1(0), numMET1e2mu_1(0), numMET3mu_1(0);;
  float errorAxe3e(0), errorAxe2e1mu(0), errorAxe1e2mu(0), errorAxe3mu(0);
  float numMET3eGEN(0), numMET2e1muGEN(0), numMET1e2muGEN(0), numMET3muGEN(0), numMETSomethingGEN(0);
  float numMET3e_2(0), numMET2e1mu_2(0), numMET1e2mu_2(0), numMET3mu_2(0);
//  float numMET3eGEN2(0), numMET2e1muGEN2(0), numMET1e2muGEN2(0), numMET3muGEN2(0), numMETSomethingGEN2(0);
  float numTau3e(0), numTau2e1mu(0), numTau1e2mu(0), numTau3mu(0);
  //0=EEE, 1=EEM, 2=EMM, 3= MMM, 4=all
  float gen[5], genq[5], fact[5], factq[5], mixed[5], error[5], errorWholeZagreb[5], errorRest[5], errorSpanish[5];
  //  float sumPuWall[5], sumPuWsel[5], sumPuWfail[5], sumPuWallq[5], sumPuWselq[5], sumPuWfailq[5];
  for (int t=0; t<5; t++){
    gen[t]=0;
    genq[t]=0;    
    fact[t]=0;
    factq[t]=0;
    mixed[t]=0;
    errorWholeZagreb[t]=0;
    errorRest[t]=0;
    errorSpanish[t]=0;
  }


  int Wdecay[6]={0,0,0,0,0,0};
  int Zdecay[9]={0,0,0,0,0,0,0,0,0};
  
  float tauFactor;
  //type: 0=EEE, 1=EEM, 2=EMM, 3=MMM
  int type(-999);

  int testEl(0), testMu(0);
  double numberPileUp(0.0);

  std::vector<TString>inputName;
  TChain wz("latino");

  readFileFromList("wzNoSkim.files", &inputName);
  //  readFileFromList("wzNoFilter10.files", &inputName);
  
  for (int input=0; input< inputName.size(); input++){
    wz.Add(inputName[input]);
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZGenEvent *cWZ= new WZGenEvent(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;
  

  for  (Int_t k = 0; k<events /*&& k<1000*/;k++) {
    wz_tTree->GetEntry(k);


    //*******various factors
    float pileUpWeight=cWZ->puW;

    numberPileUp+=pileUpWeight;
    //rejecting run 201191   not needed for MC
    //    if (cWZ->run==201191) continue;

    //if (!(cWZ->trigger)) continue;
    
    ///////////////GEN YIELDS/////////////////
 
    //    double gen_channels[3][3]={{0,0,0}, {0,0,0}};
    bool isTau(false), is3e(false), is2e1mu(false), is1e2mu(false), is3mu(false);

    //WZchan== 0/eee, 1/eemu, 2/emumu, 3/mumumu, 4/tau
    if (((cWZ->MZ)>71.1876) && ((cWZ->MZ<111.1876)))
    {
      if ((cWZ->WZchan)== 0){
	numMET3eGEN+=pileUpWeight;
	gen[0]+=pileUpWeight;
	genq[0]+=(pileUpWeight*pileUpWeight);
	is3e=true;
      }
      if ((cWZ->WZchan)== 1){
	numMET2e1muGEN+=pileUpWeight;
	gen[1]+=pileUpWeight;
	genq[1]+=(pileUpWeight*pileUpWeight);
	is2e1mu=true;
      }
      if ((cWZ->WZchan)== 2){
	numMET1e2muGEN+=pileUpWeight;
      	gen[2]+=pileUpWeight;
	genq[2]+=(pileUpWeight*pileUpWeight);
	is1e2mu=true;
      }
      if ((cWZ->WZchan)== 3){
	  numMET3muGEN+=pileUpWeight;
	  gen[3]+=pileUpWeight;
	  genq[3]+=(pileUpWeight*pileUpWeight);
	  is3mu=true;
      }
      if ((cWZ->WZchan)==4){
	numMETSomethingGEN+=pileUpWeight;
	gen[4]+=pileUpWeight;
	genq[4]+=(pileUpWeight*pileUpWeight);
	isTau=true;
      }

    }
    //    int dummy=genInformation(cWZ, Wdecay, Zdecay);
    float pt[leptonNumber]={cWZ->pt1, cWZ->pt2, cWZ->pt3, cWZ->pt4};
    float bdt[leptonNumber]={cWZ->bdt1, cWZ->bdt2, cWZ->bdt3, cWZ->bdt4};
    float ch[leptonNumber]={cWZ->ch1, cWZ->ch2, cWZ->ch3, cWZ->ch4};
    float eta[leptonNumber]={cWZ->eta1, cWZ->eta2, cWZ->eta3, cWZ->eta4};
    float phi[leptonNumber]={cWZ->phi1, cWZ->phi2, cWZ->phi3, cWZ->phi4};
    int pass2012ICHEP[leptonNumber]={cWZ->pass2012ICHEP1, cWZ->pass2012ICHEP2, cWZ->pass2012ICHEP3, cWZ->pass2012ICHEP4};
    float iso[leptonNumber]={cWZ->iso1, cWZ->iso2, cWZ->iso3, cWZ->iso4};
    float isomva[leptonNumber]={cWZ->isomva1, cWZ->isomva2, cWZ->isomva3, cWZ->isomva4};
    float pdgid[leptonNumber]={cWZ->pdgid1, cWZ->pdgid2, cWZ->pdgid3, cWZ->pdgid4};
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
      if ((pt[i]>10) && (pt[i]!=-9999))
	lepNum++;
      if (((fabs(pdgid[i]))==11) && (pt[i]>10) && (!pass2012ICHEP[i]))
	testEl++;
      if ((fabs(pdgid[i])==13) && (pt[i]>10) && (!pass2012ICHEP[i]))
	testMu++;
      if (((fabs(pdgid[i]))==11) && (pt[i]>10) && (pass2012ICHEP[i])){
	good_electrons.push_back(i);
	v_nizEl[i].SetPtEtaPhiM(pt[i],eta[i], phi[i], electronMass);
	v_3Lepton=v_3Lepton+v_nizEl[i];
      }

      if ((fabs(pdgid[i])==13) && (pt[i]>10) && (pass2012ICHEP[i])){
      good_muons.push_back(i);
	v_nizMu[i].SetPtEtaPhiM(pt[i],eta[i], phi[i], muonMass);
	v_3Lepton=v_3Lepton+v_nizMu[i];
      }
    }
    if (lepNum!=3) continue;

    if (v_3Lepton.M()<100) continue;

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

    numZ+=pileUpWeight;

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
    numW+=pileUpWeight;  
    
    if (ev3e){
      num3e+=pileUpWeight;
      type=0; 
    }
    if (ev2e1mu){
      num2e1mu+=pileUpWeight;
      type=1;
    }
    if (ev1e2mu){
      num1e2mu+=pileUpWeight;
      type=2;
    }

    if (ev3mu){
      num3mu+=pileUpWeight;
      type=3;
    }
    //////////////////////////////////////MET CUT//////////////////////////
    
    if ((cWZ->pfmet)<30) continue;
    //    if ((cWZ->pfmetTypeI)<30) continue;   ///CHANGE THIS

    float weight;
    if (ev3e){
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 0, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 0, pt, eta);

      numMET3e+=weight;
      numMET3e_1++;
      if (isTau) numTau3e+=weight;
      if (is3e) {
	numMET3e_2+=weight;
	fact[0]+=weight;
	factq[0]+=(weight*weight);
	mixed[0]+=(weight*pileUpWeight);
      }
      float triggerEff=TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 0, pt, eta);
      errorAxe3e+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 0, pt, eta, triggerEff);
      errorRest[0]+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 0, pt, eta, triggerEff);
      myfile3e<<cWZ->run<<":"<<cWZ->event<<":"<<pileUpWeight<<":"<<ScaleFactors(MuonSF, ElecSF, 0, WZcandidates, pt, eta)<<":"<<triggerEff<<":"<<std::endl;      
    }
    
    if (ev2e1mu){
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 1, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 1, pt, eta);
      
      numMET2e1mu+=weight;
      numMET2e1mu_1++;;

      if (isTau) numTau2e1mu+=weight;

      float triggerEff=TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 1, pt, eta);
	//      myfile2e1mu<<cWZ->run<<":"<<cWZ->event<<":"<<pileUpWeight<<":"<<ScaleFactors(MuonSF, ElecSF, 1, WZcandidates, pt, eta, f1, f2, f3)<<":"<<f1<<":"<<f2<<":"<<f3<<":"<<std::endl;

      if (is2e1mu){
	numMET2e1mu_2+=weight;
	fact[1]+=weight;
	factq[1]+=(weight*weight);
	mixed[1]+=(weight*pileUpWeight);   
      }
      errorAxe2e1mu+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 1, pt, eta, triggerEff);
      errorRest[1]+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 1, pt, eta, triggerEff);
      myfile2e1mu<<cWZ->run<<":"<<cWZ->event<<":"<<pileUpWeight<<":"<<ScaleFactors(MuonSF, ElecSF, 1, WZcandidates, pt, eta)<<":"<<triggerEff<<":"<<std::endl;
    }
    
    if (ev1e2mu){
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 2, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 2, pt, eta);

      numMET1e2mu+=weight;
      numMET1e2mu_1++;
      if (isTau) numTau1e2mu+=weight;
      if (is1e2mu) {
	numMET1e2mu_2+=weight;
	fact[2]+=weight;
	factq[2]+=(weight*weight);
	mixed[2]+=(weight*pileUpWeight);
      }
      float triggerEff=TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 2, pt, eta);
      errorAxe1e2mu+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 2, pt, eta, triggerEff);
      errorRest[2]+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 2, pt, eta, triggerEff);
      myfile1e2mu<<cWZ->run<<":"<<cWZ->event<<":"<<pileUpWeight<<":"<<ScaleFactors(MuonSF, ElecSF, 2, WZcandidates, pt, eta)<<":"<<triggerEff<<":"<<std::endl;
    }

    if (ev3mu){
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 3, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 3, pt, eta);
      numMET3mu+=weight;
      numMET3mu_1++;
      if (isTau) numTau3mu+=weight;
      if (is3mu) {
	numMET3mu_2+=weight;
	fact[3]+=weight;
	factq[3]+=(weight*weight);
	mixed[3]+=(weight*pileUpWeight);
      }
      float triggerEff=TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 3, pt, eta);
      errorAxe3mu+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 3, pt, eta, triggerEff);
      errorRest[3]+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, 3, pt, eta, triggerEff);
      myfile3mu<<cWZ->run<<":"<<cWZ->event<<":"<<pileUpWeight<<":"<<ScaleFactors(MuonSF, ElecSF, 3, WZcandidates, pt, eta)<<":"<<triggerEff<<":"<<std::endl;
    }
 
    
  }
  std::cout<<"Number of events pile up weighted: "<<numberPileUp<<std::endl;

  float denominator=numMET3eGEN+numMET2e1muGEN+numMET1e2muGEN+numMET3muGEN+numMETSomethingGEN;
  gen[5]=gen[0]+gen[1]+gen[2]+gen[3]+gen[4];

  std::cout<<"**********************GEN**INFORMATION**********************"<<std::endl;
  std::cout<<"****W decays*******"<<std::endl;
  for (int w=0;w<6; w++){
    std::cout<<w<<"  : "<<Wdecay[w]<<std::endl;
  }




  for (int e=0; e<4; e++){
    error[e]=sqrt((gen[e]*gen[e]*factq[e]-2*gen[e]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[e])/pow(gen[e],4));
    errorWholeZagreb[e]=sqrt(errorRest[e]/(gen[e]*gen[e])+(gen[e]*gen[e]*factq[e]-2*gen[e]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[e])/pow(gen[e],4));
    errorSpanish[e]=sqrt((gen[4]*gen[4]*factq[e]-2*gen[4]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[4])/pow(gen[4],4));
  }
  std::cout<<"Acceptance efficiancy: "<<std::endl;
  std::cout<<"3e     "<< numMET3e/denominator    <<" nWZ: "<<numMET3e     <<" nGEN: "<<denominator<<std::endl;
  std::cout<<"2e1mu  "<< numMET2e1mu/denominator <<" nWZ: "<<numMET2e1mu  <<" nGEN: "<<denominator<<std::endl;
  std::cout<<"1e2mu  "<< numMET1e2mu/denominator <<" nWZ: "<<numMET1e2mu  <<" nGEN: "<<denominator<<std::endl;
  std::cout<<"3mu    "<< numMET3mu/denominator   <<" nWZ: "<<numMET3mu    <<" nGEN: "<<denominator<<std::endl;


  std::cout<<"test!!: "<<gen[5]<<"="<<denominator<<std::endl;
  std::cout<<"errorZagreb"<<std::endl;
  std::cout<<"3e:   "<<error[0]<<std::endl;
  std::cout<<"2e1mu "<<error[1]<<std::endl;
  std::cout<<"1e2mu "<<error[2]<<std::endl;
  std::cout<<"3mu   "<<error[3]<<std::endl;
  std::cout<<"complete..."<<std::endl;
  
  std::cout<<"3e:   "<<errorWholeZagreb[0]<<std::endl;
  std::cout<<"2e1mu "<<errorWholeZagreb[1]<<std::endl;
  std::cout<<"1e2mu "<<errorWholeZagreb[2]<<std::endl;
  std::cout<<"3mu   "<<errorWholeZagreb[3]<<std::endl;

  std::cout<<"errorSpanish"<<std::endl;
  std::cout<<"3e:   "<<errorSpanish[0]<<std::endl;
  std::cout<<"2e1mu "<<errorSpanish[1]<<std::endl;
  std::cout<<"1e2mu "<<errorSpanish[2]<<std::endl;
  std::cout<<"3mu   "<<errorSpanish[3]<<std::endl;
  std::cout<<"complete..."<<std::endl;

  float normalization=(1.058*luminosity)/numberPileUp;
  std::cout<<"YIELDS:"<<std::endl;
  std::cout<<numMET3e*normalization<<" , "<<numMET3e<<" , "<<numMET3e_1<<std::endl;
  std::cout<<numMET2e1mu*normalization<<" , "<<numMET2e1mu<<" , "<<numMET2e1mu_1<<std::endl;
  std::cout<<numMET1e2mu*normalization<<" , "<<numMET1e2mu<<" , "<<numMET1e2mu_1<<std::endl;
  std::cout<<numMET3mu*normalization<<" , "<<numMET3mu<<" , "<<numMET3mu_1<<std::endl;

  //writing in file
  if (writeOutputNumbers){
    fileNumGEN<<"#define dAxe3eS "<< numMET3e/denominator<<std::endl;
    fileNumGEN<<"#define dAxe2e1muS "<<numMET2e1mu/denominator<<std::endl;
    fileNumGEN<<"#define dAxe1e2muS "<<numMET1e2mu/denominator<<std::endl;
    fileNumGEN<<"#define dAxe3muS "<<numMET3mu/denominator<<std::endl;
    
    fileNumGEN<<"#define dAxe3eZ "<< numMET3e/numMET3eGEN<<std::endl;
    fileNumGEN<<"#define dAxe2e1muZ "<<numMET2e1mu/numMET2e1muGEN<<std::endl;
    fileNumGEN<<"#define dAxe1e2muZ "<<numMET1e2mu/numMET1e2muGEN<<std::endl;
    fileNumGEN<<"#define dAxe3muZ "<<numMET3mu/numMET3muGEN<<std::endl;
    
    fileNumGEN<<"#define dAxe3eZ_2 "<< numMET3e_2/numMET3eGEN<<std::endl;
    fileNumGEN<<"#define dAxe2e1muZ_2 "<<numMET2e1mu_2/numMET2e1muGEN<<std::endl;
    fileNumGEN<<"#define dAxe1e2muZ_2 "<<numMET1e2mu_2/numMET1e2muGEN<<std::endl;
    fileNumGEN<<"#define dAxe3muZ_2 "<<numMET3mu_2/numMET3muGEN<<std::endl;
    
    fileNumGEN<<"#define dtauFactor3e "<<numTau3e/numMET3e<<std::endl;
    fileNumGEN<<"#define dtauFactor2e1mu "<<numTau2e1mu/numMET2e1mu<<std::endl;
    fileNumGEN<<"#define dtauFactor1e2mu "<<numTau1e2mu/numMET1e2mu<<std::endl;
    fileNumGEN<<"#define dtauFactor3mu "<<numTau3mu/numMET3mu<<std::endl;
  
    fileNumGEN<<"#define dsAxe3eS "<< errorSpanish[0]<<std::endl;
    fileNumGEN<<"#define dsAxe2e1muS "<<errorSpanish[1]<<std::endl;
    fileNumGEN<<"#define dsAxe1e2muS "<<errorSpanish[2]<<std::endl;
    fileNumGEN<<"#define dsAxe3muS "<<errorSpanish[3]<<std::endl;
    
    fileNumGEN<<"#define dsAxe3eZ "<< error[0]<<std::endl;
    fileNumGEN<<"#define dsAxe2e1muZ "<<error[1]<<std::endl;
    fileNumGEN<<"#define dsAxe1e2muZ "<<error[2]<<std::endl;
    fileNumGEN<<"#define dsAxe3muZ "<<error[3]<<std::endl;
    
    fileNumGEN<<"#define dstauFactor3e "<<0<<std::endl;
    fileNumGEN<<"#define dstauFactor2e1mu "<<0<<std::endl;
    fileNumGEN<<"#define dstauFactor1e2mu "<<0<<std::endl;
    fileNumGEN<<"#define dstauFactor3mu "<<0<<std::endl;
  
  

  fileNumGEN.close();
  }

}
