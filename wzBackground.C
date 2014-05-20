#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "UnfoldingHistogramFactory.h"
// Replace this with the new tree
//#include "WZEvent.h"
//#include "WZ.h"
//#include "WZ2012Data.h"
#include "WZEventMCOld.h"

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
void readChainFromList(TString fileList, TChain * chain);

bool Z_muons(WZ *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, float* pt, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, float* phi, float *eta);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

float ScaleFactors(TH2F* MuonSF, TH2F* ElecSF,int type, int * WZCandidates, float *pt, float * eta);

float trigger3sameLeptons(int* eL, int* eT);

float trigger2sameLeptons(int* eL, int* eT);

float TriggerWeight(int* WZcandidates, TH2F* DoubleElLead, TH2F* DoubleMuLead, TH2F* DoubleElTrail, TH2F* DoubleMuTrail, int type,float* pt, float* eta);

void readFileFromList(TString fileList,std::vector<TString> * inputFile);

int main()
{
  using namespace std;
  
  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;
  ofstream fileNumMC;
  //***write output numbers
  //  bool writeOutputNumbers(false);
  bool writeOutputNumbers(false);

  if (writeOutputNumbers){
    fileNumMC.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numMC.h");
    fileNumMC<<"#ifndef numMC_h"<<std::endl;
    fileNumMC<<"#define numMC_h"<<std::endl;
  }
    /*
  myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3e_Lucija.txt");
  myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts2e1mu_Lucija.txt");
  myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts1e2mu_Lucija.txt");
  myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3mu_Lucija.txt");
  myfileAll.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/all_Lucija.txt");
  */

  //TFile * fout= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data.root", "RECREATE");



  const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);
  const float luminosity(19.602);

  //check!!!!!


  //this had to be read only once
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

  std::vector<TString>name;
  std::vector<TChain*> chain;
  std::vector<TString> files;

  
  //  files.push_back("WZ.files");
  files.push_back("ZZ.files");
  files.push_back("Zgamma.files");
  files.push_back("WV.files");
  files.push_back("VVV.files");
  //  files.push_back("top.files");
  //files.push_back("Zjets.files");

  //  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WZ.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/ZZ.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zgamma.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WV.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/VVV.root");
  //  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/top.root");
  //name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zjets.root");





  for (int file=0; file<files.size() /*&& file<1*/; file++){

    float numZ(0), numW(0), numMET(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0.0), numMET2e1mu(0.0), numMET1e2mu(0.0), numMET3mu(0.0);
    float numMET3ectrl(0), numMET2e1muctrl(0), numMET1e2muctrl(0), numMET3muctrl(0);
    float numMET3eGEN(0), numMET2e1muGEN(0), numMET1e2muGEN(0), numMET3muGEN(0);
    float error3e(0), error2e1mu(0), error1e2mu(0), error3mu(0);
    TFile * fout = new TFile(name[file],"RECREATE");

    TH1F * hZmassMu1         = new TH1F ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
    TH1F * hZmassEl1         = new TH1F ("hZmassEl1", "hZmassEl1", 100, 60, 120);  
    TH1F * hZmassMuWel3      = new TH1F ("hZmassMuWel3", "hZmassMuWel3", 40, 60, 120);
    TH1F * hZmassMuWmu3      = new TH1F ("hZmassMuWmu3", "hZmassMuWmu3", 40, 60, 120);
    TH1F * hZmassElWmu3      = new TH1F ("hZmassElWmu3", "hZmassElWmu3", 40, 60, 120);
    TH1F * hZmassElWel3      = new TH1F ("hZmassElWel3", "hZmassElWel3", 40, 60, 120);
    TH1F * hMETMu1                = new TH1F("hMETMu1", "hMETMu1",150, 0, 300);
    TH1F * hMETEl1                = new TH1F("hMETEl1", "hMETEl1",150, 0, 300);
    TH1F * hMETMu1_s         = new TH1F("hMETMu1_s", "hMETMu1_s",150, 0, 150);
    TH1F * hMETEl1_s         = new TH1F("hMETEl1_s", "hMETEl1_s",150, 0, 150);


    const int nChannels(4);
    TH1D * hZpt[nChannels];
    TH1D * hLeadingJetPt[nChannels]; 

    //TH1D * hZpt= UnfoldingHistogramFactory::createZPtHistogram("name", "name");
    //    LeadingJetname<<"LeadingJetPt_";
    //Zptname<<"Zpt_";
    
    for (int hist=0; hist<nChannels; hist++){
      std::ostringstream Zptname, LeadingJetname;
      LeadingJetname<<"LeadingJetPt_"<<(hist+1);
      Zptname<<"Zpt_"<<(hist+1);
      hZpt[hist]     = UnfoldingHistogramFactory::createZPtHistogram(Zptname.str().c_str(), Zptname.str().c_str());
      hLeadingJetPt[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetname.str().c_str(),LeadingJetname.str().c_str());
  }
    
//type: 0=EEE, 1=EEM, 2=EMM, 3=MMM
    
    int type(-999);

    std::vector<TString> inputName;
    readFileFromList(files[file], &inputName);  
    TChain wz("latino");
    
    for (int input=0; input< inputName.size(); input++){
      wz.Add(inputName[input]);
    }



    TTree *wz_tTree=(TTree*)&wz;
    //WZ *cWZ= new WZ(wz_tTree);
    WZEventMCOld *cWZ= new WZEventMCOld(wz_tTree);
    Int_t events= wz_tTree->GetEntries();
    
    //    WZEventMCOld * wzevt;

    std::cout<<"number of events: "<<events << std::endl;
    
    std::cout<<name[file]<<std::endl;

    


    for  (Int_t k = 0; k<events /*&& k<10000*/;k++) {
      wz_tTree->GetEntry(k);    
      cWZ->ReadEvent();
      
      float xs_weight(0);    

    //*******various factors
    float pileUpWeight=cWZ->puW;
    
    //rejecting run 201191
    //    if (cWZ->run==201191) continue;

    //if (!(cWZ->trigger)) continue;

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
    
    int WZcandidates[3]; 

float weight;
    //xs_weight=(cWZ->baseW)*luminosity;
    xs_weight=(1.0 + 0.6 * ((cWZ->dataset) >= 82 && (cWZ->dataset) <= 84)) * (cWZ->baseW) * luminosity;
    if ((cWZ->dataset)==89) xs_weight *= (0.01968  / 0.0192);
    if ((cWZ->dataset)==90) xs_weight *= (0.005527 / 0.00459);
    if ((cWZ->dataset)==91) xs_weight *= (0.05795  / 0.0633);
    if ((cWZ->dataset)==92) xs_weight *= (0.08058  / 0.0822);
    if ((cWZ->dataset)==93) xs_weight *= (0.232    / 0.232);
    if ((cWZ->dataset)==94) xs_weight *= (0.2057 / 0.174);

    
    //here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
    int lepNum(0);
    for (int i=0; i<leptonNumber; i++){
      if ((pt[i]>10) && (pt[i]!=-9999))
	lepNum++;
      if ((bdt[i]<100) && (pass2012ICHEP[i]) && (pt[i]>10)){
	good_electrons.push_back(i);
	v_nizEl[i].SetPtEtaPhiM(pt[i],eta[i], phi[i], electronMass);
	v_3Lepton=v_3Lepton+v_nizEl[i];
      }
      if ((bdt[i]>100) && (pass2012ICHEP[i]) && (pt[i]>10)){
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
    //    if (foundZmu){
    //hZmassMu1->Fill(massMu, pileUpWeight*xs_weight);
    //}

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
 
    if (foundZmu){
      hZmassMu1->Fill(massMu, pileUpWeight*xs_weight);
      hMETMu1->Fill(cWZ->pfmet, pileUpWeight*xs_weight);
      hMETMu1_s->Fill(cWZ->pfmet, pileUpWeight*xs_weight);   
    }    

    if (foundZel){
      hZmassEl1->Fill(massEl, pileUpWeight*xs_weight);
      hMETEl1->Fill(cWZ->pfmet, pileUpWeight*xs_weight);
      hMETEl1_s->Fill(cWZ->pfmet, pileUpWeight*xs_weight);
    }
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

    //improving this
    bool ev[4]={false, false, false, false};
    
    if (ev3e){
      ev[0]=true;
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 0, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 0, pt, eta)*xs_weight;
      numMET3e+=weight;
      error3e+=(weight*weight);
      numMET3ectrl+=pileUpWeight;
      hZmassElWel3->Fill(massEl, pileUpWeight*xs_weight);

    }
    
    if (ev2e1mu){
      ev[1]=true;
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 1, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 1, pt, eta)*xs_weight;
      numMET2e1mu+=weight;
      error2e1mu+=(weight*weight);
      numMET2e1muctrl+=pileUpWeight;
      hZmassElWmu3->Fill(massEl, pileUpWeight*xs_weight);
    }
    
    if (ev1e2mu){
      ev[2]=true;
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 2, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 2, pt, eta)*xs_weight;
      numMET1e2mu+=weight;
      error1e2mu+=(weight*weight);
      numMET1e2muctrl+=pileUpWeight;
      hZmassMuWel3->Fill(massMu,pileUpWeight*xs_weight);
    }

    if (ev3mu){
      ev[3]=true;
      weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, 3, WZcandidates, pt, eta)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, 3, pt, eta)*xs_weight;
      numMET3mu+=weight;
      error3mu+=(weight*weight);
      numMET3muctrl+=pileUpWeight;
      hZmassMuWmu3->Fill(massMu,pileUpWeight*xs_weight);
    }
        

    

    int nRecoJets = 0;
    float leadingRecoJetPt = -9999.;
    int leadingRecojet = -1;
    
    
    for (int i=0; i<cWZ->recoJets.size(); i++) {

      if (cWZ->recoJets[i].Pt() > 30 && fabs(cWZ->recoJets[i].Eta()) < 2.5) {

	bool closeToLepton = false;
	float drMin = 3.;
	for (int il=0; il<cWZ->leptons.size(); il++) {
	  if (cWZ->recoJets[i].DeltaR(cWZ->leptons[il])<0.5) {
	    closeToLepton = true;
	  }
	}
	if (closeToLepton) continue;
	
	nRecoJets++;
	if (cWZ->recoJets[i].Pt() > leadingRecoJetPt) {

	  leadingRecoJetPt = cWZ->recoJets[i].Pt();
	  leadingRecojet = i;
	}
      }
      
    }
    
    //std::cout<<leadingRecoJetPt<<std::endl;
    
    for (int filH=0; filH<nChannels; filH++){
      if (ev[filH]){
	hZpt[filH]->Fill(Zpt, weight);
	hLeadingJetPt[filH]->Fill(leadingRecoJetPt, weight);
      }
    }
    

    }
    
    
    

    std::cout<<"****************************(normalized to the luminosity, with scale factors and trigger eff)"<<std::endl;
    std::cout<<"3e:     "<<numMET3e<<"+/-"<<sqrt(error3e)<<std::endl;
    std::cout<<"2e1mu:  "<<numMET2e1mu<<"+/-"<<sqrt(error2e1mu)<<std::endl;
    std::cout<<"1e2mu:  "<<numMET1e2mu<<"+/-"<<sqrt(error1e2mu)<<std::endl;
    std::cout<<"3mu:    "<<numMET3mu<<"+/-"<<sqrt(error3mu)<<std::endl;


     fout->cd();
     fout->Write();

     /*
     for (int histDel=0; histDel<nChannels<4; histDel++){
       delete hZpt[histDel];
       delete hLeadingJetPt[histDel];
     }
     */
     delete hZmassMu1;
     delete hZmassEl1;
     delete hZmassMuWel3;
     delete hZmassMuWmu3;
     delete hZmassElWmu3;
     delete hZmassElWel3;
     delete hMETMu1;
     delete hMETEl1;
     delete hMETMu1_s;
     delete hMETEl1_s;
     fout->Close();    
    //writing in file:
    if (writeOutputNumbers){
      if (file==1){
	fileNumMC<<"#define dNZZ_3e "<<numMET3e<<std::endl;
	fileNumMC<<"#define dNZZ_2e1mu "<<numMET2e1mu<<std::endl;
	fileNumMC<<"#define dNZZ_1e2mu "<<numMET1e2mu<<std::endl;
	fileNumMC<<"#define dNZZ_3mu "<<numMET3mu<<std::endl;
	fileNumMC<<"#define dsNZZ_3e "<<error3e<<std::endl;
	fileNumMC<<"#define dsNZZ_2e1mu "<<error2e1mu<<std::endl;
	fileNumMC<<"#define dsNZZ_1e2mu "<<error1e2mu<<std::endl;
	fileNumMC<<"#define dsNZZ_3mu "<<error3mu<<std::endl;
      }
      
      if (file==2){
	fileNumMC<<"#define dNZgamma_3e "<<numMET3e<<std::endl;
	fileNumMC<<"#define dNZgamma_2e1mu "<<numMET2e1mu<<std::endl;
	fileNumMC<<"#define dNZgamma_1e2mu "<<numMET1e2mu<<std::endl;
	fileNumMC<<"#define dNZgamma_3mu "<<numMET3mu<<std::endl;
	fileNumMC<<"#define dsNZgamma_3e "<<error3e<<std::endl;
	fileNumMC<<"#define dsNZgamma_2e1mu "<<error2e1mu<<std::endl;
	fileNumMC<<"#define dsNZgamma_1e2mu "<<error1e2mu<<std::endl;
	fileNumMC<<"#define dsNZgamma_3mu "<<error3mu<<std::endl;
      }
      if (file==3){
	fileNumMC<<"#define dNWV_3e "<<numMET3e<<std::endl;
	fileNumMC<<"#define dNWV_2e1mu "<<numMET2e1mu<<std::endl;
	fileNumMC<<"#define dNWV_1e2mu "<<numMET1e2mu<<std::endl;
	fileNumMC<<"#define dNWV_3mu "<<numMET3mu<<std::endl;
	fileNumMC<<"#define dsNWV_3e "<<error3e<<std::endl;
	fileNumMC<<"#define dsNWV_2e1mu "<<error2e1mu<<std::endl;
	fileNumMC<<"#define dsNWV_1e2mu "<<error1e2mu<<std::endl;
	fileNumMC<<"#define dsNWV_3mu "<<error3mu<<std::endl;
      }
      if (file==4){
	fileNumMC<<"#define dNVVV_3e "<<numMET3e<<std::endl;
	fileNumMC<<"#define dNVVV_2e1mu "<<numMET2e1mu<<std::endl;
	fileNumMC<<"#define dNVVV_1e2mu "<<numMET1e2mu<<std::endl;
	fileNumMC<<"#define dNVVV_3mu "<<numMET3mu<<std::endl;
	fileNumMC<<"#define dsNVVV_3e "<<error3e<<std::endl;
	fileNumMC<<"#define dsNVVV_2e1mu "<<error2e1mu<<std::endl;
	fileNumMC<<"#define dsNVVV_1e2mu "<<error1e2mu<<std::endl;
	fileNumMC<<"#define dsNVVV_3mu "<<error3mu<<std::endl;
      }
    }
  }
  fileNumMC.close();
    //fout->Write();
  //fout->Close();
}
