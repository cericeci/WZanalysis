#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "HistogramFactory.h"
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
bool Z_muons(WZ2012Data *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, TLorentzVector * analysisLepton, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, TLorentzVector * analysisLepton);
int main()
{
  using namespace std;
  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;
  ofstream fileNumData;

  fileNumData.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numData_met.h");
  //  myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts3e_Lucija.txt");
  myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts3e_Lucija.txt");
  myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts2e1mu_Lucija.txt");
  //  myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts2e1mu_Lucija.txt");
  myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts1e2mu_Lucija.txt");
  //  myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts1e2mu_Lucija.txt");
  myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts3mu_Lucija.txt");
  //  myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/METcomparison/WZ_evts3mu_Lucija.txt");
  myfileAll.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/all_Lucija.txt");

  TFile * fout= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data.root", "RECREATE");

  //  bool writeOutputNumbers(false);
  bool writeOutputNumbers(true);
  if (writeOutputNumbers){
    fileNumData<<"#ifndef numData_h"<<std::endl;
    fileNumData<<"#define numData_h"<<std::endl;
  }

  TH1F * hZmassFinal[4];
  
  //  for (int histos=0; histos<4; histos++)
  //numbering convention: 2-after W selection, 3-after MET cut(final selection)

  const int nChannels1(4);

  TH1F * hZmass1[nChannels1];
  TH1F * hZmass2[nChannels1];
  TH1F * hMET1[nChannels1];
  TH1F * hMET2[nChannels1];
  TH1F * hZpt_my1[nChannels1];
  TH1F * hZpt_my2[nChannels1];
  TH1F * hLeadingJetPt_my1[nChannels1];
  TH1F * hLeadingJetPt_my2[nChannels1];
  TH1F * hNjets[nChannels1];
  TH1F * hDeltaPhi[nChannels1];
  
  for (int myhist=0; myhist<nChannels1; myhist++){
    std::ostringstream nZmass1, nMET1, nZpt1, nLeadingJetPt1, nNjets, nDeltaPhi;
    std::ostringstream nZmass2, nMET2, nZpt2, nLeadingJetPt2;
    nZmass1<<"hZmass1_"<<myhist;
    nZmass2<<"hZmass2_"<<myhist;
    nMET1<<"hMET1_"<<myhist;
    nMET2<<"hMET2_"<<myhist;
    nZpt1<<"hZpt1_"<<myhist;
    nZpt2<<"hZpt2_"<<myhist;
    nLeadingJetPt1<<"hLeadingJetPt1_"<<myhist;
    nLeadingJetPt2<<"hLeadingJetPt2_"<<myhist;
    nNjets<<"hNjets_"<<myhist;
    nDeltaPhi<<"hDeltaPhi"<<myhist;
    
    hZmass1[myhist]           = HistogramFactory::createZmassHisto(nZmass1.str().c_str(), nZmass1.str().c_str());
    hZmass2[myhist]           = HistogramFactory::createZmassHisto(nZmass2.str().c_str(), nZmass2.str().c_str());
    hMET1[myhist]             = HistogramFactory::createMETHisto(nMET1.str().c_str(), nMET1.str().c_str());
    hMET2[myhist]             = HistogramFactory::createMETHisto(nMET2.str().c_str(), nMET2.str().c_str());
    hZpt_my1[myhist]          = HistogramFactory::createZptHisto(nZpt1.str().c_str(), nZpt1.str().c_str());
    hZpt_my2[myhist]          = HistogramFactory::createZptHisto(nZpt2.str().c_str(), nZpt2.str().c_str());
    hLeadingJetPt_my1[myhist] = HistogramFactory::createLeadingJetptHisto(nLeadingJetPt1.str().c_str(), nLeadingJetPt1.str().c_str());
    hLeadingJetPt_my2[myhist] = HistogramFactory::createLeadingJetptHisto(nLeadingJetPt2.str().c_str(), nLeadingJetPt2.str().c_str());
    hNjets[myhist]           =  HistogramFactory::createNjetsHisto(nNjets.str().c_str(), nNjets.str().c_str());
    hDeltaPhi[myhist]        =  HistogramFactory::createDeltaPhi(nDeltaPhi.str().c_str(), nDeltaPhi.str().c_str());
    
  }

  /*
  TH1F * hZmassMu1         = new TH1F ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
  TH1F * hZmassEl1         = new TH1F ("hZmassEl1", "hZmassEl1", 100, 60, 120);  
  TH1F * hZmassMuWel2      = new TH1F ("hZmassMuWel2", "hZmassMuWel2", 40, 60, 120);
  TH1F * hZmassMuWmu2      = new TH1F ("hZmassMuWmu2", "hZmassMuWmu2", 40, 60, 120);
  TH1F * hZmassElWmu2      = new TH1F ("hZmassElWmu2", "hZmassElWmu2", 40, 60, 120);
  TH1F * hZmassElWel2      = new TH1F ("hZmassElWel2", "hZmassElWel2", 40, 60, 120);
  
  TH1F * hZmassMuWel3      = new TH1F ("hZmassMuWel3", "hZmassMuWel3", 40, 60, 120);
  TH1F * hZmassMuWmu3      = new TH1F ("hZmassMuWmu3", "hZmassMuWmu3", 40, 60, 120);
  TH1F * hZmassElWmu3      = new TH1F ("hZmassElWmu3", "hZmassElWmu3", 40, 60, 120);
  TH1F * hZmassElWel3      = new TH1F ("hZmassElWel3", "hZmassElWel3", 40, 60, 120);
  TH1F * hMETMu1_s         = new TH1F("hMETMu1_s", "hMETMu1_s",150, 0, 150);
  TH1F * hMETEl1_s         = new TH1F("hMETEl1_s", "hMETEl1_s",150, 0, 150);

  TH1F * hZptMu            = new TH1F("hZptMu", "hZptMu",40, 0, 400);
  TH1F * hZptEl            = new TH1F("hZptEl", "hZptEl",40, 0, 400);
  TH1F * hZpt3e            = new TH1F("hZpt3e", "hZpt3e",40, 0, 400);
  TH1F * hZpt2e1mu         = new TH1F("hZpt2e1mu", "hZpt2e1mu",40, 0, 400);
  TH1F * hZpt1e2mu         = new TH1F("hZpt1e2mu", "hZpt1e2mu",40, 0, 400);
  TH1F * hZpt3mu           = new TH1F("hZpt3mu", "hZpt3mu",40, 0, 400);
  */
    
    const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);
  int numZ(0), numW(0), numMET(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0), numMET2e1mu(0), numMET1e2mu(0), numMET3mu(0);
  

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
  


  
  TFile * forSenka1= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data_aTGC_8TeV_3e_1.root", "RECREATE");
  TFile * forSenka2= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data_aTGC_8TeV_2e_1.root", "RECREATE");
  TFile * forSenka3= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data_aTGC_8TeV_2mu_1.root", "RECREATE");  
  TFile * forSenka4= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data_aTGC_8TeV_3mu_1.root", "RECREATE");

  //  TTree* wz_aTGC = new TTree("T", "wz_aTGC");
  forSenka1->cd();
  TTree* wz_aTGC1 = new TTree("tgcTree", "wz_aTGC1");
  forSenka2->cd();
  TTree* wz_aTGC2 = new TTree("tgcTree", "tgcTree");
  forSenka3->cd();
  TTree* wz_aTGC3 = new TTree("tgcTree", "tgcTree");
  forSenka4->cd();
  TTree* wz_aTGC4 = new TTree("tgcTree", "tgcTree");
  // double PtZ(0);
  float PtZ_eee(0), PtZ_eem(0), PtZ_emm(0), PtZ_mmm(0);
  wz_aTGC1->Branch("PtZ", &PtZ_eee, "PtZ/F");
  wz_aTGC2->Branch("PtZ", &PtZ_eem, "PtZ/F");
  wz_aTGC3->Branch("PtZ", &PtZ_emm, "PtZ/F");
  wz_aTGC4->Branch("PtZ", &PtZ_mmm, "PtZ/F");
  
  for (int input=0; input< inputName.size(); input++){
    wz.Add(inputName[input]);
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZ2012Data *cWZ= new WZ2012Data(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  

  std::cout<<"number of events: "<<events << std::endl;


  for  (Int_t k = 0; k<events /*&& k<5000000*/;k++) {
    wz_tTree->GetEntry(k);

    double met= cWZ->pfmetTypeI;
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

    TLorentzVector analysisLepton[leptonNumber];    
    int WZcandidates[3]; //->here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
    for (int i1=0; i1<leptonNumber; i1++){
      if ((bdt[i1]<100) && (pt[i1]>10)){    
	analysisLepton[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], electronMass);
      }
      if ((bdt[i1]>100) && (pt[i1]>10)){
	analysisLepton[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], muonMass);
      }
    }

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

    //this maybe should be in the end
    if (lepNum!=3) continue;

    if (v_3Lepton.M()<100) continue;

    justCount++;
    bool foundZel(false), foundZmu(false);
    double massMu(-999), massEl(0), Zpt(0);



    foundZmu=  Z_muons(cWZ, &good_muons, WZcandidates, v_nizMu, analysisLepton, ch, massMu, Zpt);


    foundZel= Z_muons(cWZ, &good_electrons, WZcandidates, v_nizEl, analysisLepton, ch, massEl, Zpt);

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

    double Zmass;
    if (foundZel) Zmass=massEl;
    else Zmass=massMu;
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
    bool ev[4]={false, false, false, false};    
    if (ev3e) ev[0]=true;
    if (ev2e1mu) ev[1]=true;
    if (ev1e2mu) ev[2]=true;
    if (ev3mu) ev[3]=true;
    


    if ((!ev3e) && (!ev3mu) && (!ev1e2mu) && (!ev2e1mu)) continue;

    //deltaR condition
    if (!passDeltaRWleptZlept(WZcandidates, analysisLepton)) continue;
    numW++;  

    for (int fill1=0; fill1<4; fill1++){
      if (!ev[fill1]) continue;
      hZmass1[fill1]->Fill(Zmass);
      hMET1[fill1]->Fill(met);
      hZpt_my1[fill1]->Fill(Zpt);
      //      hLeadingJetPt_my1[fill1]->Fill(leadingRecoJetPt);
    }    


    if (ev3e){
      num3e++;
    }
    if (ev2e1mu){
      num2e1mu++;
    }
    if (ev1e2mu){
      num1e2mu++;
    }

    if (ev3mu){
      num3mu++;
    }
    //////////////////////////////////////MET CUT//////////////////////////
    
    //    if ((cWZ->pfmet)<30) continue;
    if ((cWZ->pfmetTypeI)<30) continue;   ///CHANGE THIS
    //    PtZ=Zpt;
    for (int fill2=0; fill2<4; fill2++){
      if (!ev[fill2]) continue;
      hZmass2[fill2]->Fill(Zmass);
      hMET2[fill2]->Fill(met);
      hZpt_my2[fill2]->Fill(Zpt);
      //      hLeadingJetPt_my2[fill2]->Fill(leadingRecoJetPt);
    }    

    if (ev3e){
      numMET3e++;
      myfile3e<<cWZ->run<<":"<<cWZ->lumi<<":"<<cWZ->event<<":"<<cWZ->pfmet<<":"<<cWZ->pfmetTypeI<<":"<<std::endl;
      //      hZmassElWel3->Fill(massEl);
      forSenka1->cd();
      PtZ_eee=Zpt;
      wz_aTGC1->Fill();    
    }
    if (ev2e1mu){
      myfile2e1mu<<cWZ->run<<":"<<cWZ->lumi<<":"<<cWZ->event<<":"<<cWZ->pfmet<<":"<<cWZ->pfmetTypeI<<":"<<std::endl;
      numMET2e1mu++;
      //     hZmassElWmu3->Fill(massEl);
      forSenka2->cd();
      PtZ_eem=Zpt;
      wz_aTGC2->Fill();    

    }
    if (ev1e2mu){
      myfile1e2mu<<cWZ->run<<":"<<cWZ->lumi<<":"<<cWZ->event<<":"<<cWZ->pfmet<<":"<<cWZ->pfmetTypeI<<":"<<std::endl;
      numMET1e2mu++;
      //      hZmassMuWel3->Fill(massMu);
      forSenka3->cd();
      PtZ_emm=Zpt;
      wz_aTGC3->Fill();    
    }

    if (ev3mu){
      myfile3mu<<cWZ->run<<":"<<cWZ->lumi<<":"<<cWZ->event<<":"<<cWZ->pfmet<<":"<<cWZ->pfmetTypeI<<":"<<std::endl;
      numMET3mu++;
      //      hZmassMuWmu3->Fill(massMu);
      PtZ_mmm=Zpt;
      forSenka4->cd();
      wz_aTGC4->Fill();    
    }
        
    myfileAll<<cWZ->run<<":"<<cWZ->lumi<<":"<<cWZ->event<<":"<<std::endl;

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
  

  

  fout->cd();
  fout->Write();
  fout->Close();


  forSenka1->cd();
  wz_aTGC1->Write();
  forSenka1->Close();
  

  forSenka2->cd();
  wz_aTGC2->Write();
  forSenka2->Close();
  

  forSenka3->cd();
  wz_aTGC3->Write();
  forSenka3->Close();
  

  forSenka4->cd();
  wz_aTGC4->Write();
  forSenka4->Close();
  
  if (writeOutputNumbers){
    fileNumData<<"#define dN_data3e "<<numMET3e<<std::endl;
    fileNumData<<"#define dN_data2e1mu "<<numMET2e1mu<<std::endl;
    fileNumData<<"#define dN_data1e2mu "<<numMET1e2mu<<std::endl;
    fileNumData<<"#define dN_data3mu "<<numMET3mu<<std::endl;
    fileNumData.close();
    }
}
