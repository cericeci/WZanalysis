#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

// Replace this with the new tree
#include "WZEvent.h"
//#include "WZGenEvent.h"
#include "UnfoldingHistogramFactory.h"
#include "HistogramFactory.h"
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
bool Z_muons(WZGenEvent *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, TLorentzVector* analysisLepton, float * ch,double & massMu, double & Zpt);

//bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, TLorentzVector* analysisLepton);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax=-999.);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);
double ScaleFactors(TH2F* MuonSF, TH2F* ElecSF,int type, int * WZCandidates, TLorentzVector* analysisLepton, double syst=0.0);

float trigger3sameLeptons(int* eL, int* eT);

float trigger2sameLeptons(int* eL, int* eT);

float TriggerWeight(int* WZcandidates, TH2F* DoubleElLead, TH2F* DoubleMuLead, TH2F* DoubleElTrail, TH2F* DoubleMuTrail, int type,TLorentzVector * analysisLepton);

void readFileFromList(TString fileList, std::vector<TString>* inputFile);

float AxeError(float weight, TH2F * MuonSF, TH2F * ElecSF, float pileUpWeight, int * WZCandidates,int type, TLorentzVector * analysisLepton, float TriggerEff, float triggerError=0);

double ReturnBranchingWeight(int type);

int determineGenType(WZEvent * cWZ);

TH1F * GetHistogramFromGraph(TString hname, TString gname);

TLorentzVector GetMET(Float_t metModule, Float_t metPhi);

float GetError(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);

//int genInformation(WZ *cWZ, int* Wdecay, int* Zdecay);

int main()
{
  using namespace std;

  //write output numbers
  //  bool writeOutputNumbers(false);
  bool writeOutputNumbers(true);

  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;
  ofstream fileNumGEN; 
  if (writeOutputNumbers){
    fileNumGEN.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numGEN_met.h");
    //    fileNumGEN.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numGEN.h");
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

  TFile * fout= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WZ.root", "RECREATE");
  //  TFile * forSenka= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/wz-8TeV.root", "RECREATE");
  // !!! TFile * forInv= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/scaleFact.root", "RECREATE");
  //TTree* wz_aTGC = new TTree("tgcTree", "tgcTree");
  //TTree* wz_scale = new TTree("scaleFactors", "scaleFactors");
  float gen_channel(-999), gen_mZ(0), gen_ptZ(0), reco_channel(-999), reco_mZ(0), reco_ptZ(0), efficiency_weight(0), xs_weight(0), pu_weight(0), br_weight(0); 
  //creating tree for investigation
  /*
  float muonScaleF(0), eleScaleF(0), muonError(0), eleError(0), muonPt(0), elePt(0), muonEta(0), eleEta(0);
  wz_scale->Branch("muonScaleF", &muonScaleF, "muonScaleF/F");
  wz_scale->Branch("eleScaleF", &eleScaleF, "eleScaleF/F");
  wz_scale->Branch("muonError", &muonError, "muonError/F");
  wz_scale->Branch("eleError", &eleError, "eleError/F");
  wz_scale->Branch("muonPt", &muonPt, "muonPt/F");
  wz_scale->Branch("elePt", &elePt, "elePt/F");
  wz_scale->Branch("muonEta", &muonEta, "muonEta/F");
  wz_scale->Branch("eleEta", &eleEta, "eleEta/F");
  */
  //creating tree for Senka
  /*
  wz_aTGC->Branch("gen_channel", &gen_channel, "gen_channel/F");
  wz_aTGC->Branch("gen_mZ", &gen_mZ, "gen_mZ/F");
  wz_aTGC->Branch("gen_ptZ", &gen_ptZ, "gen_ptZ/F");
  wz_aTGC->Branch("reco_channel", &reco_channel, "reco_channel/F");
  wz_aTGC->Branch("reco_mZ", &reco_mZ, "reco_mZ/F");
  wz_aTGC->Branch("reco_ptZ", &reco_ptZ, "reco_ptZ/F");
  wz_aTGC->Branch("efficiency_weight", &efficiency_weight, "efficiency_weight/F");
  wz_aTGC->Branch("pu_weight", &pu_weight, "pu_weight/F");
  wz_aTGC->Branch("xs_weight", &xs_weight, "xs_weight/F");
  wz_aTGC->Branch("br_weight", &br_weight, "br_weight/F");
  */

  const int nChannels(4);
  const int nChannels1(4);
  
  TH1D * hZptCh[nChannels];
  TH1D * hZptTau[nChannels];
  TH1D * hZptTauFraction[nChannels];
  TH1D * hLeadingJetCh[nChannels]; 
  TH1D * hLeadingJetTau[nChannels]; 
  TH1D * hLeadingJetTauFraction[nChannels]; 

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
  

  for (int hist=0; hist<nChannels; hist++){
    std::ostringstream ZptChname,ZptTauname, LeadingJetChname,LeadingJetTau, LeadingJetTauFraction, ZptTauFraction;
    LeadingJetTau<<"LeadingJetTau_"<<(hist+1);
    LeadingJetChname<<"LeadingJetCh_"<<(hist+1);
    LeadingJetTauFraction<<"LeadingJetPt_"<<(hist+1);
    ZptChname<<"ZptCh_"<<(hist+1);
    ZptTauname<<"ZptTau_"<<(hist+1);
    ZptTauFraction<<"Zpt_"<<(hist+1);
    
    hZptCh[hist]     = UnfoldingHistogramFactory::createZPtHistogram(ZptChname.str().c_str(), ZptChname.str().c_str());
    hZptTau[hist]     = UnfoldingHistogramFactory::createZPtHistogram(ZptTauname.str().c_str(), ZptTauname.str().c_str());
    hZptTauFraction[hist]     = UnfoldingHistogramFactory::createZPtHistogram(ZptTauFraction.str().c_str(),ZptTauFraction.str().c_str());
    hLeadingJetCh[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetChname.str().c_str(),LeadingJetChname.str().c_str());
    hLeadingJetTau[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetTau.str().c_str(),LeadingJetTau.str().c_str());
    hLeadingJetTauFraction[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetTauFraction.str().c_str(),LeadingJetTauFraction.str().c_str());
  }
  
  TH2F* MuonSF;
  TH2F* ElecSF;
  TH2F* DoubleElLead;
  TH2F* DoubleElTrail;
  TH2F* DoubleMuLead;
  TH2F* DoubleMuTrail;
  TH1F* hScaleInEB;
  TH1F* hScaleOutEB;
  TH1F* hScaleEE;



  MuonSF=LoadHistogram("auxiliaryFiles/MuSF_2012.root", "h2inverted", "MuonSF");
  ElecSF=LoadHistogram("auxiliaryFiles/EleSF_2012.root", "h2inverted", "ElecSF");
  DoubleElLead=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleElLead", "DoubleElLead");
  DoubleMuLead=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleMuLead", "DoubleMuLead");
  DoubleElTrail=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleElTrail", "DoubleElTrail");
  DoubleMuTrail=LoadHistogram("auxiliaryFiles/triggerEfficiencies.root", "DoubleMuTrail", "DoubleMuTrail");

  hScaleInEB  = GetHistogramFromGraph("hScaleInEB",  "gScaleInEB");
  hScaleOutEB = GetHistogramFromGraph("hScaleOutEB", "gScaleOutEB");
  hScaleEE    = GetHistogramFromGraph("hScaleEE",    "gScaleEE");

  const int leptonNumber(4);
  const double electronMass(0.000511);
  const double muonMass(0.106);
  const double luminosity(19602);
  double numZ(0), numW(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0), numMET2e1mu(0), numMET1e2mu(0), numMET3mu(0),numMET3e_1(0), numMET2e1mu_1(0), numMET1e2mu_1(0), numMET3mu_1(0), numMET3e_brCorr(0), numMET2e1mu_brCorr(0), numMET1e2mu_brCorr(0), numMET3mu_brCorr(0);
  double errorAxe3e(0), errorAxe2e1mu(0), errorAxe1e2mu(0), errorAxe3mu(0);
  double numMET3eGEN(0), numMET2e1muGEN(0), numMET1e2muGEN(0), numMET3muGEN(0), numMETSomethingGEN(0);
  double numMET3eGEN_brCorr(0), numMET2e1muGEN_brCorr(0), numMET1e2muGEN_brCorr(0), numMET3muGEN_brCorr(0), numMETSomethingGEN_brCorr(0);
  double numMET3e_2(0), numMET2e1mu_2(0), numMET1e2mu_2(0), numMET3mu_2(0);
  int numGEN3e_test(0);

  double errorTau[4],numMET[4], numMET_brCorr[4], numMET_brCorr2[4], numMET_1[4], numTau[4], numTau2[4], numMET_2[4], errorAxe[4], numMET_test[4], numTau_test[4], numMET2_test[4];
  int outOfZmassWindow[4];

  double numTau3e(0), numTau2e1mu(0), numTau1e2mu(0), numTau3mu(0);
  //0=EEE, 1=EEM, 2=EMM, 3= MMM, 4=all
  for (int ini=0; ini<4; ini++){
    numMET[ini]=0;
    numMET_brCorr[ini]=0;
    numMET_brCorr2[ini]=0;
    numMET_1[ini]=0;
    numTau[ini]=0;
    numTau2[ini]=0;
    numMET_2[ini]=0;
    errorAxe[ini]=0;
    outOfZmassWindow[ini]=0;
    numMET_test[ini]=0;
    numTau_test[ini]=0;
    numMET2_test[ini]=0;
  }  

double gen[5], genq[5], fact[5], factq[5], mixed[5], error[5], errorWholeZagreb[5], errorRest[5], errorSpanish[5];
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
  double numberPileupDown(0.0);
  int numGenChannels(9);
  //type: 0-(eee), 1-(eem), 2-(mme), 3-(mmm), 4-(tte), 5-(ttm), 6-(ttt), 7-(eet), 8-(mmt), 9-all(denominator)
  int numbersGen[10]={0,0,0,0,0,0,0,0,0,0};
  std::vector<TString>inputName;
  TChain wz("latino");

  readFileFromList("wzNewPu.files", &inputName);
  //    readFileFromList("wzNoSkim.files", &inputName);
  //  readFileFromList("wzNoFilter10.files", &inputName);

  for (int input=0; input< inputName.size(); input++){
    wz.Add(inputName[input]);
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZEvent *cWZ= new WZEvent(wz_tTree);
  //  WZGenEvent *cWZ= new WZGenEvent(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;
  
  for  (Int_t k = 0; k<events /*&& k<2000*/;k++) {
    wz_tTree->GetEntry(k);
    
    cWZ->ReadEvent();
    //*******various factors
    float pileUpWeight=cWZ->puW_new;
    //    float pileUpWeight=1;
    //float pileUpWeight=1;
    
    numberPileUp+=pileUpWeight;
    numberPileupDown++;
    
    xs_weight= (1.058* luminosity)/1.82755e+06;    //this is pileup weighted number of events
    //xs_weight= (1.058* luminosity)/1.82798e+06; //this is number of events   

    //    std::cout<<"XS WEIGHT: "<<xs_weight<<std::endl;
    ///////////////GEN YIELDS/////////////////
 
    //    double gen_channels[3][3]={{0,0,0}, {0,0,0}};
    bool isTau(false), is3e(false), is2e1mu(false), is1e2mu(false), is3mu(false);
    bool genKind[5]={false,false, false,false};
    int typeGen=-9999;


   
 
   
    //wz    std::cout<<typeGen<<" : "<<weightBr<<std::endl;
    
    double inZmassWindow(false);
    //    double weightBr=1;
    if (((cWZ->MZ)>71.1876) && ((cWZ->MZ<111.1876)))
      {
	inZmassWindow=true;
	typeGen= determineGenType(cWZ);
      }
    double weightBr(1);
    
    if (typeGen>=0){
      numbersGen[typeGen]++;
      numbersGen[9]++;
      weightBr=ReturnBranchingWeight(typeGen);
    }    
    else{
      int typeGenOut=determineGenType(cWZ);
      weightBr=ReturnBranchingWeight(typeGenOut);
    }

    
    //    weightBr=1;
    //float weightBr2=cWZ->GetBrWeight();
    //    std::cout<<weightBr<<" =? "<<weightBr2<<std::endl;
     if ((cWZ->WZchan)== 0){
       numGEN3e_test++;
     }
    //WZchan== 0/eee, 1/eemu, 2/emumu, 3/mumumu, 4/tau
    if (((cWZ->MZ)>71.1876) && ((cWZ->MZ<111.1876)))
      {
	//	typeGen= determineGenType(cWZ);
      if ((cWZ->WZchan)== 0){
	numMET3eGEN+=pileUpWeight;
	numMET3eGEN_brCorr+=pileUpWeight*weightBr;
       gen[0]+=pileUpWeight;
       genq[0]+=(pileUpWeight*pileUpWeight);
       //       typeGen=0;
       is3e=true;
      }
      if ((cWZ->WZchan)== 1){
	numMET2e1muGEN+=pileUpWeight;
	numMET2e1muGEN_brCorr+=pileUpWeight*weightBr;
	gen[1]+=pileUpWeight;
	genq[1]+=(pileUpWeight*pileUpWeight);
	is2e1mu=true;
      }
      if ((cWZ->WZchan)== 2){
	numMET1e2muGEN+=pileUpWeight;
	numMET1e2muGEN_brCorr+=pileUpWeight*weightBr;
      	gen[2]+=pileUpWeight;
	genq[2]+=(pileUpWeight*pileUpWeight);
	is1e2mu=true;
      }
      if ((cWZ->WZchan)== 3){
	numMET3muGEN+=pileUpWeight;
	numMET3muGEN_brCorr+=pileUpWeight*weightBr;
	gen[3]+=pileUpWeight;
	genq[3]+=(pileUpWeight*pileUpWeight);
	is3mu=true;
      }
      if ((cWZ->WZchan)==4){
	numMETSomethingGEN+=pileUpWeight;
	numMETSomethingGEN_brCorr+=pileUpWeight*weightBr;
	gen[4]+=pileUpWeight;
	genq[4]+=(pileUpWeight*pileUpWeight);

      }
	      
    }
    if ((cWZ->WZchan)==4){
      isTau=true;
    }
    int vrsta=cWZ->WZchan;
    genKind[vrsta]=true;
	
    bool kindOfEvent[5]={false, false, false, false};
    for (int kind=0; kind<5; kind++){
      if ((cWZ->WZchan)==kind){
	kindOfEvent[kind]=true;
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
    
    //here goes upgrade to TLorentz vector for analysis leptons
    TLorentzVector analysisLepton[leptonNumber];
    TLorentzVector analysisLeptonOld[leptonNumber];
    TLorentzVector EventMET;
    //    std::vector<TLorentzVector> analysisLepton;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////scale syst:///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool muScaleSyst(false); 
    bool elScaleSyst(false);
    double muScale(0.002);
    double elScale(1.0);

    double pfmet=cWZ->pfmetTypeI;
    double pfmetphi=cWZ->pfmetTypeIphi;
    EventMET = GetMET(pfmet, pfmetphi);

    for (int i1=0; i1<leptonNumber; i1++){
      if ((fabs(pdgid[i1])==13)&& (pt[i1]>0)){
	analysisLepton[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], muonMass);
	analysisLeptonOld[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], muonMass);
	if (muScaleSyst){
	  double spt=pt[i1]+ pt[i1]*muScale;
	  double factScale=pt[i1]/spt;
	  analysisLepton[i1]*=factScale;
	  EventMET +=(analysisLepton[i1]-analysisLeptonOld[i1]);
	}
      }
      if ((fabs(pdgid[i1])==11) && pt[i1]>0){
	analysisLepton[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], electronMass);
      	analysisLeptonOld[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], electronMass);
	if (elScaleSyst){
	  double scale;
	  const Float_t InEBMax  = hScaleInEB ->GetXaxis()->GetBinCenter(hScaleInEB ->GetNbinsX());
	  const Float_t OutEBMax = hScaleOutEB->GetXaxis()->GetBinCenter(hScaleOutEB->GetNbinsX());
	  const Float_t EEMax    = hScaleEE   ->GetXaxis()->GetBinCenter(hScaleEE   ->GetNbinsX());
	  const Float_t scaleInEB  = hScaleInEB ->GetBinContent(hScaleInEB ->FindBin(min(pt[i1], InEBMax)));
	  const Float_t scaleOutEB = hScaleOutEB->GetBinContent(hScaleOutEB->FindBin(min(pt[i1], OutEBMax)));
	  const Float_t scaleEE    = hScaleEE   ->GetBinContent(hScaleEE   ->FindBin(min(pt[i1], EEMax)));


	  const Float_t aeta = fabs(eta[i1]);

	  if (aeta < 0.8)
	    {
	      scale = scaleInEB;
	    }
	  else if (aeta >= 0.8 && aeta < 1.479)
	    {
	      scale = scaleOutEB;
	    }
	  else
	    {
	      scale = scaleEE;
	    }
	  double spt=pt[i1]+ pt[i1]*scale*elScale;
	  double factScale=pt[i1]/spt;
	  analysisLepton[i1]*=factScale;
	}
	EventMET +=(analysisLepton[i1]-analysisLeptonOld[i1]);
      }	
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    int WZcandidates[3]; 
    
    //here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
    int lepNum(0);
    for (int i=0; i<leptonNumber; i++){
      if ((analysisLepton[i].Pt()>10) && (analysisLepton[i].Pt()!=-9999))
	lepNum++;
      if (((fabs(pdgid[i]))==11) && (analysisLepton[i].Pt()>10) && (!pass2012ICHEP[i]))
	testEl++;
      if ((fabs(pdgid[i])==13) && (analysisLepton[i].Pt()>10) && (!pass2012ICHEP[i]))
	testMu++;
      if (((fabs(pdgid[i]))==11) && (analysisLepton[i].Pt()>10) && (pass2012ICHEP[i])){
	good_electrons.push_back(i);
	v_nizEl[i].SetPtEtaPhiM(analysisLepton[i].Pt(),analysisLepton[i].Eta(), analysisLepton[i].Phi(), electronMass);
	v_3Lepton=v_3Lepton+v_nizEl[i];
      }

      if ((fabs(pdgid[i])==13) && (analysisLepton[i].Pt()>10) && (pass2012ICHEP[i])){
      good_muons.push_back(i);
      v_nizMu[i].SetPtEtaPhiM(analysisLepton[i].Pt(),analysisLepton[i].Eta(), analysisLepton[i].Phi(), muonMass);
	v_3Lepton=v_3Lepton+v_nizMu[i];
      }
    }
    if (lepNum!=3) continue;

    if (v_3Lepton.M()<100) continue;

    bool foundZel(false), foundZmu(false);
    double massMu(-999), massEl(0), Zpt(0);

    foundZmu=  Z_muons(cWZ, &good_muons, WZcandidates, v_nizMu, analysisLepton, ch, massMu, Zpt);
    if (foundZmu){
      //    hZmassMu1->Fill(massMu);
    }

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
    
    double Zmass;
    if (foundZel) Zmass=massEl;
    else Zmass=massMu;


    numZ+=pileUpWeight;

    ///////////////////////JETS//////////////////////////////////
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
	if (analysisLepton[elIndex1].Pt()>20){
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
	if (analysisLepton[muIndex1].Pt()>20){
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
	if (analysisLepton[elIndex2].Pt()>20)
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
	if (analysisLepton[muIndex2].Pt()>20)
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
    
    double weights[4]={0,0,0,0};
    

    if ((!ev3e) && (!ev3mu) && (!ev1e2mu) && (!ev2e1mu)) continue;

    //deltaR condition
    if (!passDeltaRWleptZlept(WZcandidates, analysisLepton)) continue;
    numW+=pileUpWeight;  

    for (int fill1=0; fill1<4; fill1++){
      //weights
      weights[fill1]=pileUpWeight*ScaleFactors(MuonSF, ElecSF, fill1, WZcandidates, analysisLepton)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, fill1, analysisLepton)*xs_weight;
      if (!ev[fill1]) continue;
      hZmass1[fill1]->Fill(Zmass, weights[fill1]);
      hMET1[fill1]->Fill(EventMET.Et(), weights[fill1]);
      hZpt_my1[fill1]->Fill(Zpt, weights[fill1]);
      hLeadingJetPt_my1[fill1]->Fill(leadingRecoJetPt, weights[fill1]);
    }


    
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
    
    //    if ((cWZ->pfmet)<30) continue;
    //    if ((cWZ->pfmetTypeI)<30) continue;   ///CHANGE THIS
    if ((EventMET.Et())<30) continue;   ///CHANGE THIS

    double weight;
    //-systematics:
    //    double syst=1.0;
    double syst=0.0;
    for (int ch=0; ch<4; ch++){
      if (ev[ch]){
	reco_channel=ch;
	///////////////////////////////////part for testing scale factors////////////////////////////////////////////////////////////////////////////////////
	for (int f=0; f<3; f++){
	  int index=WZcandidates[f];
	  double mass=analysisLepton[index].M();
	  if (mass=electronMass){
	    double factor= GetFactor(ElecSF, analysisLepton[index].Pt(), analysisLepton[index].Eta());
	    double error= GetError(ElecSF, analysisLepton[index].Pt(), analysisLepton[index].Eta());
	    /*
	    eleScaleF=factor;
	    eleError=error/factor;
	    elePt= analysisLepton[index].Pt();
	    eleEta=analysisLepton[index].Eta();
	    hScaleEl->Fill(error/factor,analysisLepton[index].Eta());
	    */
	    }
	  if (mass=muonMass){
	    double factor= GetFactor(MuonSF, analysisLepton[index].Pt(), analysisLepton[index].Eta());
	    double error= GetError(MuonSF, analysisLepton[index].Pt(), analysisLepton[index].Eta());
	    /*
	    muonScaleF=factor;
	    muonError=error/factor;
	    muonPt= analysisLepton[index].Pt();
	    muonEta=analysisLepton[index].Eta();
	    hScaleMu->Fill(error/factor,analysisLepton[index].Eta());
	    */
	    }
	  //wz_scale->Fill();
	}
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, ch, WZcandidates, analysisLepton, syst)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, ch, analysisLepton);
	efficiency_weight=ScaleFactors(MuonSF, ElecSF, ch, WZcandidates, analysisLepton, syst)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, ch, analysisLepton);
	numMET_test[ch]+=weightBr;
	numMET[ch]+=weight;
	numMET_brCorr[ch]+=weight*weightBr;
	numMET_brCorr2[ch]+=weight*weight*weightBr*weightBr;
	numMET_1[ch]++;
	if (!inZmassWindow){
	  outOfZmassWindow[ch]++;
	}
	if (isTau) 
	  {
	    numTau[ch]+=weight*weightBr;
	    numTau2[ch]+=weight*weight*weightBr*weightBr;
	    numTau_test[ch]+=weightBr;
	  }
	if (genKind[ch]) {
	  numMET_2[ch]+=weight;
	  numMET2_test[ch]+=weightBr;
	  fact[ch]+=weight;
	  factq[ch]+=(weight*weight);
	  mixed[ch]+=(weight*pileUpWeight);
      }
	float triggerEff=TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, ch, analysisLepton);
	errorAxe[ch]+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, ch, analysisLepton, triggerEff);
	errorRest[ch]+=AxeError(weight, MuonSF, ElecSF, pileUpWeight, WZcandidates, ch, analysisLepton, triggerEff);
      //      myfile3e<<cWZ->run<<":"<<cWZ->event<<":"<<pileUpWeight<<":"<<ScaleFactors(MuonSF, ElecSF, 0, WZcandidates, pt, eta)<<":"<<triggerEff<<":"<<std::endl;      
      }
    }
    //filling aTGC tree
    gen_channel=cWZ->WZchan;
    gen_mZ= cWZ->MZ;
    gen_ptZ=cWZ->PtZ;
    reco_mZ=Zmass;
    reco_ptZ=Zpt;
    pu_weight=cWZ->puW_new;
    //xs_weight:

    br_weight= weightBr;
    //    wz_aTGC->Fill();

    for (int fill2=0; fill2<4; fill2++){
      if (!ev[fill2]) continue;
      hZmass2[fill2]->Fill(Zmass, weights[fill2]);
      hMET2[fill2]->Fill(EventMET.Et(), weights[fill2]);
      hZpt_my2[fill2]->Fill(Zpt, weights[fill2]);
      hLeadingJetPt_my2[fill2]->Fill(leadingRecoJetPt, weights[fill2]);
    }
    

    


    //filling histograms
     
    for (int filH=0; filH<nChannels; filH++){
      if (ev[filH]){
	hZptCh[filH]->Fill(Zpt, weight);
	hLeadingJetCh[filH]->Fill(leadingRecoJetPt, weight);
	if (isTau){
	  hZptTau[filH]->Fill(Zpt, weight);
	  hLeadingJetTau[filH]->Fill(leadingRecoJetPt, weight);
	}
      }
      
    }
  }
  

  //dividing histograms:

  for (int divide=0; divide<nChannels; divide++){
    hZptTau[divide]->Sumw2();
    hZptCh[divide]->Sumw2();
    hZptTauFraction[divide]->Divide(hZptTau[divide],hZptCh[divide],1,1,"B");

    hLeadingJetTau[divide]->Sumw2();
    hLeadingJetCh[divide]->Sumw2();
    hLeadingJetTauFraction[divide]->Divide(hLeadingJetTau[divide],hLeadingJetCh[divide],1,1,"B");

  }

  float denominator=numMET3eGEN_brCorr+numMET2e1muGEN_brCorr+numMET1e2muGEN_brCorr+numMET3muGEN_brCorr+numMETSomethingGEN_brCorr;
  gen[5]=gen[0]+gen[1]+gen[2]+gen[3]+gen[4];

  for (int e=0; e<4; e++){
    std::cout<<errorAxe[e]<<std::endl;
    error[e]=sqrt((gen[e]*gen[e]*factq[e]-2*gen[e]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[e])/pow(gen[e],4));
    //    errorWholeZagreb[e]=sqrt(errorRest[e]/(gen[e]*gen[e])+(gen[e]*gen[e]*factq[e]-2*gen[e]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[e])/pow(gen[e],4));
    errorWholeZagreb[e]=sqrt((errorRest[e]/pow(gen[e],4))+(gen[e]*gen[e]*factq[e]-2*gen[e]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[e])/pow(gen[e],4));
    errorSpanish[e]=sqrt((gen[4]*gen[4]*factq[e]-2*gen[4]*fact[e]*mixed[e]+fact[e]*fact[e]*genq[4])/pow(gen[4],4));
    errorTau[e]=sqrt((numMET_brCorr[e]*numMET_brCorr[e]*numTau2[e]- 2*numMET_brCorr[e]*numTau[e]*numTau[e] + numTau[e]*numTau[e]*numMET_brCorr2[e])/(pow(numMET_brCorr[e],4)));
    
  }

  
  ///*********OUTPUT
  //some tests:
  std::cout<<"FOR NORMALIZATION: "<<numberPileUp<<std::endl<<" , "<<numberPileupDown<<std::endl;
  std::cout<<"3e: "<< numMET_2[0]<<std::endl;
  std::cout<<"2e1mu: "<<numMET_2[1]<<std::endl;
  std::cout<<"1e2mu: "<<numMET_2[2]<<std::endl;
  std::cout<<"3mu: "<<numMET_2[3]<<std::endl;
  std::cout<<"***************"<<std::endl;;

  /*  std::cout<<"3e: "<< numMET3eGEN<<std::endl;
  std::cout<<"2e1mu: "<<numMET2e1muGEN<<std::endl;
  std::cout<<"1e2mu: "<<numMET1e2muGEN<<std::endl;
  std::cout<<"3mu: "<<numMET3muGEN<<std::endl;
  */
  std::cout<<"Normalized: "<<std::endl;
  std::cout<<"3e: "<< numMET_2[0]*xs_weight<<std::endl;
  std::cout<<"2e1mu: "<<numMET_2[1]*xs_weight<<std::endl;
  std::cout<<"1e2mu: "<<numMET_2[2]*xs_weight<<std::endl;
  std::cout<<"3mu: "<<numMET_2[3]*xs_weight<<std::endl;

  std::cout<<"NumMET_brCorr"<<std::endl;
  std::cout<<"3e: "<< numMET_brCorr[0]*xs_weight<<std::endl;
  std::cout<<"2e1mu: "<<numMET_brCorr[1]*xs_weight<<std::endl;
  std::cout<<"1e2mu: "<<numMET_brCorr[2]*xs_weight<<std::endl;
  std::cout<<"3mu: "<<numMET_brCorr[3]*xs_weight<<std::endl;

  std::cout<<"without corr"<<std::endl;
  std::cout<<"3e: "<< numMET[0]*xs_weight<<std::endl;
  std::cout<<"2e1mu: "<<numMET[1]*xs_weight<<std::endl;
  std::cout<<"1e2mu: "<<numMET[2]*xs_weight<<std::endl;
  std::cout<<"3mu: "<<numMET[3]*xs_weight<<std::endl;



  std::cout<<"Acceptance efficiency:"<<std::endl;
  std::cout<<"Zagreb: "<<std::endl;
  std::cout<<"3e: "<< numMET_2[0]/numMET3eGEN<<"+/-"<<error[0]<<" or "<<errorWholeZagreb[0]<<std::endl;
  std::cout<<"2e1mu: "<<numMET_2[1]/numMET2e1muGEN<<"+/-"<<error[1]<<" or "<<errorWholeZagreb[1]<<std::endl;
  std::cout<<"1e2mu: "<<numMET_2[2]/numMET1e2muGEN<<"+/-"<<error[2]<<" or "<<errorWholeZagreb[2]<<std::endl;
  std::cout<<"3mu: "<<numMET_2[3]/numMET3muGEN<<"+/-"<<error[3]<<" or "<<errorWholeZagreb[3]<<std::endl;
  std::cout<<"************"<<std::endl;
  /*
  std::cout<<"Spanish: "<<std::endl;
  std::cout<<"3e: "<< numMET_brCorr[0]/denominator<<std::endl;
  std::cout<<"2e1mu: "<<numMET_brCorr[1]/denominator<<std::endl;
  std::cout<<"1e2mu:"<<numMET_brCorr[2]/denominator<<std::endl;
  std::cout<<"3mu: "<<numMET_brCorr[3]/denominator<<std::endl;
  */

  std::cout<<numMET_brCorr[0]/denominator<<std::endl;
//writing in file
  if (writeOutputNumbers){
    //spnish numbers
    fileNumGEN<<"#define dAxe3eS "<< numMET_brCorr[0]/denominator<<std::endl;
    fileNumGEN<<"#define dAxe2e1muS "<<numMET_brCorr[1]/denominator<<std::endl;
    fileNumGEN<<"#define dAxe1e2muS "<<numMET_brCorr[2]/denominator<<std::endl;
    fileNumGEN<<"#define dAxe3muS "<<numMET_brCorr[3]/denominator<<std::endl;
    

    fileNumGEN<<"#define dAxe3eZ_2 "<< numMET_2[0]/numMET3eGEN<<std::endl;
    fileNumGEN<<"#define dAxe2e1muZ_2 "<<numMET_2[1]/numMET2e1muGEN<<std::endl;
    fileNumGEN<<"#define dAxe1e2muZ_2 "<<numMET_2[2]/numMET1e2muGEN<<std::endl;
    fileNumGEN<<"#define dAxe3muZ_2 "<<numMET_2[3]/numMET3muGEN<<std::endl;
    
    fileNumGEN<<"#define dtauFactor3e "<<numTau[0]/numMET_brCorr[0]<<std::endl;
    fileNumGEN<<"#define dtauFactor2e1mu "<<numTau[1]/numMET_brCorr[1]<<std::endl;
    fileNumGEN<<"#define dtauFactor1e2mu "<<numTau[2]/numMET_brCorr[2]<<std::endl;
    fileNumGEN<<"#define dtauFactor3mu "<<numTau[3]/numMET_brCorr[3]<<std::endl;
  
    fileNumGEN<<"#define dsAxe3eS "<< errorSpanish[0]<<std::endl;
    fileNumGEN<<"#define dsAxe2e1muS "<<errorSpanish[1]<<std::endl;
    fileNumGEN<<"#define dsAxe1e2muS "<<errorSpanish[2]<<std::endl;
    fileNumGEN<<"#define dsAxe3muS "<<errorSpanish[3]<<std::endl;
    
    fileNumGEN<<"#define dsAxe3eZ "<< error[0]<<std::endl;
    fileNumGEN<<"#define dsAxe2e1muZ "<<error[1]<<std::endl;
    fileNumGEN<<"#define dsAxe1e2muZ "<<error[2]<<std::endl;
    fileNumGEN<<"#define dsAxe3muZ "<<error[3]<<std::endl;
    
    fileNumGEN<<"#define dstauFactor3e "<<errorTau[0]<<std::endl;
    fileNumGEN<<"#define dstauFactor2e1mu "<<errorTau[1]<<std::endl;
    fileNumGEN<<"#define dstauFactor1e2mu "<<errorTau[2]<<std::endl;
    fileNumGEN<<"#define dstauFactor3mu "<<errorTau[3]<<std::endl;
  
  

  fileNumGEN.close();
  }
    
  fout->cd();
  fout->Write();
  fout->Close();
  
  /*
  forInv->cd();
  wz_scale->Write();
  forInv->Close();
  */
  ////part for TGC
  /*
  forSenka->cd();
  wz_aTGC->Write();
  forSenka->Close();
  */
}
