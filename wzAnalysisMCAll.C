#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "UnfoldingHistogramFactory.h"
#include "HistogramFactory.h"
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

bool Z_muons(WZ *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, TLorentzVector* analysisLepton, float * ch,double & massMu, double & Zpt);

bool passMVAiso(float isomva, float pt, float eta);

bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz);

bool passDeltaRWleptZlept(int * WZcandidates, TLorentzVector * analysisLepton);

float GetFactor(TH2F* h2, float leptonPt, float leptonEta);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

double ScaleFactors(TH2F* MuonSF, TH2F* ElecSF,int type, int * WZCandidates, TLorentzVector* analysisLepton, double syst=0.0);

float trigger3sameLeptons(int* eL, int* eT);

float trigger2sameLeptons(int* eL, int* eT);

float TriggerWeight(int* WZcandidates, TH2F* DoubleElLead, TH2F* DoubleMuLead, TH2F* DoubleElTrail, TH2F* DoubleMuTrail, int type,TLorentzVector* analysisLepton);

void readFileFromList(TString fileList,std::vector<TString> * inputFile);

TH1F * GetHistogramFromGraph(TString hname, TString gname);

TLorentzVector GetMET(Float_t metModule, Float_t metPhi);

double deltaPhiWMET(int * WZcandidates, TLorentzVector* analysisLepton, TLorentzVector EventMET);

double wTransverseMass(int index, TLorentzVector* analysisLepton, TLorentzVector EventMET);

int main()
{
  using namespace std;
  
  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;
  ofstream fileNumMC;
  //***write output numbers
  //  bool writeOutputNumbers(false);
  bool writeOutputNumbers(true);

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

  std::vector<TString>name;
  std::vector<TChain*> chain;
  std::vector<TString> files;
  std::vector<double> errorMC;
  
  errorMC.push_back(0);
  errorMC.push_back(0.15);
  errorMC.push_back(0.15);
  errorMC.push_back(1);
  errorMC.push_back(0.2);
  errorMC.push_back(0);
  errorMC.push_back(0);
  
  files.push_back("WZ.files");
  files.push_back("ZZpu.files");
  files.push_back("ZgammaPu.files");
  files.push_back("WVpu.files");
  files.push_back("VVVpu.files");
  files.push_back("top.files");
  files.push_back("Zjets.files");

  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WZ.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/ZZ.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zgamma.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WV.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/VVV.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/top.root");
  name.push_back("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zjets.root");



  for (int file=1; file<(files.size()-2) /*&& file<1*/; file++){

    float numZ(0), numW(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0.0), numMET2e1mu(0.0), numMET1e2mu(0.0), numMET3mu(0.0);
    float numMET3ectrl(0), numMET2e1muctrl(0), numMET1e2muctrl(0), numMET3muctrl(0);
    float numMET3eGEN(0), numMET2e1muGEN(0), numMET1e2muGEN(0), numMET3muGEN(0);
    float error3e(0), error2e1mu(0), error1e2mu(0), error3mu(0);
    TFile * fout = new TFile(name[file],"RECREATE");
    int numTest(0);
    double numMET[4], error[4];
    
    //initialization:
    for (int ini=0; ini<4; ini++){
      numMET[ini]=0;
      error[ini]=0;
    }
    //old (erase?)

    const int nChannels(4);
    const int nChannels1(5);
    
    TH1D * hZpt[nChannels];
    TH1D * hLeadingJetPt[nChannels]; 
    TH1D * hZptAll[nChannels];
    TH1D * hZptAllError[nChannels];
    TH1D * hZptAll3sigmaUp[nChannels];
    TH1D * hZptAll3sigmaDown[nChannels];

    //naming convention: 0-after Z candidate, 1-after W candidate 2- after whole selection    
    TH1F * hZmass1[nChannels1];
    TH1F * hZmass2[nChannels1];
    TH1F * hMET1[nChannels1];
    TH1F * hMET2[nChannels1];
    TH1F * hZpt_my1[nChannels1];
    TH1F * hZpt_my2[nChannels1];
    TH1F * hLeadingJetPt_my1[nChannels1];
    TH1F * hLeadingJetPt_my2[nChannels1];
    TH1F * hNjets1[nChannels1];
    TH1F * hNjets2[nChannels1];
    TH1F * hDeltaPhi1[nChannels1];
    TH1F * hDeltaPhi2[nChannels1];
    TH1F * hZlepton1pt1[nChannels1];
    TH1F * hZlepton2pt1[nChannels1];
    TH1F * hZlepton1pt2[nChannels1]; 
    TH1F * hZlepton2pt2[nChannels1]; 
    TH1F * hWleptonpt1[nChannels1];
    TH1F * hWleptonpt2[nChannels1];    
    TH1F * hMTW1[nChannels1];
    TH1F * hMTW2[nChannels1];
    
    for (int myhist=0; myhist<nChannels1; myhist++){
      std::ostringstream nZmass1, nMET1, nZpt1, nLeadingJetPt1, nNjets1, nDeltaPhi1, nZlepton1pt1, nZlepton2pt1, nWleptonpt1, nMTW1;
      std::ostringstream nZmass2, nMET2, nZpt2, nLeadingJetPt2, nNjets2, nDeltaPhi2, nZlepton1pt2, nZlepton2pt2, nWleptonpt2, nMTW2;
      nZmass1<<"hZmass1_"<<myhist;
      nZmass2<<"hZmass2_"<<myhist;
      nMET1<<"hMET1_"<<myhist;
      nMET2<<"hMET2_"<<myhist;
      nZpt1<<"hZpt1_"<<myhist;
      nZpt2<<"hZpt2_"<<myhist;
      nLeadingJetPt1<<"hLeadingJetPt1_"<<myhist;
      nLeadingJetPt2<<"hLeadingJetPt2_"<<myhist;
      nNjets1<<"hNjets1_"<<myhist;
      nNjets2<<"hNjets2_"<<myhist;
      nDeltaPhi1<<"hDeltaPhi1_"<<myhist;
      nDeltaPhi2<<"hDeltaPhi2_"<<myhist;
      nZlepton1pt1<<"hZlepton1pt1_"<<myhist;
      nZlepton2pt1<<"hZlepton2pt1_"<<myhist;
      nZlepton1pt2<<"hZlepton1pt2_"<<myhist;
      nZlepton2pt2<<"hZlepton2pt2_"<<myhist;
      nWleptonpt1<<"hWleptonpt1_"<<myhist;
      nWleptonpt2<<"hWleptonpt2_"<<myhist;
      nMTW1<<"hMTW1_"<<myhist;
      nMTW2<<"hMTW2_"<<myhist;

    
      hZmass1[myhist]           = HistogramFactory::createZmassHisto(nZmass1.str().c_str(), nZmass1.str().c_str());
      hZmass2[myhist]           = HistogramFactory::createZmassHisto(nZmass2.str().c_str(), nZmass2.str().c_str());
      hMET1[myhist]             = HistogramFactory::createMETHisto(nMET1.str().c_str(), nMET1.str().c_str());
      hMET2[myhist]             = HistogramFactory::createMETHisto(nMET2.str().c_str(), nMET2.str().c_str());
      hZpt_my1[myhist]          = HistogramFactory::createZptHisto(nZpt1.str().c_str(), nZpt1.str().c_str());
      hZpt_my2[myhist]          = HistogramFactory::createZptHisto(nZpt2.str().c_str(), nZpt2.str().c_str());
      hLeadingJetPt_my1[myhist] = HistogramFactory::createLeadingJetptHisto(nLeadingJetPt1.str().c_str(), nLeadingJetPt1.str().c_str());
      hLeadingJetPt_my2[myhist] = HistogramFactory::createLeadingJetptHisto(nLeadingJetPt2.str().c_str(), nLeadingJetPt2.str().c_str());
      hNjets1[myhist]           =  HistogramFactory::createNjetsHisto(nNjets1.str().c_str(), nNjets1.str().c_str());
      hNjets2[myhist]           =  HistogramFactory::createNjetsHisto(nNjets2.str().c_str(), nNjets2.str().c_str());
      hDeltaPhi1[myhist]        =  HistogramFactory::createDeltaPhi(nDeltaPhi1.str().c_str(), nDeltaPhi1.str().c_str());
      hDeltaPhi2[myhist]        =  HistogramFactory::createDeltaPhi(nDeltaPhi2.str().c_str(), nDeltaPhi2.str().c_str());
      hZlepton1pt1[myhist]       = HistogramFactory::createZptHisto(nZlepton1pt1.str().c_str(), nZlepton1pt1.str().c_str());
      hZlepton2pt1[myhist]       = HistogramFactory::createZptHisto(nZlepton2pt1.str().c_str(), nZlepton2pt1.str().c_str());
      hZlepton1pt2[myhist]       = HistogramFactory::createZptHisto(nZlepton1pt2.str().c_str(), nZlepton1pt2.str().c_str());
      hZlepton2pt2[myhist]       = HistogramFactory::createZptHisto(nZlepton2pt2.str().c_str(), nZlepton2pt2.str().c_str());
      hWleptonpt1[myhist]       = HistogramFactory::createZptHisto(nWleptonpt1.str().c_str(), nWleptonpt1.str().c_str());
      hWleptonpt2[myhist]       = HistogramFactory::createZptHisto(nWleptonpt2.str().c_str(), nWleptonpt2.str().c_str());
      hMTW1[myhist]              = HistogramFactory::createMTW(nMTW1.str().c_str(), nMTW1.str().c_str());
      hMTW2[myhist]              = HistogramFactory::createMTW(nMTW2.str().c_str(), nMTW2.str().c_str());
    }

    for (int hist=0; hist<nChannels; hist++){
      std::ostringstream Zptname, LeadingJetname, ZptfSname, ZptfSnameError, ZptfSname3Up, ZptfSname3Down ;
      LeadingJetname<<"LeadingJetPt_"<<(hist+1);
      Zptname<<"Zpt_"<<(hist+1);
      ZptfSname<<"total_bkg_rebined_"<<(hist);
      ZptfSnameError<<"total_bkg_rebined_error_"<<(hist);
      ZptfSname3Up<<"total_bkg_rebined_3sigmaUp_"<<(hist);
      ZptfSname3Down<<"total_bkg_rebined_3sigmaDown_"<<(hist);

      hZpt[hist]     = UnfoldingHistogramFactory::createZPtHistogram(Zptname.str().c_str(), Zptname.str().c_str());
      hLeadingJetPt[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetname.str().c_str(),LeadingJetname.str().c_str());
      hZptAll[hist ]=UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptfSname.str().c_str(), ZptfSname.str().c_str());
      hZptAllError[hist ]=UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptfSnameError.str().c_str(), ZptfSnameError.str().c_str());
      hZptAll3sigmaUp[hist ]=UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptfSname3Up.str().c_str(), ZptfSname3Up.str().c_str());
      hZptAll3sigmaDown[hist ]=UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptfSname3Down.str().c_str(), ZptfSname3Down.str().c_str());
      
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

    


    for  (Int_t k = 0; k<events /*&& k<1000*/;k++) {
      wz_tTree->GetEntry(k);    
      cWZ->ReadEvent();
      
      float xs_weight(0);    

      //*******various factors
      float pileUpWeight=cWZ->puW_new;
      
      //float pileUpWeight=1.0;
    
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
    
    TLorentzVector analysisLepton[leptonNumber];
    TLorentzVector analysisLeptonOld[leptonNumber];
    TLorentzVector EventMET;

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
    
    bool muScaleSyst(false);
    bool elScaleSyst(false);
    double muScale(0.002);
    double elScale(1.0);

    double pfmet=cWZ->pfmet;
    double pfmetphi=cWZ->pfmetphi;
    EventMET = GetMET(pfmet, pfmetphi);

    for (int i1=0; i1<leptonNumber; i1++){
      if ((bdt[i1]<100) && (pt[i1]>10)){
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
	  EventMET +=(analysisLepton[i1]-analysisLeptonOld[i1]);
	}
      }
      if ((bdt[i1]>100) && (pt[i1]>10)){
	analysisLepton[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], muonMass);
	analysisLeptonOld[i1].SetPtEtaPhiM(pt[i1], eta[i1], phi[i1], muonMass);
	if (muScaleSyst){
	  double spt=pt[i1]+ pt[i1]*muScale;
	  double factScale=pt[i1]/spt;
	  analysisLepton[i1]*=factScale;
	  EventMET +=(analysisLepton[i1]-analysisLeptonOld[i1]);
	}
      }
    }
//here goes index of WZ candidate: 1. first Z, 2. second Z, 3. W lepton
    
    int lepNum(0);
    for (int i=0; i<leptonNumber; i++){
      if ((analysisLepton[i].Pt()>10) && (analysisLepton[i].Pt()!=-9999))
	lepNum++;
      if ((bdt[i]<100) && (pass2012ICHEP[i]) && (analysisLepton[i].Pt()>10)){
	good_electrons.push_back(i);
	v_nizEl[i].SetPtEtaPhiM(analysisLepton[i].Pt(),analysisLepton[i].Eta(), analysisLepton[i].Phi(), electronMass);
	v_3Lepton=v_3Lepton+v_nizEl[i];
      }
      if ((bdt[i]>100) && (pass2012ICHEP[i]) && (analysisLepton[i].Pt()>10)){
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
    //    if (foundZmu){
    //hZmassMu1->Fill(massMu, pileUpWeight*xs_weight);
    //}

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

    numZ+=pileUpWeight;

    double Zmass;
    if (foundZel) Zmass=massEl;
    else Zmass=massMu;
    
    /////////JETS//////////////////////////////////
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
    
    /////////////FILLING HISTOGRAMS//////////////////////////////
    bool ev[5]={false, false, false, false, true};
    if (ev3e) ev[0]=true;
    if (ev2e1mu) ev[1]=true;
    if (ev1e2mu) ev[2]=true;
    if (ev3mu) ev[3]=true;
    
    double weights[4]={0,0,0,0};
    


    if ((!ev3e) && (!ev3mu) && (!ev1e2mu) && (!ev2e1mu)) continue;


    //deltaR condition
    if (!passDeltaRWleptZlept(WZcandidates, analysisLepton)) continue;
    numW+=pileUpWeight;  

    double deltaPhi=deltaPhiWMET(WZcandidates,analysisLepton, EventMET);
    int indW=WZcandidates[2];
    int indZ1=WZcandidates[0];
    int indZ2=WZcandidates[1];
    double WleptPt= pt[indW];
    double Zlept1Pt= pt[indZ1];
    double Zlept2Pt= pt[indZ2];
    double MTW= wTransverseMass(indW, analysisLepton,EventMET);
    double weights_h(0);

    for (int w=0; w<4; w++){
      weights[w]=pileUpWeight*ScaleFactors(MuonSF, ElecSF, w, WZcandidates, analysisLepton)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, w, analysisLepton)*xs_weight;
      weights_h=pileUpWeight*ScaleFactors(MuonSF, ElecSF, w, WZcandidates, analysisLepton)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, w, analysisLepton)*xs_weight;
    }

    /////////////FILLING HISTOGRAMS//////////////////////////////
    for (int fill1=0; fill1<5; fill1++){
      //weights

      if (!ev[fill1]) continue;
      hZmass1[fill1]->Fill(Zmass, weights_h);
      hMET1[fill1]->Fill(EventMET.Et(), weights_h);
      hZpt_my1[fill1]->Fill(Zpt, weights_h);
      hLeadingJetPt_my1[fill1]->Fill(leadingRecoJetPt, weights_h);
      hNjets1[fill1]->Fill(nRecoJets, weights_h);
      hDeltaPhi1[fill1]->Fill(deltaPhi, weights_h);
      hZlepton1pt1[fill1]->Fill(Zlept1Pt,weights_h);
      hZlepton2pt1[fill1]->Fill(Zlept2Pt,weights_h);
      hWleptonpt1[fill1]->Fill(WleptPt,weights_h);
      hMTW1[fill1]->Fill(MTW, weights_h);
    }
    


    //////////////////////////////////////MET CUT//////////////////////////
    
    if (EventMET.Et()<30) continue;
    //    if ((cWZ->pfmet)<30) continue;
    //    if ((cWZ->pfmetTypeI)<30) continue;   ///CHANGE THIS

    //improving this
    
    for (int iEv=0; iEv<4; iEv++){
      if (ev[iEv]){
	weight=pileUpWeight*ScaleFactors(MuonSF, ElecSF, iEv, WZcandidates, analysisLepton)*TriggerWeight(WZcandidates, DoubleElLead, DoubleMuLead, DoubleElTrail, DoubleMuTrail, iEv, analysisLepton)*xs_weight;
	numMET[iEv]+=weight;
	error[iEv]+=(weight*weight);
	//	numMETcrtl[iEv]+=pileUpWeight;
      }
    }

    for (int fill2=0; fill2<5; fill2++){
      if (!ev[fill2]) continue;
      hZmass2[fill2]->Fill(Zmass, weights_h);
      hMET2[fill2]->Fill(EventMET.Et(), weights_h);
      hZpt_my2[fill2]->Fill(Zpt, weights_h);
      hLeadingJetPt_my2[fill2]->Fill(leadingRecoJetPt, weights_h);
      hNjets2[fill2]->Fill(nRecoJets, weights_h);
      hDeltaPhi2[fill2]->Fill(deltaPhi, weights_h);
      hZlepton1pt2[fill2]->Fill(Zlept1Pt,weights_h);
      hZlepton2pt2[fill2]->Fill(Zlept2Pt,weights_h);
      hWleptonpt2[fill2]->Fill(WleptPt,weights_h);
      hMTW2[fill2]->Fill(MTW, weights_h);
    }    


    
    for (int filH=0; filH<nChannels; filH++){
      if (ev[filH]){
	//hZpt[filH]->Fill(Zpt, weight);
	hLeadingJetPt[filH]->Fill(leadingRecoJetPt, weight);
	hZptAll[filH]->Fill(Zpt, weight);
      }
    }


    } //end of loop over all events
    
    //filling other histos and errors
    
    for (int er=0; er< nChannels; er++){
      //double binError(0);
      for (int zpt2=0; zpt2<(hZptAll[er]->GetNbinsX()+1); zpt2++){
	double binContent=hZptAll[er]->GetBinContent(zpt2);	
	double binError=errorMC[file]*binContent;
	hZptAll[er]->SetBinError(zpt2, binError);
	hZptAll3sigmaUp[er]-> SetBinContent(zpt2, (binContent+3*binError));
	//	hZptAll3sigmaUp[er]->SetBinError(zpt2, binError);
	if ((binContent-3*binError)>0)
	  hZptAll3sigmaDown[er]-> SetBinContent(zpt2, (binContent-3*binError));
	else
	  hZptAll3sigmaDown[er]-> SetBinContent(zpt2, 0);
	//	hZptAll3sigmaDown[er]->SetBinError(zpt2, binError);
      }
  }
   
  std::cout<<"****************************(normalized to the luminosity, with scale factors and trigger eff)"<<std::endl;
  std::cout<<"3e:     "<<numMET[0]<<"+/-"<<sqrt(error[0])<<std::endl;
  std::cout<<"2e1mu:  "<<numMET[1]<<"+/-"<<sqrt(error[1])<<std::endl;
  std::cout<<"1e2mu:  "<<numMET[2]<<"+/-"<<sqrt(error[2])<<std::endl;
  std::cout<<"3mu:    "<<numMET[3]<<"+/-"<<sqrt(error[3])<<std::endl;
  std::cout<<"numTest:"<< numTest<<std::endl;
  
  fout->cd();
  fout->Write();
  
     /*
     for (int histDel=0; histDel<nChannels<4; histDel++){
       delete hZpt[histDel];
       delete hLeadingJetPt[histDel];
     }
     */
  /*
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
  */    
 fout->Close();    
    //writing in file:
    if (writeOutputNumbers){
      if (file==1){
	fileNumMC<<"#define dNZZ_3e "<<numMET[0]<<std::endl;
	fileNumMC<<"#define dNZZ_2e1mu "<<numMET[1]<<std::endl;
	fileNumMC<<"#define dNZZ_1e2mu "<<numMET[2]<<std::endl;
	fileNumMC<<"#define dNZZ_3mu "<<numMET[3]<<std::endl;
	fileNumMC<<"#define dsNZZ_3e "<<error[0]<<std::endl;
	fileNumMC<<"#define dsNZZ_2e1mu "<<error[1]<<std::endl;
	fileNumMC<<"#define dsNZZ_1e2mu "<<error[2]<<std::endl;
	fileNumMC<<"#define dsNZZ_3mu "<<error[3]<<std::endl;
      }
      
      if (file==2){
	fileNumMC<<"#define dNZgamma_3e "<<numMET[0]<<std::endl;
	fileNumMC<<"#define dNZgamma_2e1mu "<<numMET[1]<<std::endl;
	fileNumMC<<"#define dNZgamma_1e2mu "<<numMET[2]<<std::endl;
	fileNumMC<<"#define dNZgamma_3mu "<<numMET[3]<<std::endl;
	fileNumMC<<"#define dsNZgamma_3e "<<error[0]<<std::endl;
	fileNumMC<<"#define dsNZgamma_2e1mu "<<error[1]<<std::endl;
	fileNumMC<<"#define dsNZgamma_1e2mu "<<error[2]<<std::endl;
	fileNumMC<<"#define dsNZgamma_3mu "<<error[3]<<std::endl;
      }
      if (file==3){
	fileNumMC<<"#define dNWV_3e "<<numMET[0]<<std::endl;
	fileNumMC<<"#define dNWV_2e1mu "<<numMET[1]<<std::endl;
	fileNumMC<<"#define dNWV_1e2mu "<<numMET[2]<<std::endl;
	fileNumMC<<"#define dNWV_3mu "<<numMET[3]<<std::endl;
	fileNumMC<<"#define dsNWV_3e "<<error[0]<<std::endl;
	fileNumMC<<"#define dsNWV_2e1mu "<<error[1]<<std::endl;
	fileNumMC<<"#define dsNWV_1e2mu "<<error[2]<<std::endl;
	fileNumMC<<"#define dsNWV_3mu "<<error[3]<<std::endl;
      }
      if (file==4){
	fileNumMC<<"#define dNVVV_3e "<<numMET[0]<<std::endl;
	fileNumMC<<"#define dNVVV_2e1mu "<<numMET[1]<<std::endl;
	fileNumMC<<"#define dNVVV_1e2mu "<<numMET[2]<<std::endl;
	fileNumMC<<"#define dNVVV_3mu "<<numMET[3]<<std::endl;
	fileNumMC<<"#define dsNVVV_3e "<<error[0]<<std::endl;
	fileNumMC<<"#define dsNVVV_2e1mu "<<error[1]<<std::endl;
	fileNumMC<<"#define dsNVVV_1e2mu "<<error[2]<<std::endl;
	fileNumMC<<"#define dsNVVV_3mu "<<error[3]<<std::endl;
      }
    }
    
  }


  fileNumMC.close();

}
