#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
//#include "ThirdMuons.h"
#include "TString.h"
#include "constants.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>

// DATA is defined in the makefile for the case of Data analysis

#ifdef DATA
#define WZTREE WZ2012Data
#include "WZ2012Data.h"
// <<<<<<< HEAD
// #else
// #define WZTREE WZGenEvent
// #include "WZGenEvent.h"
// =======
#endif
#ifdef OLDMC
#define WZTREE WZ
#include "WZ.h"
// >>>>>>> 70fa253f99ebb78c063c3656e7c573b6bc6f7b8b
#endif
#ifdef NEWMC
//#define WZTREE WZEvent
//#include "WZEvent.h"
#define WZTREE WZGenEvent
#include "WZGenEvent.h"
#endif

#define WZTEST WZEvent
#include "WZEvent.h"

//#include "WZ.h"
void readChainFromList(TString fileList, TChain * chain)
{
  const int linesize=1024;
  
  std::ifstream  inputChain(fileList);
  char inbufChain[linesize];
  
  while (inputChain.getline(inbufChain,linesize)) {
    if (inbufChain[0] == '#') {
      continue;
    }
    chain->Add(inbufChain);
  }
}

void readFileFromList(TString fileList, std::vector<TString>* inputFile)
{
  const int linesize=1024;
  
  std::ifstream  inputChain(fileList);
  char inbufChain[linesize];
  
  while (inputChain.getline(inbufChain,linesize)) {
    if (inbufChain[0] == '#') {
      continue;
    }
    inputFile->push_back(inbufChain);
  }
}

bool Z_muons(WZTREE *cWZ, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz, TLorentzVector* analysisLepton, float * ch,double & massMu, double & Zpt)
{
  TLorentzVector v;  
  float diff= 500;
  double mass1=0;
  double ZtransIm=0;
  std::vector<float> mass;
  std::vector<int> Z_muons;
  int Z_muonIndex1=-100, Z_muonIndex2=-100;
  float Z_muon_charge1=0, Z_muon_charge2=0;
  bool found_Z=false;

  //computing invariant mass
  for (int imu1=0; imu1<good_muons->size(); imu1++) {
    for (int imu2=imu1+1; imu2<good_muons->size(); imu2++) {
      int mu1=(*good_muons)[imu1];
      int mu2=(*good_muons)[imu2];
      float pt1= analysisLepton[mu1].Pt();
      float pt2= analysisLepton[mu2].Pt();
      v= v_niz[mu1]+v_niz[mu2];

      mass.push_back(v.M());
      if (((fabs(v.M()-91.1876))<diff) && ((ch[mu1])*(ch[mu2])==-1)) {
      	diff=(fabs(v.M()-91.1876));
	mass1= v.M();
	ZtransIm=v.Pt();
	Z_muonIndex1= (*good_muons)[imu1];
	Z_muonIndex2= (*good_muons)[imu2];
	Z_muon_charge1=ch[Z_muonIndex1]; 
       	Z_muon_charge2=ch[Z_muonIndex2]; 
      }
    }
  }


  if ((Z_muon_charge1==0) || (Z_muon_charge2==0)) return false;
  //conditions for Z

  //saving index of Z muons
  if (((diff<20)) && (Z_muon_charge1 != 0) 
      && (Z_muon_charge2 != 0)
      && ((-1)*Z_muon_charge1== Z_muon_charge2)
      && (((analysisLepton[Z_muonIndex1].Pt()>20 ) && (analysisLepton[Z_muonIndex2].Pt())>10 ) 
	  || ((analysisLepton[Z_muonIndex2].Pt()>20 ) && ((analysisLepton[Z_muonIndex1].Pt())>10 )))) 
    {
      //      if (diff>20) std::cout<<diff<<std::endl;
      massMu=mass1;
      Zpt=ZtransIm;
      //first goes lepton with bigger pt
      if (analysisLepton[Z_muonIndex1].Pt()>analysisLepton[Z_muonIndex2].Pt()){
	WZcandidates[0] = Z_muonIndex1;
	WZcandidates[1] = Z_muonIndex2;
      }
      else{
	WZcandidates[0] = Z_muonIndex2;
	WZcandidates[1] = Z_muonIndex1;
      }
      return true;
    }   
  else 
    {
      return false;
    }
}



bool passMVAiso(float isomva, float pt, float eta){
  bool pass(false);
  if (pt<=20){
    if (fabs(eta)<1.479){
      if (isomva>=0.86) pass=true; 
    }
    else{
      if (isomva>=0.82) pass=true;
    }
  }
  if (pt>20){
    if (fabs(eta)<1.479){
      if (isomva>=0.82) pass=true;
    }
    else{
      if (isomva>=0.86) pass=true;
    }
  }
  return pass;
}


bool Z_independent(float * ch, std::vector<int>* good_muons,int * WZcandidates, TLorentzVector *v_niz)
{
  TLorentzVector v;  
  float diff = 500;
  float mass1=0;
  std::vector<float> mass;
  std::vector<float> differences;
  std::vector<int> index1;
  std::vector<int> index2;
  std::vector<int> Z_muons;
  int Z_muonIndex1=-100, Z_muonIndex2=-100;
  float Z_muon_charge1=0, Z_muon_charge2=0;
  bool foundZ(false);
  
  if (good_muons->size()<4) return false;
  else{
    for (int imu1=0; imu1<good_muons->size(); imu1++) {
      for (int imu2=imu1+1; imu2<good_muons->size(); imu2++) {
	v= v_niz[(*good_muons)[imu1]]+v_niz[(*good_muons)[imu2]];
	Z_muon_charge1=ch[(*good_muons)[imu1]]; 
	Z_muon_charge2=ch[(*good_muons)[imu2]]; 
	diff=fabs(v.M()-91.1876);
	if ((diff<20) && ((-1)*Z_muon_charge1== Z_muon_charge2)){
	mass.push_back(v.M());
	differences.push_back(diff);
	index1.push_back(imu1);
	index2.push_back(imu2);
     
	}
      }
    }
    int counter(0);
    for (int m=0; m< mass.size(); m++){
      for (int n=m+1; n<mass.size(); n++){
	if (index1[n]!= index1[m] &&
	    index1[n]!= index2[m] &&
	    index2[n]!= index2[m] &&
	    index2[n]!= index1[m]){
	  counter++; 
	}
      }
    }
    
    if (counter!=0){
      //std::cout<<"Number of independent Z muon candidates in one event: "<<counter<<std::endl;
      return true;
    }
    else return false;
  }
}

bool passDeltaRWleptZlept(int * WZcandidates, TLorentzVector* analysisLepton){
  int z1=WZcandidates[0];
  int z2=WZcandidates[1];
  int w=WZcandidates[2];
  float phi1=analysisLepton[z1].Phi();
  float phi2=analysisLepton[z2].Phi();
  float phi3=analysisLepton[w].Phi();
  float eta1=analysisLepton[z1].Eta();
  float eta2=analysisLepton[z2].Eta();
  float eta3=analysisLepton[w].Eta();
  float deltaPhiZ1W=acos(cos(phi1-phi3));
  float deltaEtaZ1W= eta1-eta3;
  float deltaPhiZ2W=acos(cos(phi2-phi3));
  float deltaEtaZ2W=eta2-eta3;
  float deltaRZ1W= sqrt(pow(deltaPhiZ1W,2)+pow(deltaEtaZ1W,2));
  float deltaRZ2W= sqrt(pow(deltaPhiZ2W,2)+pow(deltaEtaZ2W,2));
  if ((deltaRZ1W<0.1) || (deltaRZ2W<0.1))
    return false;
  else return true;
}

double deltaPhiWMET(int * WZcandidates, TLorentzVector* analysisLepton, TLorentzVector EventMET){
  int w=WZcandidates[2];
  float phiW=analysisLepton[w].Phi();
  float phiMET=EventMET.Phi();
  float deltaPhiWMET=acos(cos(phiW-phiMET));
  return deltaPhiWMET;
}

double wTransverseMass(int index, TLorentzVector* analysisLepton, TLorentzVector EventMET)
{
  float metPhi= EventMET.Phi();
  float phi=analysisLepton[index].Phi();
  float deltaPhi=acos(cos(metPhi-phi));
  float met=EventMET.Et();
  float pt=analysisLepton[index].Pt(); 
  float tm=sqrt(2*pt*met*(1-TMath::Cos(deltaPhi)));
  return tm;
}



float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.){

  float aeta=fabs(leptonEta);
  int nbins= h2->GetNbinsX();
  float ptMax = (leptonPtMax > 0) ? leptonPtMax : h2->GetXaxis()->GetBinCenter(nbins);
  //  float ptMax= h2->GetXaxis()->GetBinCenter(nbins);
  float factor= h2->GetBinContent(h2->FindBin(TMath::Min(leptonPt, ptMax), aeta));
  return factor;
}

float GetError(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.){

  float aeta=fabs(leptonEta);
  int nbins= h2->GetNbinsX();
  float ptMax = (leptonPtMax > 0) ? leptonPtMax : h2->GetXaxis()->GetBinCenter(nbins);
  //  float ptMax= h2->GetXaxis()->GetBinCenter(nbins);
  float factor= h2->GetBinError(h2->FindBin(TMath::Min(leptonPt, ptMax), aeta));
  return factor;
}

TH2F* LoadHistogram(TString filename, TString hname, TString cname){
  
  TFile * inputFile= TFile::Open(filename);
  TH2F* hist= (TH2F*)inputFile->Get(hname)->Clone(cname);
  
  hist->SetDirectory(0);
  inputFile->Close();
  
  return hist;
}

double ScaleFactors(TH2F* MuonSF, TH2F* ElecSF,int type, int * WZCandidates, TLorentzVector* analysisLepton, double syst=0.0){
  int index1=WZCandidates[0];
  int index2=WZCandidates[1];
  int index3=WZCandidates[2];
  double factor1(-999), factor2(-999), factor3(-999);
  double errorEle(0), errorMu(0);
  syst=0;
  //double error1(0), error2(0), error3(0);
  
  //  std::cout<<pt[index1]<<":"<<pt[index2]<<":"<<pt[index3]<<std::endl;
  //  std::cout<<syst<<std::endl;
  if (type==0) {
    //    error1= GetError(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    //error2= GetError(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    //error3= GetError(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    //    factor1=1000;
    //factor2=1000;
    //factor3=1000;
    factor1= GetFactor(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta()) + syst*errorEle;
    factor2= GetFactor(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta()) + syst*errorEle;
    factor3= GetFactor(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta()) + syst*errorEle;
  }

  if (type==1){
    //error1= GetError(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    //error2= GetError(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    //error3= GetError(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());

    factor1= GetFactor(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta()) + syst*errorEle;
    factor2= GetFactor(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta()) + syst*errorEle;
    factor3= GetFactor(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta()) + syst*errorMu;
  }
  if (type==2){
    //error1= GetError(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    //error2= GetError(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    //    error3= GetError(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    
    factor1= GetFactor(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta()) + syst*errorMu;
    factor2= GetFactor(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta()) + syst*errorMu;
    factor3= GetFactor(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta()) + syst*errorEle;
  }
  if (type==3){
    //    error1= GetError(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    //error2= GetError(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    //error3= GetError(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());

    factor1= GetFactor(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta()) + syst*errorMu;
    factor2= GetFactor(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta()) + syst*errorMu;
    factor3= GetFactor(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta()) + syst*errorMu;
    //    std::cout<<factor1<<":"<<factor2<<":"<<factor3<<std::endl;
    // std::cout<<error1<<":"<<error2<<":"<<error3<<std::endl;  
    //std::cout<<"*********"<<std::endl;
  }
  
  //  std::cout<<"********************"<<std::endl;
  //  std::cout<<factor1*factor2*factor3<<std::endl;
  //  double factor=factor1*factor2*factor3;

  return (factor1*factor2*factor3);
}

float trigger3sameLeptons(float* eL, float* eT){
  //  float r1= 1-eL[0]*eL[1]*eL[2];
  float r1= (1-eL[0])*(1-eL[1])*(1-eL[2]);
  float r2= eL[0]*(1-eT[1])*(1-eT[2]);
  float r3= eL[1]*(1-eT[0])*(1-eT[2]);
  float r4= eL[2]*(1-eT[0])*(1-eT[1]);
  return (1-(r1+r2+r3+r4));
}

float trigger2sameLeptons(float* eL, float* eT){
  //  float r1=1-eL[0]*eL[1];
  float r1=(1-eL[0])*(1-eL[1]);
  float r2=eL[0]*(1-eT[1]);
  float r3=eL[1]*(1-eT[0]);
  return (1-(r1+r2+r3));
}

// Cross trigger efficiency: flavour of 3rd lepton is trailing
float triggerDifferentLeptons(float* eL, float* eT){
  //  float r1=1-eL[0]*eL[1];
  float r1=(1-eT[2]);
  float r2=eT[2]*(1-eL[0])*(1-eL[1]);
  return (1-(r1+r2));
}


float TriggerWeight(int* WZcandidates, TH2F* DoubleElLead, TH2F* DoubleMuLead, TH2F* DoubleElTrail, TH2F* DoubleMuTrail, int type,TLorentzVector * analysisLepton){
  int index1=WZcandidates[0];
  int index2=WZcandidates[1];
  int index3=WZcandidates[2];
  float factor;
  float eL[3]={0,0,0};
  float eT[3]={0,0,0};
  if (type==0){
    eL[0]=GetFactor(DoubleElLead, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eL[1]=GetFactor(DoubleElLead, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eL[2]=GetFactor(DoubleElLead, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    eT[0]=GetFactor(DoubleElTrail, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eT[1]=GetFactor(DoubleElTrail, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eT[2]=GetFactor(DoubleElTrail, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    factor=trigger3sameLeptons(eL, eT);
    
  }
  if (type==1){
    //Jonatan's way
    eL[0]=GetFactor(DoubleElLead, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eL[1]=GetFactor(DoubleElLead, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eL[2]=GetFactor(DoubleMuLead, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    eT[0]=GetFactor(DoubleElTrail, analysisLepton[index1].Pt(),analysisLepton[index1].Eta());
    eT[1]=GetFactor(DoubleElTrail, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eT[2]=GetFactor(DoubleMuTrail, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    
    //    factor=trigger2sameLeptons(eL, eT);
    factor=trigger3sameLeptons(eL, eT);
  }
  if (type==2){
    //Jonatan's way
    eL[0]=GetFactor(DoubleMuLead, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eL[1]=GetFactor(DoubleMuLead, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eL[2]=GetFactor(DoubleElLead, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    eT[0]=GetFactor(DoubleMuTrail, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eT[1]=GetFactor(DoubleMuTrail, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eT[2]=GetFactor(DoubleElTrail, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    
    //factor=trigger2sameLeptons(eL, eT);
    factor=trigger3sameLeptons(eL, eT);
  }
  if (type==3){
    eL[0]=GetFactor(DoubleMuLead, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eL[1]=GetFactor(DoubleMuLead, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eL[2]=GetFactor(DoubleMuLead, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    eT[0]=GetFactor(DoubleMuTrail, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    eT[1]=GetFactor(DoubleMuTrail, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    eT[2]=GetFactor(DoubleMuTrail, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    factor=trigger3sameLeptons(eL, eT);
  }
  return factor;
}

double weight(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label)
{
  int index1=WZcandidates[0];
  int index2=WZcandidates[1];
  int index3=WZcandidates[2];
  float pt1=pt[index1];
  float pt2=pt[index2];
  float pt3=pt[index3];
  float eta1=eta[index1];
  float eta2=eta[index2];
  float eta3=eta[index3];
  float p1(0), p2(0), p3(0), epsilon1(0), epsilon2(0), epsilon3(0);
  float ep1(0), ep2(0), ep3(0), eepsilon1(0), eepsilon2(0), eepsilon3(0);

  //eee
  if (type==0){
    p1=GetFactor(ElectronFR, pt1, eta1);
    ep1=GetError(ElectronFR, pt1, eta1);
    p2=GetFactor(ElectronFR, pt2, eta2);
    ep2=GetError(ElectronFR, pt2, eta2);
    p3=GetFactor(ElectronFR, pt3, eta3);
    ep3=GetError(ElectronFR, pt3, eta3);
    epsilon1=GetFactor(ElectronPR, pt1, eta1);
    //    eepsilon1=GetError(ElectronPR, pt1, eta1);
    epsilon2=GetFactor(ElectronPR, pt2, eta2);
    //    eepsilon2=GetError(ElectronPR, pt2, eta2);
    epsilon3=GetFactor(ElectronPR, pt3, eta3);
    //    eepsilon3=GetError(ElectronPR, pt3, eta3);
  }
  //eem
  if (type==1){
    p1=GetFactor(ElectronFR, pt1, eta1);
    ep1=GetError(ElectronFR, pt1, eta1);
    p2=GetFactor(ElectronFR, pt2, eta2);
    ep2=GetError(ElectronFR, pt2, eta2);
    p3=GetFactor(MuonFR, pt3, eta3, 34.);
    ep3=GetError(MuonFR, pt3, eta3, 34.);
    epsilon1=GetFactor(ElectronPR, pt1, eta1);
    //    eepsilon1=GetError(ElectronPR, pt1, eta1);
    epsilon2=GetFactor(ElectronPR, pt2, eta2);
    //    eepsilon2=GetError(ElectronPR, pt2, eta2);
    epsilon3=GetFactor(MuonPR, pt3, eta3);
    //    eepsilon3=GetError(MuonPR, pt3, eta3);
  }
  //emm
  if (type==2){
    p1=GetFactor(MuonFR, pt1, eta1, 34.);
    ep1=GetError(MuonFR, pt1, eta1, 34.);
    p2=GetFactor(MuonFR, pt2, eta2, 34.);
    ep2=GetError(MuonFR, pt2, eta2, 34.);
    p3=GetFactor(ElectronFR, pt3, eta3);
    ep3=GetError(ElectronFR, pt3, eta3);
    epsilon1=GetFactor(MuonPR, pt1, eta1);
    //    eepsilon1=GetError(MuonPR, pt1, eta1);
    epsilon2=GetFactor(MuonPR, pt2, eta2);
    //    eepsilon2=GetError(MuonPR, pt2, eta2);
    epsilon3=GetFactor(ElectronPR, pt3, eta3);
    //    eepsilon3=GetError(ElectronPR, pt3, eta3);
  }
  //mmm
  if (type==3){
    p1=GetFactor(MuonFR, pt1, eta1, 34.);
    ep1=GetError(MuonFR, pt1, eta1, 34.);
    p2=GetFactor(MuonFR, pt2, eta2, 34.);
    ep2=GetError(MuonFR, pt2, eta2, 34.);
    p3=GetFactor(MuonFR, pt3, eta3, 34.);
    ep3=GetError(MuonFR, pt3, eta3, 34.);
    epsilon1=GetFactor(MuonPR, pt1, eta1);
    //    eepsilon1=GetError(MuonPR, pt1, eta1);
    epsilon2=GetFactor(MuonPR, pt2, eta2);
    //    eepsilon2=GetError(MuonPR, pt2, eta2);
    epsilon3=GetFactor(MuonPR, pt3, eta3);
    //    eepsilon3=GetError(MuonPR, pt3, eta3);
  }
  
  double Afactor=1 / ((epsilon1-p1)*(epsilon2-p2)*(epsilon3-p3));
  
  double Bfactor=epsilon1*epsilon2*epsilon3;
  
  double Cfactor(0);
  if (label==1) Cfactor=(1-p1)*(1-p2)*(1-p3);
  if (label==2) Cfactor=-(1-p1)*(1-p2)*p3;
  if (label==3) Cfactor=-(1-p1)*p2*(1-p3);
  if (label==4) Cfactor=(1-p1)*p2*p3;
  if (label==5) Cfactor=-p1*(1-p2)*(1-p3);
  if (label==6) Cfactor=p1*(1-p2)*p3;
  if (label==7) Cfactor=p1*p2*(1-p3);
  if (label==8) Cfactor=p1*p2*p3;
  
  double w= Afactor*Bfactor*Cfactor;
  
  return w;
}

int determineLabel(int * pass2012ICHEP, int * WZcandidates){
  //1= TTT
  //2= TTF
  //3= TFT
  //4= TFF
  //5= FTT
  //6= FTF
  //7= FFT
  //8= FFF
  int index1=WZcandidates[0];
  int index2=WZcandidates[1];
  int index3=WZcandidates[2];
  if ((pass2012ICHEP[index1]==true)  && (pass2012ICHEP[index2]==true)  && (pass2012ICHEP[index3]==true))  return 1;
  if ((pass2012ICHEP[index1]==true)  && (pass2012ICHEP[index2]==true)  && (pass2012ICHEP[index3]==false)) return 2;
  if ((pass2012ICHEP[index1]==true)  && (pass2012ICHEP[index2]==false) && (pass2012ICHEP[index3]==true))  return 3;
  if ((pass2012ICHEP[index1]==true)  && (pass2012ICHEP[index2]==false) && (pass2012ICHEP[index3]==false)) return 4;
  if ((pass2012ICHEP[index1]==false) && (pass2012ICHEP[index2]==true)  && (pass2012ICHEP[index3]==true))  return 5;
  if ((pass2012ICHEP[index1]==false) && (pass2012ICHEP[index2]==true)  && (pass2012ICHEP[index3]==false)) return 6;
  if ((pass2012ICHEP[index1]==false) && (pass2012ICHEP[index2]==false) && (pass2012ICHEP[index3]==true))  return 7;
  if ((pass2012ICHEP[index1]==false) && (pass2012ICHEP[index2]==false) && (pass2012ICHEP[index3]==false)) return 8;
}

float MMerror(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label, float weight){
  int index1=WZcandidates[0];
  int index2=WZcandidates[1];
  int index3=WZcandidates[2];
  float pt1=pt[index1];
  float pt2=pt[index2];
  float pt3=pt[index3];
  float eta1=eta[index1];
  float eta2=eta[index2];
  float eta3=eta[index3];
  float p1(0), p2(0), p3(0), e1(0), e2(0), e3(0);
  float ee1(0), ee2(0), ee3(0), ep1(0), ep2(0), ep3(0);
  double errorEvent(0);

  //eee
  if (type==0){
    p1=GetFactor(ElectronFR, pt1, eta1);
    ep1=GetError(ElectronFR, pt1, eta1);
    p2=GetFactor(ElectronFR, pt2, eta2);
    ep2=GetError(ElectronFR, pt2, eta2);
    p3=GetFactor(ElectronFR, pt3, eta3);
    ep3=GetError(ElectronFR, pt3, eta3);
    e1=GetFactor(ElectronPR, pt1, eta1);
    ee1=GetError(ElectronPR, pt1, eta1);
    e2=GetFactor(ElectronPR, pt2, eta2);
    ee2=GetError(ElectronPR, pt2, eta2);
    e3=GetFactor(ElectronPR, pt3, eta3);
    ee3=GetError(ElectronPR, pt3, eta3);
  }
  //eem
  if (type==1){
    p1=GetFactor(ElectronFR, pt1, eta1);
    ep1=GetError(ElectronFR, pt1, eta1);
    p2=GetFactor(ElectronFR, pt2, eta2);
    ep2=GetError(ElectronFR, pt2, eta2);
    p3=GetFactor(MuonFR, pt3, eta3, 34.);
    ep3=GetError(MuonFR, pt3, eta3, 34.);
    e1=GetFactor(ElectronPR, pt1, eta1);
    ee1=GetError(ElectronPR, pt1, eta1);
    e2=GetFactor(ElectronPR, pt2, eta2);
    ee2=GetError(ElectronPR, pt2, eta2);
    e3=GetFactor(MuonPR, pt3, eta3);
    ee3=GetError(MuonPR, pt3, eta3);
  }
  //emm
  if (type==2){
    p1=GetFactor(MuonFR, pt1, eta1, 34.);
    ep1=GetError(MuonFR, pt1, eta1, 34.);
    p2=GetFactor(MuonFR, pt2, eta2, 34.);
    ep2=GetError(MuonFR, pt2, eta2, 34.);
    p3=GetFactor(ElectronFR, pt3, eta3);
    ep3=GetError(ElectronFR, pt3, eta3);
    e1=GetFactor(MuonPR, pt1, eta1);
    ee1=GetError(MuonPR, pt1, eta1);
    e2=GetFactor(MuonPR, pt2, eta2);
    ee2=GetError(MuonPR, pt2, eta2);
    e3=GetFactor(ElectronPR, pt3, eta3);
    ee3=GetError(ElectronPR, pt3, eta3);
  }
  //mmm
  if (type==3){
    p1=GetFactor(MuonFR, pt1, eta1, 34.);
    ep1=GetError(MuonFR, pt1, eta1, 34.);
    p2=GetFactor(MuonFR, pt2, eta2, 34.);
    ep2=GetError(MuonFR, pt2, eta2, 34.);
    p3=GetFactor(MuonFR, pt3, eta3, 34.);
    ep3=GetError(MuonFR, pt3, eta3, 34.);
    e1=GetFactor(MuonPR, pt1, eta1);
    ee1=GetError(MuonPR, pt1, eta1);
    e2=GetFactor(MuonPR, pt2, eta2);
    ee2=GetError(MuonPR, pt2, eta2);
    e3=GetFactor(MuonPR, pt3, eta3);
    ee3=GetError(MuonPR, pt3, eta3);
  }

    //second part= epsilon derivation
  if (label==1){
    float dl1de1=(e2*e3*(-1+p1)*p1*(-1+p2)*(-1+p3))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl1de2=(e1*e3*(-1+p1)*(-1+p2)*p2*(-1+p3))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl1de3=(e1*e2*(-1+p1)*(-1+p2)*(-1+p3)*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl1dp1=-((-1+e1)*e1*e2*e3*(-1+p2)*(-1+p3))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl1dp2=-(e1*(-1+e2)*e2*e3*(-1+p1)*(-1+p3))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl1dp3=-(e1*e2*(-1+e3)*e3*(-1+p1)*(-1+p2))/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl1de1*ee1),2)+ pow((dl1de2*ee2),2)+pow((dl1de3*ee3),2)+pow((dl1dp1*ep1),2)+pow((dl1dp2*ep2),2)+ pow((dl1dp3*ep3),2);
  }
  if (label==2){
    float dl2de1=(e2*e3*(p1-1)*p1*(p2-1)*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl2de2=(e1*e3*(p1-1)*(p2-1)*p2*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl2de3=(e1*e2*(p1-1)*(p2-1)*p3*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl2dp1=(-(e1-1)*e1*e2*e3*(p2-1)*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl2dp2=(-e1*(e2-1)*e2*e3*(p1-1)*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl2dp3=(-e1*e2*e3*e3*(p1-1)*(p2-1))/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl2de1*ee1),2)+ pow((dl2de2*ee2),2)+pow((dl2de3*ee3),2)+pow((dl2dp1*ep1),2)+pow((dl2dp2*ep2),2)+ pow((dl2dp3*ep3),2);
  }
  if (label==3){
    float dl3de1=(e2*e3*(p1-1)*p1*p2*(p3-1))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl3de2=(e1*e3*(p1-1)*p2*p2*(p3-1))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl3de3=(e1*e2*(p1-1)*p2*(p3-1)*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl3dp1=(-(e1-p1)*e1*e2*e3*p2*(p3-1))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl3dp2=(-e1*e2*e2*e3*(p1-1)*(p3-1))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl3dp3=(-e1*e2*(e3-1)*e3*(p1-1)*p2)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl3de1*ee1),2)+ pow((dl3de2*ee2),2)+pow((dl3de3*ee3),2)+pow((dl3dp1*ep1),2)+pow((dl3dp2*ep2),2)+ pow((dl3dp3*ep3),2);
  }
  if (label==4){
    float dl4de1=(e2*e3*(p1-1)*p1*p2*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl4de2=(e1*e2*(p1-1)*p2*p2*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl4de3=(e1*e2*(p1-1)*p2*p3*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl4dp1=(-(e1-1)*e1*e2*e3*p2*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl4dp2=(-e1*e2*e2*e3*(p1-1)*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl4dp3=(-e1*e2*e3*e3*(p1-1)*p2)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl4de1*ee1),2)+ pow((dl4de2*ee2),2)+pow((dl4de3*ee3),2)+pow((dl4dp1*ep1),2)+pow((dl4dp2*ep2),2)+ pow((dl4dp3*ep3),2);
  }
  if (label==5){
    float dl5de1=(-e2*e3*p1*p1*(p2-1)*(p3-1))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl5de2=(-e1*e3*p1*(p2-1)*p2*(p3-1))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl5de3=(-e1*e2*p1*(p2-1)*(p3-1)*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl5dp1=(e1*e1*e2*e3*(p2-1)*(p3-1))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl5dp2=(e1*(e2-1)*e2*e3*p1*(p3-1))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl5dp3=(e1*e2*(e3-1)*e3*p1*(p2-1))/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl5de1*ee1),2)+ pow((dl5de2*ee2),2)+pow((dl5de3*ee3),2)+pow((dl5dp1*ep1),2)+pow((dl5dp2*ep2),2)+ pow((dl5dp3*ep3),2);
  }
  if (label==6){
    float dl6de1=(e2*e3*p1*p1*(p2-1)*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl6de2=(e1*e3*p1*(p2-1)*p2*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl6de3=(e1*e2*p1*(p2-1)*p3*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl6dp1=(-e1*e1*e2*e3*(p2-1)*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl6dp2=(-e1*(e2-1)*e2*e3*p1*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl6dp3=(-e1*e2*e3*e3*p1*(p2-1))/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl6de1*ee1),2)+ pow((dl6de2*ee2),2)+pow((dl6de3*ee3),2)+pow((dl6dp1*ep1),2)+pow((dl6dp2*ep2),2)+ pow((dl6dp3*ep3),2);
  }
  if (label==7){
    float dl7de1=(e2*e3*p1*p1*p2*(p3-1))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl7de2=(e1*e3*p1*p2*p2*(p3-1))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl7de3=(e1*e2*p1*p2*(p3-1)*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl7dp1=(-e1*e1*e2*e3*p2*(p3-1))/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl7dp2=(-e1*e2*e2*e3*p1*(p3-1))/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl7dp3=(-e1*e2*(e3-1)*e3*p1*p2)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl7de1*ee1),2)+ pow((dl7de2*ee2),2)+pow((dl7de3*ee3),2)+pow((dl7dp1*ep1),2)+pow((dl7dp2*ep2),2)+ pow((dl7dp3*ep3),2);
  }
  if (label==8){
    float dl8de1=(-e2*e3*p1*p1*p2*p3)/((p1-e1)*(p1-e1)*(p2-e2)*(p3-e3));
    float dl8de2=(-e1*e3*p1*p2*p2*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl8de3=(-e1*e2*p1*p2*p3*p3)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    float dl8dp1=(e1*e1*e2*e3*p2*p3)/((e1-p1)*(e1-p1)*(e2-p2)*(e3-p3));
    float dl8dp2=(e1*e2*e2*e3*p1*p3)/((e1-p1)*(e2-p2)*(e2-p2)*(e3-p3));
    float dl8dp3=(e1*e2*e3*e3*p1*p2)/((e1-p1)*(e2-p2)*(e3-p3)*(e3-p3));
    errorEvent=pow((dl8de1*ee1),2)+ pow((dl8de2*ee2),2)+pow((dl8de3*ee3),2)+pow((dl8dp1*ep1),2)+pow((dl8dp2*ep2),2)+ pow((dl8dp3*ep3),2);
  }
  //std::cout<<errorEvent<<std::endl;
  errorEvent+=(weight*weight);
  /*  std::cout<<weight*weight<<std::endl;
  std::cout<<errorEvent<<std::endl;
  std::cout<<"********"<<std::endl;
  */
  return errorEvent;
}

float AxeError(float weight, TH2F * MuonSF, TH2F * ElecSF, float pileUpWeight, int * WZCandidates, int type,TLorentzVector* analysisLepton, float TriggerEff, float triggerError=0){
  float error(0);
  int index1=WZCandidates[0];
  int index2=WZCandidates[1];
  int index3=WZCandidates[2];
  float factor1(-999), factor2(-999), factor3(-999);
  float factorError1(100000), factorError2(100000), factorError3(100000);

  if (type==0) {
    factor1= GetFactor(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    //factorError1= GetError(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factorError1=10;
    factor2= GetFactor(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    //    factorError2= GetError(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factorError2=10;
    factor3= GetFactor(ElecSF, analysisLepton[index3].Pt(),analysisLepton[index3].Eta());
    factorError3= GetError(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    
  }
  if (type==1){
    factor1= GetFactor(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factorError1= GetError(ElecSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factor2= GetFactor(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factorError2= GetError(ElecSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factor3= GetFactor(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    factorError3= GetError(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
  }
  if (type==2){
    factor1= GetFactor(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factorError1= GetError(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factor2= GetFactor(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factorError2= GetError(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factor3= GetFactor(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    factorError3= GetError(ElecSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
  }
  if (type==3){
    factor1= GetFactor(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factorError1= GetError(MuonSF, analysisLepton[index1].Pt(), analysisLepton[index1].Eta());
    factor2= GetFactor(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factorError2= GetError(MuonSF, analysisLepton[index2].Pt(), analysisLepton[index2].Eta());
    factor3= GetFactor(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
    factorError3= GetError(MuonSF, analysisLepton[index3].Pt(), analysisLepton[index3].Eta());
  }
  error= weight*weight+ pow(pileUpWeight*pileUpWeight*factorError1*factor2*factor3*TriggerEff,2)+ pow(pileUpWeight*pileUpWeight*factor1*factorError2*factor3*TriggerEff,2)+ pow(pileUpWeight*pileUpWeight*factor1*factor2*factorError3*TriggerEff,2);//+ pow(pileUpWeight*factor1*factor2*factor3*triggerError,2);
  return error;
}

double ReturnBranchingWeight(int type){
  
  //PDG values:
  
    //type: 0-(eee), 1-(eem), 2-(mme), 3-(mmm), 4-(tte), 5-(ttm), 6-(ttt), 7-(eet), 8-(mmt)
  
  if (type==-999) return 0;
  double ratioPDG[9], ratioMad[9];
  double numMad[9];
  /*
  ratioPDG[0] = 0.03363*0.1075/ (0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[1] = 0.03363*0.1057/ (0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[2] = 0.03366*0.1075/ (0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[3] = 0.03366*0.1057/ (0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[4] = 0.03370*0.1075/(0.033658*3*(0.1075+0.1057+0.1125)); 
  ratioPDG[5] = 0.03370*0.1057/(0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[6] = 0.03370*0.1125/(0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[7] = 0.03363*0.1125/(0.033658*3*(0.1075+0.1057+0.1125));
  ratioPDG[8] = 0.03366*0.1125/(0.033658*3*(0.1075+0.1057+0.1125));

  double srecko[9];
  srecko[0]=131706/1210189;
  srecko[1]=131036.0/1210189.0;
  srecko[2]=131745.0/1210189.0;
  srecko[3]=131742.0/1210189.0;
  srecko[4]=136564.0/1210189.0;
  srecko[5]=136984.0/1210189.0;
  srecko[6]=135828.0/1210189.0;
  srecko[7]=137090.0/1210189.0;
  srecko[8]=137090.0/1210189.0;

    
  numMad[0]=135772;
  numMad[1]=135591;
  numMad[2]=134010;
  numMad[3]=137417;
  numMad[4]=152809;
  numMad[5]=153569;
  numMad[6]=151913;
  numMad[7]=154520;
  numMad[8]=154395;
  double all=1309996;
  */

  ratioPDG[0]= eee_pdg;
  ratioPDG[1]= eem_pdg;
  ratioPDG[2]= mme_pdg;
  ratioPDG[3]= mmm_pdg;
  ratioPDG[4]= tte_pdg;
  ratioPDG[5]= ttm_pdg;
  ratioPDG[6]= ttt_pdg;
  ratioPDG[7]= eet_pdg;
  ratioPDG[8]= mmt_pdg;

  ratioMad[0]= eee_madgraph;
  ratioMad[1]= eem_madgraph;
  ratioMad[2]= mme_madgraph;
  ratioMad[3]= mmm_madgraph;
  ratioMad[4]= tte_madgraph;
  ratioMad[5]= ttm_madgraph;
  ratioMad[6]= ttt_madgraph;
  ratioMad[7]= eet_madgraph;
  ratioMad[8]= mmt_madgraph;
  
  

  //double BrMad=numMad[type]/all;
    

  /*   
  for (int t=0; t<9; t++){
    std::cout<<"PDG: "<<ratioPDG[t]<<std::endl;
    std::cout<<"kod: "<<numMad[t]/all<<std::endl;
//    std::cout<<"sre: "<<srecko[t]<<std::endl;
    std::cout<<"************"<<std::endl;
  }

  std::cout<<"********************"<<std::endl;
  */
  //  double weightBr= ratioPDG[type]/BrMad;
  double weightBr= ratioPDG[type]/ratioMad[type];
  
  //  if (type==6) weightBr=0.11528/0.1122376;
  //  if (type==7) std::cout<<ratioPDG[type]<<"/"<<BrMad<<"="<<weightBr<<std::endl;  
  return weightBr;
}

TLorentzVector GetMET(Float_t metModule, Float_t metPhi)
{
  Float_t px = metModule * cos(metPhi);
  Float_t py = metModule * sin(metPhi);
  
  TLorentzVector metv(px, py, 0., metModule);

  return metv;
}




#ifdef NEWMC
int determineGenType(WZTEST * cWZ){
  int genType=-999;
  
  double Wel(false), Wmu(false), Wtau(false), Zel(false), Zmu(false), Ztau(false);
  int numW(0), numZ(0), indexW(-999), indexZ1(-999), indexZ2(-999);

  for (int igl=0; igl<cWZ->genLeptons.size() ; igl++) {
    int bosonId = cWZ->genLeptons[igl].MotherBoson();
    //W
    if (abs(bosonId) == 24) {
      numW++;
      indexW=igl;
    }
    //Z
    if (abs(bosonId) == 23) {
      numZ++;
      if (numZ==2) indexZ2=igl;
      if (numZ==1) indexZ1=igl;
    }
  }
  if (numW > 1) return genType;
  if (numZ > 2) return genType;
  /*
  int Zid1=abs(cWZ->genLeptons[indexZ1].Id());
  int Zid2=abs(cWZ->genLeptons[indexZ2].Id());
  int wid=abs(cWZ->genLeptons[indexW].Id());
  */

  //hadronic decays
  if (numW==0) Wtau=true;
  if (numZ<2) Ztau=true;

  //W lepton
  if (numW>0){
    if (((cWZ->genLeptons[indexW].ComesFromTau()))) Wtau=true;
    else {
      if ((abs(cWZ->genLeptons[indexW].Id()))==11) Wel=true;
      if ((abs(cWZ->genLeptons[indexW].Id()))==13) Wmu=true;
    }
  }
  //Z lepton
  if (numZ>1){
    if ((cWZ->genLeptons[indexZ1].ComesFromTau()) && (cWZ->genLeptons[indexZ2].ComesFromTau())) Ztau=true;
    else {
      if (((abs(cWZ->genLeptons[indexZ1].Id()))==11) && ((abs(cWZ->genLeptons[indexZ2].Id()))==11)) Zel=true;
      if (((abs(cWZ->genLeptons[indexZ1].Id()))==13) && ((abs(cWZ->genLeptons[indexZ2].Id()))==13)) Zmu=true;
    }
  }
  if (Zel && Wel) genType=0;
  if (Zel && Wmu) genType=1;
  if (Zmu && Wel) genType=2;
  if (Zmu && Wmu) genType=3;
  if (Ztau && Wel) genType=4;
  if (Ztau && Wmu) genType=5;
  if (Ztau && Wtau) genType=6;
  if (Zel && Wtau) genType=7;
  if (Zmu && Wtau) genType=8;
  

  return genType;
}


#endif

TH1F* GetHistogramFromGraph(TString hname, TString gname)
{
  TString path = "../AuxiliaryFilesWZXS8TeV/";

  TFile* inputfile = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/auxiliaryFiles/gScaleSyst-8TeV.root");

  TGraphErrors* graph = (TGraphErrors*)inputfile->Get(gname);

  UInt_t nbins = graph->GetN();

  Double_t* xx = graph->GetX();
  Double_t* yy = graph->GetY();

  Float_t* range = new Float_t[nbins+1];

  for (UInt_t i=0; i<nbins; i++)
    {
      range[i+1] = (xx[i] + xx[i+1]) / 2.;
    }

  range[0]     = xx[0]       - (xx[1]       - xx[0]      ) / 2.;
  range[nbins] = xx[nbins-1] + (xx[nbins-1] - xx[nbins-2]) / 2.;

  TH1F* hist = new TH1F(hname, hname, nbins, range);
  
  hist->SetDirectory(0);

  for (UInt_t i=1; i<=nbins; i++)
    {
      hist->SetBinContent(i, fabs(yy[i-1]));
    }

  inputfile->Close();

  delete [] range;

  return hist;
}



