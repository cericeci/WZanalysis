#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "UnfoldingHistogramFactory.h"
#include "HistogramFactory.h"
// Replace this with the new tree
//#include "WZEvent.h"
#include "WZEventMCOld.h"
//#include "WZ2012Data.h"

#include "MatrixTools.h"

#include <iostream>
#include <iomanip>
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

int determineLabel(int * pass2012ICHEP, int * WZcandidates);

double weight(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt,float * eta, int label);

double promptWeight(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, 
		    int * WZcandidates, int * pass2012ICHEP, 
		    int channel, float* pt, float* eta);

float promptError(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, 
		  int * WZcandidates, int * pass2012ICHEP, 
		  int channel, float* pt, float* eta, float weight,
		  std::vector<MatrixTerm> * xterms_pf=0,
		  std::vector<MatrixTerm> * xterms_eff=0);


float GetFactor(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);

float GetError(TH2F* h2, float leptonPt, float leptonEta, float leptonPtMax= -999.);

TH2F* LoadHistogram(TString filename, TString hname, TString cname);

float MMerror(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label, float weight);

TLorentzVector GetMET(Float_t metModule, Float_t metPhi);

double deltaPhiWMET(int * WZcandidates, TLorentzVector* analysisLepton, TLorentzVector EventMET);
//float MMerror(TH2F* ElectronFR, TH2F* ElectronPR, TH2F* MuonFR, TH2F* MuonPR, int* WZcandidates, int type, float* pt, float* eta, int label);

double wTransverseMass(int index, TLorentzVector* analysisLepton, TLorentzVector EventMET);



// Put correct stat. uncertainties in electron FR map
//  (values in the input histogram are weird, way too large)
void fixElectronFRErrors(TH2F* ElectronFR) {

  std::map< pair<float,float>, float> errorMap;

  float pt(10.),eta(2.), error(0.1);
  errorMap.insert( make_pair( make_pair(12.,0.5), 0.005));
  errorMap.insert( make_pair( make_pair(12.,1.2), 0.004));
  errorMap.insert( make_pair( make_pair(12.,1.8), 0.002));
  errorMap.insert( make_pair( make_pair(12.,2.2), 0.005));

  errorMap.insert( make_pair( make_pair(17.,0.5), 0.003));
  errorMap.insert( make_pair( make_pair(17.,1.2), 0.003));
  errorMap.insert( make_pair( make_pair(17.,1.8), 0.001));
  errorMap.insert( make_pair( make_pair(17.,2.2), 0.002));

  errorMap.insert( make_pair( make_pair(22.,0.5), 0.002));
  errorMap.insert( make_pair( make_pair(22.,1.2), 0.003));
  errorMap.insert( make_pair( make_pair(22.,1.8), 0.002));
  errorMap.insert( make_pair( make_pair(22.,2.2), 0.002));

  errorMap.insert( make_pair( make_pair(27.,0.5), 0.003));
  errorMap.insert( make_pair( make_pair(27.,1.2), 0.005));
  errorMap.insert( make_pair( make_pair(27.,1.8), 0.003));
  errorMap.insert( make_pair( make_pair(27.,2.2), 0.003));

  errorMap.insert( make_pair( make_pair(33.,0.5), 0.006));
  errorMap.insert( make_pair( make_pair(33.,1.2), 0.009));
  errorMap.insert( make_pair( make_pair(33.,1.8), 0.006));
  errorMap.insert( make_pair( make_pair(33.,2.2), 0.005));

  std::map< pair<float,float>, float>::iterator it;

  for (it=errorMap.begin();it!=errorMap.end(); it++) {
    float pt  = it->first.first;
    float eta = it->first.second;
    float err = it->second;
    cout << "Error map: Pt = " << pt << "\t Eta = " << eta << "\t error = " << err << std::endl;

    ElectronFR->SetBinError( ElectronFR->FindBin(pt,eta), err);

  }


}



void printFakeRates(TH2F* ElectronFR, TH2F* MuonFR)
{
  // Print out Prompt rates
  float ptBinLimits [6] = {10,15,20,25,30,14000};
  float etaBinLimits[5] = {0., 1., 1.479, 2., 2.5};

  ofstream ele_fr_list("fakerates_ele.txt");
  ofstream mu_fr_list("fakerates_mu.txt");

  ele_fr_list << " Pt  ";
  mu_fr_list  << " Pt  ";

  for (int ieta = 0; ieta<4; ieta++) {
    ele_fr_list << "\t\t|\t" 
		<< etaBinLimits[ieta] << " - " 
		<< etaBinLimits[ieta+1];
    mu_fr_list << "\t\t|\t" 
		<< etaBinLimits[ieta] << "    -      " 
		<< etaBinLimits[ieta+1];
  }
  ele_fr_list << std::endl;
  mu_fr_list << std::endl;

  for (int ipt = 0; ipt<6; ipt++) {
    float pt = 0.5*(ptBinLimits[ipt]+ptBinLimits[ipt+1]);

    ele_fr_list << "Pt = " << pt;
    mu_fr_list  << "Pt = " << pt;
    for (int ieta = 0; ieta<4; ieta++) {
      float eta = 0.5*(etaBinLimits[ieta]+etaBinLimits[ieta+1]);
      float elepf  = GetFactor(ElectronFR, pt, eta);
      float eledpf = GetError(ElectronFR, pt, eta);
      float mupf  = GetFactor(MuonFR, pt, eta);
      float mudpf = GetError(MuonFR, pt, eta);

      ele_fr_list << "\t|\t" << setprecision(4) << elepf << " +/- " << setprecision(4) << eledpf;
      mu_fr_list  << "\t|\t" << setprecision(4) << mupf  << " +/- " << setprecision(4) << mudpf;
    }
    ele_fr_list << std::endl;
    mu_fr_list << std::endl;
  }
  
}


double computeXtermUncertainty(std::vector<MatrixTerm>  & xterms, bool useOnlyDiagonalTerms=false) {

    double sumXterms(0.);
    double sumDiagTerms(0.);

    int nXterms = 0;
    
    for (int ievt1=0; ievt1<xterms.size(); ievt1++) {
      for (int ievt2=ievt1; ievt2<xterms.size(); ievt2++) {
	for (int ilep1=0; ilep1<3; ilep1++) {
	  for (int ilep2=0; ilep2<3; ilep2++) {
	    if (xterms[ievt1].histoBin[ilep1] == xterms[ievt2].histoBin[ilep2]) {


	      // std::cout << "Found xterm " 
	      // 		<< " bins: " << xterms[ievt1].histoBin[ilep1] 
	      // 		<< "\t" <<  xterms[ievt2].histoBin[ilep2] 
	      // 		<< "\t errors = " << xterms[ievt1].effError[ilep1]
	      // 		<< "\t " << xterms[ievt2].effError[ilep2] << std::endl;
	      //	      if (ilep1 != 2 && ilep2 != 2) {
	      double x = xterms[ievt1].derivative[ilep1] * xterms[ievt2].derivative[ilep2];
	      x  *=  pow(xterms[ievt1].effError[ilep1],2);

	      if (ievt1 == ievt2) {
		if (ilep1 == ilep2) {
		  sumDiagTerms += x;
		} else {
		  //  std::cout << "Found same efficiency bin in different leptons in same event : " 
		  //    << ilep1 << "\t" << ilep2 << std::endl;
		}
	      } else {
		sumXterms += x;
		nXterms++;
	      }
	    }
	  }
	}
      }
    }

    sumXterms  *= 2;
    std::cout << "xterm : " << sqrt(sumXterms) 
	      << "\t diag.term : " << sqrt(sumDiagTerms)
	      << "\t # x-terms : " << nXterms
	      << "\t events : " << xterms.size()
	      << std::endl;


    double sum = sumDiagTerms;
    if (! useOnlyDiagonalTerms) sum+=sumXterms;

    return sum;

}



int main(int argc, char **argv)
{
  using namespace std;
  bool gotHistoBinning(false);
  char * binningFileName(0);

  char c;
  
  while ((c = getopt (argc, argv, "i:o:l:H:")) != -1)
    switch (c)
      {
	case 'H':
	  gotHistoBinning = true;
	  binningFileName = new char[strlen(optarg)+1];
	  strcpy(binningFileName,optarg);
	  break;
	default:
	  std::cout << "usage: [-k|-g|-l] [-v] [-b <binWidth>]   -i <input> -o <output> \n";
	  abort ();
	}  

  //  ofstream myfile3e, myfile3mu, myfile2e1mu, myfile1e2mu, myfileAll;
  ofstream fileNumMM;
   
  //fileNumMM.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numMM.h");
  //  fileNumMM.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numMM_met.h");
  fileNumMM.open("numMM_met.h");
  // myfile3e.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3e_Lucija.txt");
  // myfile2e1mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts2e1mu_Lucija.txt");
  // myfile1e2mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts1e2mu_Lucija.txt");
  // myfile3mu.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/WZ_evts3mu_Lucija.txt");
  //  myfileAll.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/comparisonWithJonatan/all_Lucija.txt");

  //  bool writeOutputNumbers(true);
  bool writeOutputNumbers(true);
  bool latexOutput(true);
  if (writeOutputNumbers){
    fileNumMM<<"#ifndef numMM_h"<<std::endl;
    fileNumMM<<"#define numMM_h"<<std::endl;
  }

  //  TFile * fout= new TFile("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data_driven.root", "RECREATE");  
  TFile * fout= new TFile("data_driven.root", "RECREATE");

  TH1F * hZmassMu1         = new TH1F ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
  TH1F * hZmassEl1         = new TH1F ("hZmassEl1", "hZmassEl1", 100, 60, 120);  

  const int nChannels1(5);

  TH1F * hZmassError[nChannels1];
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
  TH1F * hNjetsBigger1[nChannels1];
  TH1F * hNjetsBigger2[nChannels1];
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
  TH1F * h3Lmass1[nChannels1];
  TH1F * h3Lmass2[nChannels1];
  
  for (int myhist=0; myhist<nChannels1; myhist++){
    std::ostringstream nZmass1, nMET1, nZpt1, nLeadingJetPt1, nNjets1, nDeltaPhi1, nZlepton1pt1, nZlepton2pt1, nWleptonpt1, nMTW1, nNjetsBigger1, n3Lmass1;
    std::ostringstream nZmass2, nMET2, nZpt2, nLeadingJetPt2, nNjets2, nDeltaPhi2, nZlepton1pt2, nZlepton2pt2, nWleptonpt2, nMTW2, nNjetsBigger2, n3Lmass2;
    std::ostringstream nZmassError;
    nZmassError<<"hZmassError_"<<myhist;
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
    nNjetsBigger1<<"hNjetsBigger1_"<<myhist;
    nNjetsBigger2<<"hNjetsBigger2_"<<myhist;
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
    n3Lmass1<<"h3Lmass1_"<<myhist;
    n3Lmass2<<"h3Lmass2_"<<myhist;

    hZmass1[myhist]           = HistogramFactory::createZmassHisto(nZmass1.str().c_str(), nZmass1.str().c_str());
    hZmass2[myhist]           = HistogramFactory::createZmassHisto(nZmass2.str().c_str(), nZmass2.str().c_str());
    hZmassError[myhist]       = HistogramFactory::createZmassHisto(nZmassError.str().c_str(), nZmassError.str().c_str());
    hMET1[myhist]             = HistogramFactory::createMETHisto(nMET1.str().c_str(), nMET1.str().c_str());
    hMET2[myhist]             = HistogramFactory::createMETHisto(nMET2.str().c_str(), nMET2.str().c_str());
    hZpt_my1[myhist]          = HistogramFactory::createZptHisto(nZpt1.str().c_str(), nZpt1.str().c_str());
    hZpt_my2[myhist]          = HistogramFactory::createZptHisto(nZpt2.str().c_str(), nZpt2.str().c_str());
    hLeadingJetPt_my1[myhist] = HistogramFactory::createLeadingJetptHisto(nLeadingJetPt1.str().c_str(), nLeadingJetPt1.str().c_str());
    hLeadingJetPt_my2[myhist] = HistogramFactory::createLeadingJetptHisto(nLeadingJetPt2.str().c_str(), nLeadingJetPt2.str().c_str());
    hNjets1[myhist]           =  HistogramFactory::createNjetsHisto(nNjets1.str().c_str(), nNjets1.str().c_str());
    hNjets2[myhist]           =  HistogramFactory::createNjetsHisto(nNjets2.str().c_str(), nNjets2.str().c_str());
    hNjetsBigger1[myhist]     =  HistogramFactory::createNjetsHistoBigger(nNjetsBigger1.str().c_str(), nNjetsBigger1.str().c_str());
    hNjetsBigger2[myhist]     =  HistogramFactory::createNjetsHistoBigger(nNjetsBigger2.str().c_str(), nNjetsBigger2.str().c_str());
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
    h3Lmass1[myhist]            = HistogramFactory::create3LmassHisto(n3Lmass1.str().c_str(), n3Lmass1.str().c_str());          
    h3Lmass2[myhist]            = HistogramFactory::create3LmassHisto(n3Lmass2.str().c_str(), n3Lmass2.str().c_str());          

  }

  const int nChannels(4);
  TH1D * hZpt[nChannels];
  TH1D * hZpt_test[nChannels];
  TH1D * hLeadingJetPt[nChannels]; 
  TH1D * hLeadingJetPtError[nChannels]; 
  TH1D * hNjets[nChannels]; 
  TH1D * hNjetsError[nChannels]; 
  TH1D * hZpterror[nChannels];
  TH1D * hZptFake[nChannels];
  TH1D * hZptFake3sigmaUp[nChannels];
  TH1D * hZptFake3sigmaDown[nChannels];
  TH1D * hZptFakeError[nChannels];

  //HISTOGRAM NEEDED FOR UNFOLFING: "Zpt_", "LeadingJetPt_" 

  if (gotHistoBinning) {

    UnfoldingHistogramFactory * histoFac = UnfoldingHistogramFactory::GetInstance();
    histoFac->SetBinning(binningFileName);

  }

  for (int hist=0; hist<nChannels; hist++){
    std::ostringstream Zptname, LeadingJetname, Zpterrorname, Zpttest,LeadingJetError, ZptFakename,ZptFakename3Up, ZptFakename3Down, ZptFakeErrorname, Njetsname, Njetserror ;  
    LeadingJetname<<"LeadingJetPt_"<<(hist+1);
    LeadingJetError<<"LeadingJetError_"<<(hist+1);
    Zpttest<<"Zpt_test_"<<hist;
    Zptname<<"Zpt_"<<(hist+1);
    Zpterrorname<<"Zpterror_"<<(hist+1);
    Njetsname<<"Njets_"<<(hist+1);
    Njetserror<<"Njets_error"<<(hist+1);

    hZpt[hist]     = UnfoldingHistogramFactory::createZPtHistogram(Zptname.str().c_str(), Zptname.str().c_str());
    hZpt_test[hist]     = UnfoldingHistogramFactory::createZPtHistogram(Zpttest.str().c_str(), Zpttest.str().c_str());
    hZpterror[hist]     = UnfoldingHistogramFactory::createZPtHistogram(Zpterrorname.str().c_str(), Zpterrorname.str().c_str());
    hLeadingJetPt[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetname.str().c_str(),LeadingJetname.str().c_str());
    hLeadingJetPtError[hist] = UnfoldingHistogramFactory::createLeadingJetHistogram(LeadingJetError.str().c_str(),LeadingJetError.str().c_str()); 
    hNjets[hist] = UnfoldingHistogramFactory::createNjetsHistogram(Njetsname.str().c_str(),Njetsname.str().c_str());
    hNjetsError[hist] = UnfoldingHistogramFactory::createNjetsHistogram(Njetserror.str().c_str(),Njetserror.str().c_str());
    
    ZptFakename<<"fake_"<<(hist+1);
    ZptFakename3Up<<"fake_3sigmaUp_"<<(hist+1);
    ZptFakename3Down<<"fake_3sigmaDown"<<(hist+1);
    ZptFakeErrorname<<"error_fake_"<<(hist+1);
    /*
    hZptFake[hist] =  UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptFakename.str().c_str(), ZptFakename.str().c_str());

    hZptFake3sigmaUp[hist] =  UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptFakename3Up.str().c_str(), ZptFakename3Up.str().c_str());

    hZptFake3sigmaDown[hist] =  UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptFakename3Down.str().c_str(), ZptFakename3Down.str().c_str());

    hZptFakeError[hist] =  UnfoldingHistogramFactory::createZPtHistogram_aTGC(ZptFakeErrorname.str().c_str(), ZptFakeErrorname.str().c_str());
    */
    }
  TH2F* MuonFR;
  TH2F* MuonPR;
  TH2F* ElectronFR;
  TH2F* ElectronPR;


  //nominal!!!  
  MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet20_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR"); //nominal!!!
  //syst
  // MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet10_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR")
  // MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet15_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  // MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet30_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  // MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet35_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  // MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet40_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  //  MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet50_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");

  //MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet30_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  //MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet35_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  //MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet40_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");
  //MuonFR=LoadHistogram("auxiliaryFiles/MuFR_Moriond13_jet50_EWKcorr.root", "FR_pT_eta_EWKcorr", "MuonFR");

  MuonPR=LoadHistogram("auxiliaryFiles/MuPR_Moriond13_2012.root", "h2inverted", "MuonPR");
  
  ElectronFR=LoadHistogram("auxiliaryFiles/EleFR_Moriond13_jet35_EWKcorr.root", "fakeElH2", "ElectronFR"); //nominal!!!
  //syst:
  //  ElectronFR=LoadHistogram("auxiliaryFiles/EleFR_Moriond13_jet15_EWKcorr.root", "fakeElH2", "ElectronFR");
  //ElectronFR=LoadHistogram("auxiliaryFiles/EleFR_Moriond13_jet50_EWKcorr.root", "fakeElH2", "ElectronFR");

  ElectronPR=LoadHistogram("auxiliaryFiles/ElePR_Moriond13_2012.root", "h2inverted", "ElectronPR");

  fixElectronFRErrors(ElectronFR);
  printFakeRates(ElectronFR, MuonFR);
    
  const int leptonNumber(4);
  const float electronMass(0.000511);
  const float muonMass(0.106);
  int numZ(0), numW(0), numMET(0), num3e(0), num2e1mu(0), num1e2mu(0), num3mu(0), numMET3e(0), numMET2e1mu(0), numMET1e2mu(0), numMET3mu(0);
  double N_good_3e(0), N_good_2e1mu(0), N_good_1e2mu(0), N_good_3mu(0);
  double N_fake_3e(0), N_fake_2e1mu(0), N_fake_1e2mu(0), N_fake_3mu(0);
  double ferror_3e(0),ferror_2e1mu(0), ferror_1e2mu(0), ferror_3mu(0);
  double error_3e(0),error_2e1mu(0), error_1e2mu(0), error_3mu(0);
  //counters for matrix method
  double N_matrix[8][4]={0,0,0,0,0,0,0,0,
			 0,0,0,0,0,0,0,0,
			 0,0,0,0,0,0,0,0,
			 0,0,0,0,0,0,0,0,
                     };
  int passLoose[4]={0,0,0,0};
  double N_good[4]={0,0,0,0};
  double error[4]={0,0,0,0};
  double N_fake[4]={0,0,0,0};
  double ferror[4]={0,0,0,0};
  double N_fake_TEST[4]={0,0,0,0};

  double N_good_minus[4]={0,0,0,0};
  double error_minus[4]={0,0,0,0};
  double N_fake_minus[4]={0,0,0,0};
  double ferror_minus[4]={0,0,0,0};


  double N_good_plus[4]={0,0,0,0};
  double error_plus[4]={0,0,0,0};
  double N_fake_plus[4]={0,0,0,0};
  double ferror_plus[4]={0,0,0,0};

  // Vuko's additions
  double N_prompt[4]={0,0,0,0};
  double dN_prompt[4]={0,0,0,0};
  double N_fakelep[4]={0,0,0,0};
  double dN_fakelep[4]={0,0,0,0};



  std::vector<MatrixTerm> xtermsFR[4];
  std::vector<MatrixTerm> xtermsPR[4];


  double data_driven_syst[4]={0.029, 0.022, 0.032, 0.028};
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


  WZEventMCOld *cWZ= new WZEventMCOld(wz_tTree);

  //  WZEvent *cWZ= new WZEvent(wz_tTree);
  //  WZ2012Data *cWZ= new WZ2012Data(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;

  //////////////////////////////////////////////
  //   EVENT LOOP


  for  (Int_t k = 0; k<events  /* && k<1000000 */ ;k++) {
    wz_tTree->GetEntry(k);
    cWZ->ReadEvent();

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
    TLorentzVector EventMET;
    TLorentzVector analysisLepton[leptonNumber];
    
    double pfmet=cWZ->pfmetTypeI;
    double pfmetphi=cWZ->pfmetTypeIphi;
    EventMET = GetMET(pfmet, pfmetphi);

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
      if ((analysisLepton[i].Pt()>10) && (analysisLepton[i].Pt()!=-9999))
	lepNum++;
      if ((bdt[i]<100) && (analysisLepton[i].Pt()>10)){
	good_electrons.push_back(i);
	v_nizEl[i].SetPtEtaPhiM(analysisLepton[i].Pt(),eta[i], phi[i], electronMass);
	v_3Lepton=v_3Lepton+v_nizEl[i];
      }
      if ((bdt[i]>100) && (analysisLepton[i].Pt()>10)){
	good_muons.push_back(i);
	v_nizMu[i].SetPtEtaPhiM(analysisLepton[i].Pt(),eta[i], phi[i], muonMass);
	v_3Lepton=v_3Lepton+v_nizMu[i];
      }
    }
    if (lepNum!=3) continue; ///???

    if (v_3Lepton.M()<100) continue;
    
    justCount++;
    bool foundZel(false), foundZmu(false);
    double massMu(-999), massEl(0), Zpt(0);

    foundZmu=  Z_muons(cWZ, &good_muons, WZcandidates, v_nizMu, analysisLepton, ch, massMu, Zpt);
    if (foundZmu){
      hZmassMu1->Fill(massMu);
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

    numZ++;

    double Zmass;
    if (foundZel) Zmass=massEl;
    else Zmass=massMu;

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
    if (!passDeltaRWleptZlept(WZcandidates, analysisLepton)) continue;
    numW++;  
    
    int type;
    //0-EEE, 1-EEM, 2-EMM, 3-MMM

    bool ev[5]={false, false, false, false, true};
    double MMweights[4]={.0,.0,.0,.0};
    double MMerrors[4]={.0,.0,.0,.0};
    double MMweights_fake[5]={.0,.0,.0,.0};
    double MMerrors_fake[4]={.0,.0,.0,.0};
    
    // Vuko's additions
    double promptWeights[4]={.0,.0,.0,.0};
    double promptMMErrors2[4]={.0,.0,.0,.0};


    int label= determineLabel(pass2012ICHEP, WZcandidates);
 

    float Nttt(0);
    if (label==1) Nttt=1;     
    
    if (ev3e){
      num3e++;
      type=0;
      ev[0]=true;
    }
    if (ev2e1mu){
      num2e1mu++;
      type=1;
      ev[1]=true;
    }
    if (ev1e2mu){
      num1e2mu++;
      type=2;
      ev[2]=true;
    }

    if (ev3mu){
      num3mu++;
      type=3;
      ev[3]=true;
    }
    ////////////////////MATRIX METHOD PART//////////////////////////
    double MMweights_fake_h(0);
    double MMerrors_fake_h(0);
    
    for (int matrix=0; matrix<4; matrix++){
      if (!ev[matrix]) continue;
      MMweights[matrix]=weight(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, matrix, pt, eta, label);
      promptWeights[matrix] = 
	promptWeight(ElectronFR, ElectronPR, MuonFR, MuonPR, 
		     WZcandidates, pass2012ICHEP, 
		     matrix, pt, eta);
      promptMMErrors2[matrix] = promptError(ElectronFR, ElectronPR, MuonFR, MuonPR, 
					    WZcandidates, pass2012ICHEP, 
					    matrix, pt, eta, promptWeights[matrix], 
					    &xtermsFR[matrix],
					    &xtermsPR[matrix]); // xxx

      //      double xc_diff = MMweights[matrix] - promptWeights[matrix];
      //      if (fabs(xc_diff) > 0.0001) {
	  // std::cout << "Channel: " << matrix
	  // 	    << "\t cat : " << label
	  // 	    << "\t X-check weights : " << MMweights[matrix] << "\t" << promptWeights[matrix] 
	  // 	    << "\t diff = " << MMweights[matrix] - promptWeights[matrix] 
	  // 	    << std::endl;
	  //      }
      MMweights_fake[matrix]=(Nttt-MMweights[matrix]);
      MMweights_fake_h=(Nttt-MMweights[matrix]);
      MMerrors[matrix]=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, matrix, pt, eta, label, MMweights[matrix]);
      MMerrors_fake[matrix]=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, matrix, pt, eta, label, (Nttt-MMweights[matrix]));
      MMerrors_fake_h=MMerror(ElectronFR, ElectronPR, MuonFR, MuonPR, WZcandidates, matrix, pt, eta, label, (Nttt-MMweights[matrix]));

      // if (fabs(MMerrors[matrix] -  promptMMErrors2[matrix] ) > 0.00001) 
      // 	std::cout << "Channel: " << matrix
      // 		  << "\t cat : " << label
      // 		  << "\t X-check error : " << MMerrors[matrix] << "\t" << promptMMErrors2[matrix]
      // 		  << "\t diff = " << MMerrors[matrix] -  promptMMErrors2[matrix]
      // 		  << std::endl;


    }

    double deltaPhi=deltaPhiWMET(WZcandidates,analysisLepton, EventMET);
    int indW=WZcandidates[2];
    int indZ1=WZcandidates[0];
    int indZ2=WZcandidates[1];
    double WleptPt= pt[indW];  
    double Zlept1Pt= pt[indZ1];
    double Zlept2Pt= pt[indZ2];
    double MTW= wTransverseMass(indW, analysisLepton,EventMET);

    for (int fill2=0; fill2<5; fill2++){
      if (!ev[fill2]) continue;
      hZmass1[fill2]->Fill(Zmass, MMweights_fake_h);
      hMET1[fill2]->Fill(met, MMweights_fake_h);
      hZpt_my1[fill2]->Fill(Zpt, MMweights_fake_h);
      hLeadingJetPt_my1[fill2]->Fill(leadingRecoJetPt, MMweights_fake_h);
      hNjets1[fill2]->Fill(nRecoJets, MMweights_fake_h);
      hNjetsBigger1[fill2]->Fill(nRecoJets, MMweights_fake_h);
      hDeltaPhi1[fill2]->Fill(deltaPhi, MMweights_fake_h);
      hZlepton1pt1[fill2]->Fill(Zlept1Pt,MMweights_fake_h);
      hZlepton2pt1[fill2]->Fill(Zlept2Pt,MMweights_fake_h);
      hWleptonpt1[fill2]->Fill(WleptPt,MMweights_fake_h);
      hMTW1[fill2]->Fill(MTW, MMweights_fake_h);
      //      h3Lmass1[fill2]->Fill(v_3Lepton.M(), MMweights_fake_h);
    }    
    //////////////////////////////////////MET CUT//////////////////////////

    
    if ((cWZ->pfmetTypeI)<30) continue;  
    //    if ((cWZ->pfmet)<30) continue;  

    //    passLoose++;
    for (int mm=0; mm<4; mm++){
      if (!ev[mm]) continue;
      N_matrix[label-1][mm]++;
      passLoose[mm]++;
    }
    for (int fill3=0; fill3<5; fill3++){
      if (!ev[fill3]) continue;

      hZmass2[fill3]->Fill(Zmass, MMweights_fake_h);
      //      hZmassError[fill3]->Fill(Zmass, MMerrors[filH]);
      hZmassError[fill3]->Fill(Zmass, MMerrors_fake_h);
      hMET2[fill3]->Fill(met, MMweights_fake_h);
      hZpt_my2[fill3]->Fill(Zpt, MMweights_fake_h);
      hLeadingJetPt_my2[fill3]->Fill(leadingRecoJetPt, MMweights_fake_h);
      hNjets2[fill3]->Fill(nRecoJets, MMweights_fake_h);
      hNjetsBigger2[fill3]->Fill(nRecoJets, MMweights_fake_h);
      hDeltaPhi2[fill3]->Fill(deltaPhi, MMweights_fake_h);
      hZlepton1pt2[fill3]->Fill(Zlept1Pt,MMweights_fake_h);
      hZlepton2pt2[fill3]->Fill(Zlept2Pt,MMweights_fake_h);
      hWleptonpt2[fill3]->Fill(WleptPt,MMweights_fake_h);
      hMTW2[fill3]->Fill(MTW, MMweights_fake_h);
      //      h3Lmass2[fill3]->Fill(v_3Lepton.M(), MMweights_fake_h);
    }    
    ////////////////////////////////////////////////////////////////////////
    
    int Wlepton_charge= ch[WZcandidates[2]];
    

    for (int final=0; final<4; final++){
      if (!ev[final]) continue;
      N_good[final]+=MMweights[final];
      error[final]+=MMerrors[final];
      N_fake[final]+=(Nttt-MMweights[final]);
      ferror[final]+=MMerrors_fake[final];
      if (Wlepton_charge==-1){
	N_good_minus[final]+=MMweights[final];
	error_minus[final]+=MMerrors[final];
	N_fake_minus[final]+=(Nttt-MMweights[final]);
	ferror_minus[final]+=MMerrors_fake[final];
      }
      if (Wlepton_charge==1){
	N_good_plus[final]+=MMweights[final];
	error_plus[final]+=MMerrors[final];
	N_fake_plus[final]+=(Nttt-MMweights[final]);
	ferror_plus[final]+=MMerrors_fake[final];
      }
      double  w_prompt = promptWeights[final];
      //      N_prompt[final]  += promptWeights[final];
      N_prompt[final]  += w_prompt;
      dN_prompt[final] += promptMMErrors2[final];
      N_fakelep[final] += Nttt - w_prompt;
      dN_fakelep[final] += pow(Nttt - w_prompt,2);

    }


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

    //std::cout<<leadingRecoJetPt<<std::endl;
    
    for (int filH=0; filH<nChannels; filH++){
      if (ev[filH]){
	hZpt[filH]->Fill(Zpt, MMweights[filH]);
	hZpt_test[filH]->Fill(Zpt, MMweights[filH]);
	hZpterror[filH]->Fill(Zpt, MMerrors[filH]);
	hLeadingJetPt[filH]->Fill(leadingRecoJetPt, MMweights[filH]);
	hLeadingJetPtError[filH]->Fill(leadingRecoJetPt, MMerrors[filH]);
	hNjets[filH]->Fill(nRecoJets, MMweights[filH]);
	hNjetsError[filH]->Fill(nRecoJets, MMerrors[filH]);
	//	hZptFake[filH]->Fill(Zpt, MMweights_fake[filH]);
	//hZptFakeError[filH]->Fill(Zpt, MMerrors_fake[filH]);
      }
    }
    //filling Zpt and errors
    
  }
    //filling histogram errors
  
  //  std::cout<<"test"<<hZmassMu1->GetNbinsX()<<std::endl;
  double check[4]={0,0,0,0};
  
  for (int er=0; er< nChannels; er++){
    //double binError(0);
    for (int nzpt=0; nzpt<(hZpterror[er]->GetNbinsX()+1); nzpt++){
      double binError=hZpterror[er]->GetBinContent(nzpt);
      // check[er]+=binError;
      hZpt[er]->SetBinError(nzpt, sqrt(binError));
    }
    for (int ljet=0; ljet<(hLeadingJetPtError[er]->GetNbinsX()+1); ljet++){
      double binError2=hLeadingJetPtError[er]->GetBinContent(ljet);
      hLeadingJetPt[er]->SetBinError(ljet, sqrt(binError2));
    }
    for (int cnjet=0; cnjet<(hNjetsError[er]->GetNbinsX()+1); cnjet++){
      double binError2=hNjetsError[er]->GetBinContent(cnjet);
      hNjets[er]->SetBinError(cnjet, sqrt(binError2));
    }
    /*
    for (int czmass=0; czmass<(hZmassError[er]->GetNbinsX()+1); czmass++){
      double binError2=hZmassError[er]->GetBinContent(czmass);
      hZmass2[er]->SetBinError(czmass, sqrt(binError2));
    }
    */
    //for Senka TGC
    /*
    for (int zpt2=0; zpt2<(hZptFakeError[er]->GetNbinsX()+1); zpt2++){
      double binErrorStat2=hZptFakeError[er]->GetBinContent(zpt2);
      check[er]+=binErrorStat2;
      double binContent=hZptFake[er]->GetBinContent(zpt2);
      double binErrorSyst2= pow(data_driven_syst[er]*binContent,2);
      double binError2=binErrorStat2+binErrorSyst2;
      hZptFake[er]->SetBinError(zpt2, sqrt(binError2));
      //      std::cout<<"BIN CONTENT: "<<binContent<<std::endl;
      //std::cout<<"BIN ERROR: "<<sqrt(binError2)<<std::endl;
      //std::cout<<"***************"<<std::endl;;
      hZptFake3sigmaUp[er]-> SetBinContent(zpt2, (binContent+3*sqrt(binError2)));
      //      hZptFake3sigmaUp[er]->SetBinError(zpt2, sqrt(binError2));
      if ((binContent-3*sqrt(binError2))>0)
	hZptFake3sigmaDown[er]-> SetBinContent(zpt2, (binContent-3*sqrt(binError2)));
      else
	hZptFake3sigmaDown[er]-> SetBinContent(zpt2, 0);
      //      hZptFake3sigmaDown[er]->SetBinError(zpt2, sqrt(binError2));
    }
    */
  }
  
  std::cout<<"CHECK:"<<std::endl;
  std::cout<<sqrt(check[0])<<std::endl;
  std::cout<<sqrt(check[1])<<std::endl;
  std::cout<<sqrt(check[2])<<std::endl;
  std::cout<<sqrt(check[3])<<std::endl;

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
  std::cout<<"N_good_3e:     "<<N_good[0]<<std::endl;
  std::cout<<"error_3e:     "<<sqrt(error[0])<<std::endl;
  std::cout<<"N_good_2e1mu:  "<<N_good[1]<<std::endl;
  std::cout<<"error_2e1mu:  "<<sqrt(error[1])<<std::endl;
  std::cout<<"N_good_1e2mu:  "<<N_good[2]<<std::endl;
  std::cout<<"error_1e2mu:  "<<sqrt(error[2])<<std::endl;
  std::cout<<"N_good_3mu:    "<<N_good[3]<<std::endl;
  std::cout<<"error_3mu:    "<<sqrt(error[3])<<std::endl;
  std::cout<<"*****************************************************"<<std::endl;
  std::cout<<"N_fake_3e:     "<<N_fake[0]<<std::endl;
  std::cout<<"N_fake_3e:     "<<N_fake_TEST[0]<<std::endl;

  std::cout<<"ferror_3e:     "<<sqrt(ferror[0])<<std::endl;
  std::cout<<"N_fake_2e1mu:  "<<N_fake[1]<<std::endl;
  std::cout<<"N_fake_2e1mu:  "<<N_fake_TEST[1]<<std::endl;
  std::cout<<"ferror_2e1mu:  "<<sqrt(ferror[1])<<std::endl;
  std::cout<<"N_fake_1e2mu:  "<<N_fake[2]<<std::endl;
  std::cout<<"N_fake_1e2mu:  "<<N_fake_TEST[2]<<std::endl;
  std::cout<<"ferror_1e2mu:  "<<sqrt(ferror[2])<<std::endl;
  std::cout<<"N_fake_3mu:    "<<N_fake[3]<<std::endl;
  std::cout<<"N_fake_3mu:    "<<N_fake_TEST[3]<<std::endl;
  std::cout<<"ferror_3mu:    "<<sqrt(ferror[3])<<std::endl;  


  std::cout<<"/////////////////MATRIX METHOD RESULT !!!!!!!W-!!!!!!!//////////////////"<<std::endl;
  std::cout<<"N_good_3e:     "<<N_good_minus[0]<<std::endl;
  std::cout<<"error_3e:     "<<sqrt(error_minus[0])<<std::endl;
  std::cout<<"N_good_2e1mu:  "<<N_good_minus[1]<<std::endl;
  std::cout<<"error_2e1mu:  "<<sqrt(error_minus[1])<<std::endl;
  std::cout<<"N_good_1e2mu:  "<<N_good_minus[2]<<std::endl;
  std::cout<<"error_1e2mu:  "<<sqrt(error_minus[2])<<std::endl;
  std::cout<<"N_good_3mu:    "<<N_good_minus[3]<<std::endl;
  std::cout<<"error_3mu:    "<<sqrt(error_minus[3])<<std::endl;
  std::cout<<"*****************************************************"<<std::endl;
  std::cout<<"N_fake_3e:     "<<N_fake_minus[0]<<std::endl;
  std::cout<<"ferror_3e:     "<<sqrt(ferror_minus[0])<<std::endl;
  std::cout<<"N_fake_2e1mu:  "<<N_fake_minus[1]<<std::endl;
  std::cout<<"ferror_2e1mu:  "<<sqrt(ferror_minus[1])<<std::endl;
  std::cout<<"N_fake_1e2mu:  "<<N_fake_minus[2]<<std::endl;
  std::cout<<"ferror_1e2mu:  "<<sqrt(ferror_minus[2])<<std::endl;
  std::cout<<"N_fake_3mu:    "<<N_fake_minus[3]<<std::endl;
  std::cout<<"ferror_3mu:    "<<sqrt(ferror_minus[3])<<std::endl;  


  std::cout<<"/////////////////MATRIX METHOD RESULT !!!!!!!W+!!!!!!!//////////////////"<<std::endl;
  std::cout<<"N_good_3e:     "<<N_good_plus[0]<<std::endl;
  std::cout<<"error_3e:     "<<sqrt(error_plus[0])<<std::endl;
  std::cout<<"N_good_2e1mu:  "<<N_good_plus[1]<<std::endl;
  std::cout<<"error_2e1mu:  "<<sqrt(error_plus[1])<<std::endl;
  std::cout<<"N_good_1e2mu:  "<<N_good_plus[2]<<std::endl;
  std::cout<<"error_1e2mu:  "<<sqrt(error_plus[2])<<std::endl;
  std::cout<<"N_good_3mu:    "<<N_good_plus[3]<<std::endl;
  std::cout<<"error_3mu:    "<<sqrt(error_plus[3])<<std::endl;
  std::cout<<"*****************************************************"<<std::endl;
  std::cout<<"N_fake_3e:     "<<N_fake_plus[0]<<std::endl;
  std::cout<<"ferror_3e:     "<<sqrt(ferror_plus[0])<<std::endl;
  std::cout<<"N_fake_2e1mu:  "<<N_fake_plus[1]<<std::endl;
  std::cout<<"ferror_2e1mu:  "<<sqrt(ferror_plus[1])<<std::endl;
  std::cout<<"N_fake_1e2mu:  "<<N_fake_plus[2]<<std::endl;
  std::cout<<"ferror_1e2mu:  "<<sqrt(ferror_plus[2])<<std::endl;
  std::cout<<"N_fake_3mu:    "<<N_fake_plus[3]<<std::endl;
  std::cout<<"ferror_3mu:    "<<sqrt(ferror_plus[3])<<std::endl;  




  TString namesMM[8]={"N_TTT","N_TTF", "N_TFT", "N_TFF", "N_FTT", "N_FTF", "N_FFT", "N_FFF"};
  std::cout<<"Additional numbers for matrix method:"<<std::endl;
  for (int channel=0; channel<4; channel++){
    std::cout<<"***channel: "<<channel<<std::endl;
    for (int addNum=0; addNum<8; addNum++){
      std::cout<<namesMM[addNum]<<" = "<<N_matrix[addNum][channel]<<std::endl;
    } 
    std::cout<<"Loose sample: "<<passLoose[channel]<<std::endl;
  }

  if (latexOutput)

   
    {
      std::cout<<"channel & Ngood & Nfake \\\\ "<<std::endl;
      std::cout<<"\\hline"<<std::endl;
      std::cout<<"3e     & $"<<N_good[0]    <<"\\pm"<<sqrt(error[0])    <<"$ & $"<<N_fake[0]    <<"\\pm"<<sqrt(ferror[0])    <<"$ \\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
      std::cout<<"2e1mu  & $"<<N_good[1] <<"\\pm"<<sqrt(error[1]) <<"$ & $"<<N_fake[1] <<"\\pm"<<sqrt(ferror[1]) <<"$ \\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
      std::cout<<"1e2mu  & $"<<N_good[2] <<"\\pm"<<sqrt(error[2]) <<"$ & $"<<N_fake[2] <<"\\pm"<<sqrt(ferror[2]) <<"$ \\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
      std::cout<<"3mu    & $"<<N_good[3]   <<"\\pm"<<sqrt(error[3])   <<"$ & $"<<N_fake[3]   <<"\\pm"<<sqrt(ferror[3])   <<"$ \\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
  //write in file:
  if (writeOutputNumbers){
    fileNumMM<<"#define dN_fake3e "<<N_fake[0]<<std::endl;
    fileNumMM<<"#define dN_fake2e1mu "<<N_fake[1]<<std::endl;
    fileNumMM<<"#define dN_fake1e2mu "<<N_fake[2]<<std::endl;
    fileNumMM<<"#define dN_fake3mu "<<N_fake[3]<<std::endl;
    fileNumMM<<"#define dsN_fake3e "<<sqrt(ferror[0])<<std::endl;
    fileNumMM<<"#define dsN_fake2e1mu "<<sqrt(ferror[1])<<std::endl;
    fileNumMM<<"#define dsN_fake1e2mu "<<sqrt(ferror[2])<<std::endl;
    fileNumMM<<"#define dsN_fake3mu "<<sqrt(ferror[3])<<std::endl;

    fileNumMM<<"#define dN_good3e "<<N_good[0]<<std::endl;
    fileNumMM<<"#define dN_good2e1mu "<<N_good[1]<<std::endl;
    fileNumMM<<"#define dN_good1e2mu "<<N_good[2]<<std::endl;
    fileNumMM<<"#define dN_good3mu "<<N_good[3]<<std::endl;
    fileNumMM<<"#define dsN_good3e "<<sqrt(error[0])<<std::endl;
    fileNumMM<<"#define dsN_good2e1mu "<<sqrt(error[1])<<std::endl;
    fileNumMM<<"#define dsN_good1e2mu "<<sqrt(error[2])<<std::endl;
    fileNumMM<<"#define dsN_good3mu "<<sqrt(error[3])<<std::endl;

    fileNumMM<<"#define dN_fake3e_minus "<<N_fake_minus[0]<<std::endl;
    fileNumMM<<"#define dN_fake2e1mu_minus "<<N_fake_minus[1]<<std::endl;
    fileNumMM<<"#define dN_fake1e2mu_minus "<<N_fake_minus[2]<<std::endl;
    fileNumMM<<"#define dN_fake3mu_minus "<<N_fake_minus[3]<<std::endl;
    fileNumMM<<"#define dsN_fake3e_minus "<<sqrt(ferror_minus[0])<<std::endl;
    fileNumMM<<"#define dsN_fake2e1mu_minus "<<sqrt(ferror_minus[1])<<std::endl;
    fileNumMM<<"#define dsN_fake1e2mu_minus "<<sqrt(ferror_minus[2])<<std::endl;
    fileNumMM<<"#define dsN_fake3mu_minus "<<sqrt(ferror_minus[3])<<std::endl;

    fileNumMM<<"#define dN_good3e_minus "<<N_good_minus[0]<<std::endl;
    fileNumMM<<"#define dN_good2e1mu_minus "<<N_good_minus[1]<<std::endl;
    fileNumMM<<"#define dN_good1e2mu_minus "<<N_good_minus[2]<<std::endl;
    fileNumMM<<"#define dN_good3mu_minus "<<N_good_minus[3]<<std::endl;
    fileNumMM<<"#define dsN_good3e_minus "<<sqrt(error_minus[0])<<std::endl;
    fileNumMM<<"#define dsN_good2e1mu_minus "<<sqrt(error_minus[1])<<std::endl;
    fileNumMM<<"#define dsN_good1e2mu_minus "<<sqrt(error_minus[2])<<std::endl;
    fileNumMM<<"#define dsN_good3mu_minus "<<sqrt(error_minus[3])<<std::endl;

    fileNumMM<<"#define dN_fake3e_plus "<<N_fake_plus[0]<<std::endl;
    fileNumMM<<"#define dN_fake2e1mu_plus "<<N_fake_plus[1]<<std::endl;
    fileNumMM<<"#define dN_fake1e2mu_plus "<<N_fake_plus[2]<<std::endl;
    fileNumMM<<"#define dN_fake3mu_plus "<<N_fake_plus[3]<<std::endl;
    fileNumMM<<"#define dsN_fake3e_plus "<<sqrt(ferror_plus[0])<<std::endl;
    fileNumMM<<"#define dsN_fake2e1mu_plus "<<sqrt(ferror_plus[1])<<std::endl;
    fileNumMM<<"#define dsN_fake1e2mu_plus "<<sqrt(ferror_plus[2])<<std::endl;
    fileNumMM<<"#define dsN_fake3mu_plus "<<sqrt(ferror_plus[3])<<std::endl;

    fileNumMM<<"#define dN_good3e_plus "<<N_good_plus[0]<<std::endl;
    fileNumMM<<"#define dN_good2e1mu_plus "<<N_good_plus[1]<<std::endl;
    fileNumMM<<"#define dN_good1e2mu_plus "<<N_good_plus[2]<<std::endl;
    fileNumMM<<"#define dN_good3mu_plus "<<N_good_plus[3]<<std::endl;
    fileNumMM<<"#define dsN_good3e_plus "<<sqrt(error_plus[0])<<std::endl;
    fileNumMM<<"#define dsN_good2e1mu_plus "<<sqrt(error_plus[1])<<std::endl;
    fileNumMM<<"#define dsN_good1e2mu_plus "<<sqrt(error_plus[2])<<std::endl;
    fileNumMM<<"#define dsN_good3mu_plus "<<sqrt(error_plus[3])<<std::endl;


    fileNumMM.close();
  }


  fout->cd();
  fout->Write();
  fout->Close();

  ofstream texMMoutput("mmresults.tex");


  texMMoutput << "\\begin{tabular}{ccccccc} \\hline \n";

  texMMoutput << "Channel "
	      << "  &  $N_{prompt}$" 
	      << "  &  error "
	      << "   &  $\\sqrt{A}$"
	      << "   &  $\\sqrt{B}$"
	      << "   &  $\\sqrt{C}$"
	      << "   &  $\\sqrt{B_{diag}}$"
	      << "   &  $\\sqrt{C_{diag}}$"
	      << "\t  \\\\ \\hline" << std::endl;


  std::string names[4] = {"$eee$", "$ee\\mu$", "$\\mu\\mu e$", "$\\mu\\mu\\mu$"};

  for (int chan=0; chan<4; chan++) {
    std::cout << "Prompt weight = " << N_prompt[chan] 
	      << " +/- " << sqrt(dN_prompt[chan]) 
	      << "\t Fake weight = " << N_fakelep[chan] 
	      << " +/- " << sqrt(dN_fakelep[chan]) 
	      << std::endl;

    double sumXtermsFR = computeXtermUncertainty(xtermsFR[chan]);
    double sumXtermsFR_diag = computeXtermUncertainty(xtermsFR[chan],true);
    double sumXtermsPR      = computeXtermUncertainty(xtermsPR[chan]);
    double sumXtermsPR_diag = computeXtermUncertainty(xtermsPR[chan],true);




    double total_prompt_err = sqrt(dN_prompt[chan] + sumXtermsFR + sumXtermsFR );
    std::cout << "sum of weights (A) : " << sqrt(dN_prompt[chan])
	      << "\t xterm FR (B) : " << sqrt(sumXtermsFR) 
	      << "\t xterm FR (B) diag. : " << sqrt(sumXtermsFR_diag) 
	      << "\t xterm PR (C): " << sqrt(sumXtermsPR) 
	      << "\t xterm PR (C) diag: " << sqrt(sumXtermsPR_diag) 
	      << "\t total err = " << total_prompt_err
	      << "\t events : " << xtermsFR[chan].size()
	      << std::endl;

    texMMoutput << names[chan]
		<< "\t & \t" << N_prompt[chan] 
		<< "\t & \t" << total_prompt_err 
		<< "\t & \t" << sqrt(dN_prompt[chan])
		<< "\t & \t" << setprecision(5) << sqrt(sumXtermsFR) 
		<< "\t & \t" << setprecision(5) << sqrt(sumXtermsPR) 
		<< "\t & \t" << setprecision(5) << sqrt(sumXtermsFR_diag) 	
		<< "\t & \t" << setprecision(5) << sqrt(sumXtermsPR_diag) 	
	<< "\t \\\\" << std::endl;


  }

  texMMoutput << "\\end{tabular} \\hline \n";
  texMMoutput.close();

}



