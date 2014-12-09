#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include "constants.h"
#include "systematics.h"

#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "SystematicsManager.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#define DEBUG_LEVEL 0
// #define LUMINOSITY 19602

TH1D* GetTauCorrection(TH1D* hftau, std::string name="bla") {
  // Input: ftau histogram
  // Output: 1-ftau histogram

  TH1D* h = (TH1D*) hftau->Clone(name.c_str());
  h->Reset();

  for (int i=0; i<=hftau->GetNbinsX()+1; i++) {
    h->SetBinContent(i,1.- hftau->GetBinContent(i));
    h->SetBinError(i,hftau->GetBinError(i));
  }

  return h;
}



using namespace std;

TH1D* Unfold(string unfAlg, RooUnfoldResponse* response, 
	     RooUnfold::ErrorTreatment errorType,
	     TH1D* hData, TH1D* hSumBG, int Kterm, 
	     string hOutName, bool useOverFlow,
	     TH1D * h_mcTruth=0)
{

  if (hData==0 || hSumBG == 0 || response == 0) {
    std::cout << "Entering with missing data: cannot unfold \n";
    return 0;
  }

  if (useOverFlow) response->UseOverflow();
  RooUnfold* RObject;
  TH1D * hDataClone = (TH1D*) hData->Clone();
  //   hDataClone->Add(hSumBG, -1);

  //  RooUnfold::ErrorTreatment doerror = RooUnfold::kCovariance;
  //    RooUnfold::ErrorTreatment doerror = RooUnfold::kErrors;
  // RooUnfold::ErrorTreatment doerror = RooUnfold::kCovToy;

  if (unfAlg == "SVD") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kSVD,   response, hDataClone, Kterm);
  } else   if (unfAlg == "Bayes") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kBayes, response, hDataClone, Kterm);
  } else if (unfAlg=="Invert"){
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kInvert, response, hDataClone, Kterm);
  }
  else {
    std::cout << "Unknown unfolding algorithm " << unfAlg << std::endl;
    return 0;
  }
  RObject->SetVerbose(0);
  TH1D* hCorrected = (TH1D*) RObject->Hreco(errorType);
  std::cout << "Unfolded data histogram name: " << hOutName << std::endl;

  hCorrected->SetName(hOutName.c_str());

  RObject->PrintTable(cout, h_mcTruth, (RooUnfold::ErrorTreatment)errorType);
  
  return hCorrected;
}
  
// Input
// -v <variable_Name>
// -a <Unfolding_Algo_to_use>
// -d <dataFile>
// -r <responseMatrixFile>
// [-b <background_file>|-B <background_list> ]
// -t tau fraction file

int main(int argc, char **argv) {


  char * responseFileName(0);
  char * dataFileName(0);
  char * unfoldingAlgo(0);
  char * variableName(0);

  char * bgFileName(0);
  char * bgListName(0);
  char * tauFractionFileName(0);

  bool gotDataFile = false;
  bool gotResponse = false;
  bool gotUnfoldingAlgo = false;
  bool gotVarName  = false;
  bool gotBackgroundFile = false;
  bool gotBackgroundList = false;
  bool gotTauFraction = false;  
  bool gotKterm    = false;
  int  inputKterm = -1;
  bool gotErrorTreatment = false;
  char errorName;
  bool gotSystematicsConfig = false;
  char * systConfigFile(0);

  char c;

  while ((c = getopt (argc, argv, "r:d:a:v:b:B:t:k:E:S:")) != -1)
    switch (c)
      {
      case 'r':
	gotResponse = true;
	responseFileName = new char[strlen(optarg)+1];
	strcpy(responseFileName,optarg);
	break;
      case 'd':
	gotDataFile = true;
	dataFileName = new char[strlen(optarg)+1];
	strcpy(dataFileName,optarg);
	break;
	// TODOOOO: define the following strings!!!!!!!
      case 'b':
	gotBackgroundFile = true;
	bgFileName = new char[strlen(optarg)+1];
	strcpy(bgFileName,optarg);
	break;
      case 'B':
	gotBackgroundList = true;
	bgListName = new char[strlen(optarg)+1];
	strcpy(bgListName,optarg);
	break;
      case 't':
	gotTauFraction = true;
	tauFractionFileName = new char[strlen(optarg)+1];
	strcpy(tauFractionFileName,optarg);
	break;
      case 'v':
	gotVarName = true;
	variableName = new char[strlen(optarg)+1];
	strcpy(variableName,optarg);
	break;
      case 'a':
	gotUnfoldingAlgo = true;
	unfoldingAlgo = new char[strlen(optarg)+1];
	strcpy(unfoldingAlgo,optarg);
	break;
      case 'k':
	gotKterm = true;
	inputKterm = atoi(optarg);
	break;
      case 'E':
	gotErrorTreatment = true;
	errorName = optarg[0];
	break;
      case 'S':
	gotSystematicsConfig = true;
	systConfigFile = new char[strlen(optarg)+1];
	strcpy(systConfigFile,optarg);
	break;
      default:
	std::cout << "usage: -r responseFile [-d <dataFile>]   \n";
	abort ();
      }

  // Open Data distribution

  if (!gotResponse) {
    std::cout << "You need to provide a file with the response \n";
    return 0;
  }


  if (gotSystematicsConfig) {
    SystematicsManager * sysManager = SystematicsManager::GetInstance();
    //    std::string sysFileName = systConfigFile;
    sysManager->Setup(systConfigFile);
  }

  string algorithm;
  
  if (gotUnfoldingAlgo) {
    algorithm = unfoldingAlgo;
  } else {
    algorithm = "Bayes";
  }

  // Set Unfolding error type

  //  RooUnfold::kCovariance;
  //  RooUnfold::kErrors;
  //  RooUnfold::kCovToy;

  RooUnfold::ErrorTreatment errorType = RooUnfold::kCovToy;

  if (gotErrorTreatment) { 

    cout << "Eroor Name = [" << errorName << "] \n";
    if (errorName == 'E') {
      errorType = RooUnfold::kErrors;
    } else if (errorName == 'C') {
      errorType = RooUnfold::kCovariance;
    } else if (errorName == 'T') {
      errorType = RooUnfold::kCovToy;
    } else {
      cout << "UNDEFINED UNFOLDING ERROR TYPE = [" 
	   << errorName << "] \n";
      return -1;
    }

  }

  string variable;

  if (gotVarName) {
    variable = variableName;
  } else {
    variable = "LeadJetPt";
  }


  // 
  // INPUT FILES
  // 
  // 


  cout << "Using algo : " << algorithm << endl;

  // Get it from some file....
  // Now check it from MC: consistency check 
  TFile * fDataInput; // = new TFile("wzJets-full.root","READ");  
  TFile * fUnfoldingMatrix = new TFile(responseFileName, "READ");
  TFile * fTauFraction     = new TFile(responseFileName, "READ");
  std::vector<TFile *> backgroundFiles;
  std::vector<TString> backgroundNames;
  // Data 
  if (gotDataFile) {
    fDataInput = new TFile(dataFileName,"READ");
  } else {
    fDataInput = fUnfoldingMatrix;
  }
  // Background(s)
  if (gotBackgroundList) {
    std::ifstream list(bgListName);
    TString name;
    while (list>>name) {
      //      std::cout << "Opening background file: " << name << std::endl;
      backgroundFiles.push_back(new TFile(name));
      backgroundNames.push_back(name);
    } 
  }else if (gotBackgroundFile) {
    backgroundFiles.push_back(new TFile(bgFileName));
  }
  // WZ Tau fraction
  if (gotTauFraction) {
    fTauFraction  = new TFile(tauFractionFileName);
  }

  //  fUnfoldingMatrix = new TFile("wzJets-all_dR05.root","READ");

  TH1D * data[4];
  TH1D * signal[4];
  TH1D * backgrounds[4];
  TH1D * tauFraction[4];
  TH1D * tauCorrection[4];
  TH1D * unfoldedDataDistribution[4];
  TH1D * dSigma[4];
  TH1D * dSigmaDX[4];
  RooUnfoldResponse * response[4];

  // LOOP OVER 4 CHANNELS

  for (int chan=0; chan<4; chan++) {

    // Define Histogram keys
    std::ostringstream histoKey;
    std::ostringstream dataHistoKey;
    std::ostringstream signalHistoKey;
    std::ostringstream bgHistoKey;
    std::ostringstream ftauHistoKey;
    std::ostringstream tauCorrHistoKey;
    std::ostringstream responseKey;
    std::ostringstream resultKey;
    std::ostringstream dsigmaKey;
    std::ostringstream dsigmaDXKey;

    histoKey        <<  variable << "_"    << chan+1;
    dataHistoKey    << "hData" << variable << "_"     << chan+1;
    signalHistoKey  << "hSignal" << variable << "_"     << chan+1;
    bgHistoKey      << "hBgd" << variable << "_"      << chan+1;
    ftauHistoKey    << "hftau" << variable << "_"     << chan+1;
    tauCorrHistoKey << "hftaucorr" << variable << "_" << chan+1;
    resultKey       << "hresult" << variable << "_"   << chan+1;
    dsigmaKey       << "hdsigma" << variable << "_"   << chan+1;
    dsigmaDXKey     << "hdsigmadx" << variable << "_"   << chan+1;

    responseKey   << "response" << variable << "_"    << chan+1;

    // Read data distribution (is Output from Matrix Method)

    std::cout << "READING HISTO KEY: " << histoKey << std::endl;

    TH1D * dh = (TH1D*) fDataInput->Get(histoKey.str().c_str()); // ->Clone(dataHistoKey.str().c_str());
    if (!dh) {
      std::cout << "Found no data histogram named : " << histoKey.str() << std::endl;
      return -1;
    }

    data[chan] = (TH1D*) dh->Clone(dataHistoKey.str().c_str());

    data[chan]->Sumw2();
    signal[chan] = (TH1D*) data[chan]->Clone(signalHistoKey.str().c_str());
    
    //systematics for backgounds
    //    double MC_xs_uncertainty[4]={1, 0.5, 0.15, 0.15};
    // double

    // Read backgrounds and subtract them from data

    //substract backgounds with systematics
    double variation(1);    
    int ZZ_syst=0;
    int Zgamma_syst=0;
    int WV_syst=0;
    int WZZJets_syst=0;
    int ZZZJets_syst=0;
    int WWZJets_syst=0;
    int WWWJets_syst=0;
    int TTWJets_syst=0;
    int TTZJets_syst=0;
    int TTWWJets_syst=0;
    int TTGJets_syst=0;
    int WWGJets_syst=0;
       
    if (gotSystematicsConfig){
      
      SystematicsManager * sysManager = SystematicsManager::GetInstance();
      ZZ_syst=sysManager->GetValue("ZZ");
      Zgamma_syst=sysManager->GetValue("Zgamma");
      WV_syst=sysManager->GetValue("WV");
      WZZJets_syst=sysManager->GetValue("TTGJets");
      ZZZJets_syst=sysManager->GetValue("TTWWJets");
      WWZJets_syst=sysManager->GetValue("TTZJets");
      WWWJets_syst=sysManager->GetValue("TTWJets");
      TTWJets_syst=sysManager->GetValue("WWWJets");
      TTZJets_syst=sysManager->GetValue("WWZJets");
      TTWWJets_syst=sysManager->GetValue("ZZZJets");
      TTGJets_syst=sysManager->GetValue("WZZJets");
      WWGJets_syst=sysManager->GetValue("WWGJets");
    }
    std::cout<<"ZZ_syst"<<ZZ_syst<<std::endl;
    
    for (int ibg=0; ibg<backgroundFiles.size(); ibg++) {
      if (backgroundNames[ibg].Contains("ZZ.root")){
	variation=1+ZZ_syst*0.15;
      }
      if (backgroundNames[ibg].Contains("Zgamma.root")){
	variation=1+Zgamma_syst*0.15;
      }
      if (backgroundNames[ibg].Contains("WV.root")){
	variation=1+WV_syst*0.2;
      }
      if (backgroundNames[ibg].Contains("WZZJets.root")){
	variation=1+WZZJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("ZZZJets.root")){
	variation=1+ZZZJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("WWZJets.root")){
	variation=1+WWZJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("WWWJets.root")){
	variation=1+WWWJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("TTWJets.root")){
	variation=1+TTWJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("TTZJets.root")){
	  variation=1+TTZJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("TTWWJets.root")){
	variation=1+TTWWJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("TTGJets.root")){
	  variation=1+TTGJets_syst*0.5;
      }
      if (backgroundNames[ibg].Contains("WWGJets.root")){
	variation=1+WWGJets_syst*0.5;
      }
      std::cout<<"VAriation: "<<variation<<std::endl;
      if (ibg==0) {
	backgrounds[chan] = (TH1D*) backgroundFiles[ibg]->Get(histoKey.str().c_str())->Clone(bgHistoKey.str().c_str());
	backgrounds[chan]->Scale(variation);
      } else {
	TH1D * hbg = (TH1D*) backgroundFiles[ibg]->Get(histoKey.str().c_str());
	hbg->Scale(variation);
	backgrounds[chan]->Add(hbg);
      }
    }
   
    /*
    for (int ibg=0; ibg<backgroundFiles.size(); ibg++) {
      if (ibg==0) {
	backgrounds[chan] = (TH1D*) 
	      backgroundFiles[ibg]->Get(histoKey.str().c_str())->Clone(bgHistoKey.str().c_str());
      } else {
	TH1D * hbg = (TH1D*) backgroundFiles[ibg]->Get(histoKey.str().c_str());
	backgrounds[chan]->Add(hbg);
      }
    }
    */
  

  signal[chan]->Add(backgrounds[chan],-1);
  
  // Read tau fraction and orrect for it
  
  tauFraction[chan] = (TH1D*) fTauFraction->Get(histoKey.str().c_str())->Clone(ftauHistoKey.str().c_str());
  
  tauCorrection[chan] = GetTauCorrection(tauFraction[chan],tauCorrHistoKey.str());
  
  signal[chan]->Multiply(tauCorrection[chan]);


    // Open RooUnfoldResponse from file

    std::cout << "Fetching response matrix: " << responseKey.str() << std::endl;
    response[chan] = (RooUnfoldResponse * ) fUnfoldingMatrix->Get(responseKey.str().c_str());

    // Do the unfolding
    unfoldedDataDistribution[chan] = Unfold(algorithm  // "Bayes"
					    ,response[chan]
					    ,errorType
					    ,signal[chan], backgrounds[chan]
					    , 5 // kterm
					    , resultKey.str().c_str()
					    , 1 ); // userOverflow
					    // , truth[chan] );
    std::cout << "done with unfolding channel: " << chan << std::endl;

    // And now get the differential cross section: luminosity

    dSigma[chan] = (TH1D* ) unfoldedDataDistribution[chan]->Clone(dsigmaKey.str().c_str());

    dSigma[chan]->Scale(1./LUMINOSITY);

    dSigmaDX[chan] = (TH1D* ) dSigma[chan]->Clone(dsigmaDXKey.str().c_str());

    for (int i=1; i<=dSigmaDX[chan]->GetNbinsX(); i++) {
      double value = dSigmaDX[chan]->GetBinContent(i);
      double error = dSigmaDX[chan]->GetBinError(i);
      double width = dSigmaDX[chan]->GetBinWidth(i);
      double dsdx = value/width;
      double dsdx_err = dsdx*error/value;
      dSigmaDX[chan]->SetBinContent(i,dsdx);
      dSigmaDX[chan]->SetBinError(i,dsdx_err);
    }


    // STILL NEEDS TO DIVIDE WITH BIN WIDTH TO GET DSIGMA/Dx

    std::cout << "!!!!!! not divining with bin width yet!!!!! \n";



  }

  std::cout << "Closing files \n";

//   if (fDataInput)          fDataInput->Close();
//   if (fUnfoldingMatrix)    fUnfoldingMatrix->Close();
//   if (fTauFraction)        fTauFraction->Close();

//   for (int i=0; i<backgroundFiles.size(); i++) {
//     backgroundFiles[i]->Close();
//   }

  // Compute total cross sectin summing over bins

  float sigmaTot[4] = {0,0,0,0.};
  double WZbr[4]={0.03363*0.1075,0.03363*0.1057,0.03366*0.1075,0.03366*0.1057};

  for (int i=0; i<4; i++) {
    sigmaTot[i] = dSigma[i]->Integral();
    std::cout << "integral for channel " << i << " : "
	      << sigmaTot[i] << std::endl;
    std::cout << "inclusivee: "<<sigmaTot[i]/WZbr[i]<<std::endl;
  }


  // 
  TFile * fout = new TFile("unfolding.root","RECREATE");
  fout->cd();
  for (int i=0; i<4; i++) {
    std::cout << "Writing histos: " 
	      << data[i] << "\t" << backgrounds[i] << std::endl;
    data[i]->Write();
    //    truth[i]->Write();
    backgrounds[i]->Write();
    signal[i]->Write();
    tauFraction[i]->Write();
    tauCorrection[i]->Write();
    unfoldedDataDistribution[i]->Write();
    dSigma[i]->Write();
    dSigmaDX[i]->Write();
  }
  fout->Close();



}
