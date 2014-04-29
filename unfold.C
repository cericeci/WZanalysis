#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TLorentzVector.h"

#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"


#include <iostream>
#include <sstream>
#include <string>


using namespace std;

TH1D* Unfold(string unfAlg, RooUnfoldResponse* response, 
	     TH1D* hData, TH1D* hSumBG, int Kterm, 
	     string hOutName, bool useOverFlow,
	     TH1D * h_mcTruth=0)
{

  if (useOverFlow) response->UseOverflow();
  RooUnfold* RObject;
  TH1D * hDataClone = (TH1D*) hData->Clone();
  hDataClone->Add(hSumBG, -1);

  if (unfAlg == "SVD") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kSVD,   response, hDataClone, Kterm);
  } else   if (unfAlg == "Bayes") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kBayes, response, hDataClone, Kterm);
  } else {
    std::cout << "Unknown unfolding algorithm " << unfAlg << std::endl;
    return 0;
  }
  RObject->SetVerbose(0);
  TH1D* hCorrected = (TH1D*) RObject->Hreco();
  hCorrected->SetName(hOutName.c_str());

  RooUnfold::ErrorTreatment doerror = RooUnfold::kCovariance;


  RObject->PrintTable(cout, h_mcTruth, (RooUnfold::ErrorTreatment)doerror);
  
  return hCorrected;
}
  


int main(int argc, char **argv) {


  char * responseFileName(0);
  char * dataFileName(0);
  char * unfoldingAlgo(0);
  char * variableName(0);

  bool gotDataFile = false;
  bool gotResponse = false;
  bool gotUnfoldingAlgo = false;
  bool gotVarName  = false;

  char c;

  while ((c = getopt (argc, argv, "r:d:a:v:")) != -1)
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
      default:
	std::cout << "usage: -r responseFile [-d <dataFile>]   \n";
	abort ();
      }

  // Open Data distribution

  if (!gotResponse) {
    std::cout << "You need to provide a file with the response \n";
    return 0;
  }

  string algorithm;
  
  if (gotUnfoldingAlgo) {
    algorithm = unfoldingAlgo;
  } else {
    algorithm = "Bayes";
  }

  string variable;

  if (gotVarName) {
    variable = variableName;
  } else {
    variable = "LeadJetPt";
  }


  cout << "Using algo : " << algorithm << endl;

  // Get it from some file....
  // Now check it from MC: consistency check 
  TFile * fDataInput; // = new TFile("wzJets-full.root","READ");  
  TFile * fUnfoldingMatrix = new TFile(responseFileName, "READ");

  if (gotDataFile) {
    fDataInput = new TFile(dataFileName,"READ");
  } else {
    fDataInput = fUnfoldingMatrix;
  }
  

  //  fUnfoldingMatrix = new TFile("wzJets-all_dR05.root","READ");

  TH1D * data[4];
  TH1D * truth[4];
  TH1D * backgrounds[4];
  TH1D * unfoldedDataDistribution[4];
  RooUnfoldResponse * response[4];

//   TString dataPlotBase = "hControlJetPt_";
//   TString truthPlotBase = "hControlGenJetPt_";
//   TString responseBase = "hControlGenJetPt_";

  for (int chan=0; chan<4; chan++) {

    //  TH1D * data = (TH1D*) fDataInput->Get("hControlJetPt_1");

    TString key = "bla";
    key += chan;
    TString channel = chan+1;
    std::cout << "key: " << key 
	      << " channel: " << channel
	      << std::endl;

    std::ostringstream dataHistoKey;
    std::ostringstream truthHistoKey;
    std::ostringstream responseKey;
    std::ostringstream backgroundKey;

    if (true) {
      dataHistoKey  << "hControlReco" << variable << "_"    << chan+1;
      truthHistoKey << "hControlGen" << variable << "_" << chan+1;
    } else {
      dataHistoKey  << "hReco" << variable << "_"    << chan+1;
      truthHistoKey << "hGen" << variable << "_" << chan+1;
    }
//     dataHistoKey  << "hRecoLeadJetPt_"    << chan+1;
//     truthHistoKey << "hGenLeadJetPt_" << chan+1;
    responseKey   << "response" << variable << "_"    << chan+1;
    backgroundKey << "backgrounds_"    << chan+1;
    std::ostringstream dataKey;
    std::ostringstream truthKey;
    std::ostringstream resultKey;
    dataKey   << "data_" << chan+1;
    truthKey  << "truth_" << chan+1;
    resultKey  << "result_" << chan+1;

    std::cout << "Reading histogram: " << dataHistoKey.str()
	      << "\t " << truthHistoKey.str()
	      << "\t " << responseKey.str() << std::endl;


    data[chan] = (TH1D*) fDataInput->Get(dataHistoKey.str().c_str())->Clone(dataKey.str().c_str());
    truth[chan] = (TH1D*) fDataInput->Get(truthHistoKey.str().c_str())->Clone(truthKey.str().c_str());

    backgrounds[chan] = (TH1D*) data[chan]->Clone(backgroundKey.str().c_str());

    std::cout << "Integrals -  Data = " << data[chan]->Integral() 
	      << "\t background = " << backgrounds[chan]->Integral() << std::endl; 

    backgrounds[chan]->Scale(0.);

    std::cout << "\t  after rescaling: Integrals -  Data = " << data[chan]->Integral() 
	      << "\t background = " << backgrounds[chan]->Integral() << std::endl;
    
    // Open RooUnfoldResponse from file

    response[chan] = (RooUnfoldResponse * ) fUnfoldingMatrix->Get(responseKey.str().c_str());
    
    std::cout << "Opened response: " << response << std::endl;
    
    // Do the unfolding
    unfoldedDataDistribution[chan] = Unfold(algorithm  // "Bayes"
					    ,response[chan]
					    ,data[chan], backgrounds[chan]
					    , 5 // kterm
					    , resultKey.str().c_str()
					    , 1  // userOverflow
					    , truth[chan] );

  }

  TFile * fout = new TFile("unfolding.root","RECREATE");
  fout->cd();
  for (int i=0; i<4; i++) {
    data[i]->Write();
    truth[i]->Write();
    backgrounds[i]->Write();
    unfoldedDataDistribution[i]->Write();
  }
  fout->Close();


}
