#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TSVDUnfold.h"


#include "RooUnfold.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"


#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#define DEBUG_LEVEL 1
#define NR_CHANNELS 4


using namespace std;

TH1D* Unfold(string unfAlg, 
	     string key,
	     RooUnfoldResponse* response, 
	     std::vector<TObject*> * outputObjects,
	     TH1D* hData, 
	     TH1D* hSumBG, 
	     int Kterm, 
	     string hOutName, 
	     bool useOverFlow,
	     TH1D * h_mcTruth=0)
{

  if (useOverFlow) response->UseOverflow();
  RooUnfold* RObject;
  TH1D * hDataClone = (TH1D*) hData->Clone();
  hDataClone->Add(hSumBG, -1);

  std::cout << "ALGO: [" << unfAlg << "]\n";
  if (unfAlg == "SVD") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kSVD,   response, hDataClone, Kterm);
  } else   if (unfAlg == "Bayes") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kBayes, response, hDataClone, Kterm);
  } else   if (unfAlg == "Invert") {
    RObject = (RooUnfold*) RooUnfold::New( RooUnfold::kInvert, response, hDataClone);
  } else {
    std::cout << "Unknown unfolding algorithm " << unfAlg << std::endl;
    return 0;
  }
  RObject->SetVerbose(0);

  // Define error treatment to use
  //  RooUnfold::ErrorTreatment doerror = RooUnfold::kCovariance;//RooUnfold::kCovariance;
  // RooUnfold::kCovToy;
  // RooUnfold::kCovariance;
  // 
  RooUnfold::ErrorTreatment doerror = RooUnfold::kCovToy; // RooUnfold::kCovariance;

  TH1D* hCorrected = (TH1D*) RObject->Hreco(doerror);
  hCorrected->SetName(hOutName.c_str());


  if (DEBUG_LEVEL>0) {
    for (int i=0; i<=hCorrected->GetNbinsX()+1; i++) {
      std::cout << "Bin " << i << "\t Content = " << hCorrected->GetBinContent(i)
		<< "\t +/- " << hCorrected->GetBinError(i) << std::endl;

    }
  }

  RObject->PrintTable(cout, h_mcTruth, (RooUnfold::ErrorTreatment)doerror);

  // Returns d vector 
  if (unfAlg == "SVD") {
    TSVDUnfold *myTSVD = (TSVDUnfold*) (( RooUnfoldSvd*) RObject)->Impl();

    TH1D *svVector = myTSVD->GetSV();
  
    TH1D *dVector = myTSVD->GetD();

    TString dname = "dd";
    TString diname = "ddi";
    TString dzname = "ddz";
    dname += key.c_str();
    diname += key.c_str();
    dzname += key.c_str();
    dVector->SetName(dname);
    TH1D *diVector = (TH1D*) dVector->Clone(diname);
    TH1D *dzVector = (TH1D*) dVector->Clone(dzname);

    double tau = svVector->GetBinContent(Kterm+1);

    int nbin = h_mcTruth->GetNbinsX();
    for (int i=1; i < nbin; i++) {
      double Si = svVector->GetBinContent(i);  
      double scale = (Si*Si)/((Si*Si)+(tau*tau));
      double di = dVector->GetBinContent(i);
      
      dzVector->SetBinContent(i, di*scale);
      diVector->SetBinContent(i, di);
    }

    if (outputObjects) {
      outputObjects->push_back(dVector);
      outputObjects->push_back(diVector);
      outputObjects->push_back(dzVector);
    }

  }

  
  return hCorrected;
}
  


int main(int argc, char **argv) {


  char * responseFileName(0);
  char * dataFileName(0);
  char * unfoldingAlgo(0);
  char * variableName(0);
  char * dataHistoName(0);
  char * outputFileName(0);

  bool gotDataFile = false;
  bool gotResponse = false;
  bool gotUnfoldingAlgo = false;
  bool gotVarName  = false;
  bool useControlSample = false;
  bool useToySample = false;
  bool gotKterm    = false;
  int  inputKterm = -1;

  char c;

  while ((c = getopt (argc, argv, "r:d:a:v:k:o:N:CT")) != -1)
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
      case 'N':
	dataHistoName = new char[strlen(optarg)+1];
	strcpy(dataHistoName,optarg);
	break;
      case 'o':
	outputFileName = new char[strlen(optarg)+1];
	strcpy(outputFileName,optarg);
	break;
      case 'k':
	gotKterm = true;
	inputKterm = atoi(optarg);
	break;
      case 'C':
	useControlSample = true;
	break;
      case 'T':
	useToySample = true;
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
  int kTerm = 5;
  if (gotKterm) {
    std::cout << "K-term input = " << inputKterm << std::endl;
    if (inputKterm>0) kTerm = inputKterm;
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

  std::vector<TObject *> objectsToWrite;

  RooUnfoldResponse * response[4];

//   TString dataPlotBase = "hControlJetPt_";
//   TString truthPlotBase = "hControlGenJetPt_";
//   TString responseBase = "hControlGenJetPt_";

//  for (int chan=0; chan<4; chan++) {
  for (int chan=0; chan<NR_CHANNELS; chan++) {

    //  TH1D * data = (TH1D*) fDataInput->Get("hControlJetPt_1");

    TString key = "bla";
    key += chan;
    TString channel = chan+1;
    std::cout << "key: " << key 
	      << " channel: " << channel
	      << std::endl;

    std::ostringstream tempstr;
    tempstr << "_" << chan;
    std::string identifier = tempstr.str();
    std::cout << "IDDDD : " << identifier << std::endl;

    std::ostringstream dataHistoKey;
    std::ostringstream truthHistoKey;
    std::ostringstream responseKey;
    std::ostringstream backgroundKey;


    if (dataHistoName) {
      dataHistoKey  << dataHistoName  << "_"    << chan+1;      
      truthHistoKey << "hGen" << variable << "_" << chan+1;
    } else if (useControlSample) {
      dataHistoKey  << "hControlReco" << variable << "_"    << chan+1;
      truthHistoKey << "hControlGen" << variable << "_" << chan+1;
    } else if (useToySample) {
      dataHistoKey  << "hToyReco" << variable << "_"    << chan+1;
      truthHistoKey << "hGen" << variable << "_" << chan+1;
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
    //    truth[chan] = (TH1D*) fUnfoldingMatrix->Get(truthHistoKey.str().c_str())->Clone(truthKey.str().c_str());
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
					    , identifier
					    ,response[chan]
					    , &objectsToWrite
					    ,data[chan], backgrounds[chan]
					    , kTerm // kterm
					    , resultKey.str().c_str()
					    , 1  // userOverflow
					    , truth[chan] );

    

    // Study vs nr of iterations

    for (int nIter=1; nIter<10; nIter++) {

      std::ostringstream outKey;
      outKey << resultKey.str() << "_k" << nIter;
      std::cout << "key for iteration : " << outKey.str() << std::endl;

      std::ostringstream resKey;
      resKey << outKey.str() << "_resid";

      std::ostringstream errKey;
      errKey << outKey.str() << "_err";

      std::ostringstream relerrKey;
      relerrKey << outKey.str() << "_relerr";


      // Do the unfolding
      TH1D * hResultVsIter = Unfold(algorithm  // "Bayes"
				    , identifier
				    ,response[chan]
				    , 0  // point to vector of objects to write
				    ,data[chan], backgrounds[chan]
				    , nIter // kterm
				    , outKey.str().c_str()
				    , 1  // userOverflow
				    , truth[chan] );

      TH1D* residual   = (TH1D*) hResultVsIter->Clone(resKey.str().c_str());
      TH1D* errors     = (TH1D*) hResultVsIter->Clone(errKey.str().c_str());
      TH1D* errors_rel = (TH1D*) hResultVsIter->Clone(relerrKey.str().c_str());

      errors->Reset();
      errors_rel->Reset();

      int nbins = errors->GetNbinsX();
      for (int i=0; i <= nbins+1; i++) {
	double x =   hResultVsIter->GetBinContent(i);
	double dx =  hResultVsIter->GetBinError(i);

	errors->SetBinContent(i,dx);
	errors_rel->SetBinContent(i,dx/x);


      }


      residual->Add(truth[chan],-1);
      residual->Divide(hResultVsIter);

      objectsToWrite.push_back(hResultVsIter);
      objectsToWrite.push_back(residual);
      objectsToWrite.push_back(errors);
      objectsToWrite.push_back(errors_rel);

      // We also need:
      // - histo with bin errors
      // - histo with  error / content


    }


  }

  TString outputFile;

  if (outputFileName == 0) {
    outputFile = "unfolding.root";
  } else {
    outputFile = outputFileName;
  }

  TFile * fout = new TFile(outputFile,"RECREATE");
  fout->cd();
  for (int i=0; i<NR_CHANNELS; i++) {
    data[i]->Write();
    truth[i]->Write();
    backgrounds[i]->Write();
    unfoldedDataDistribution[i]->Write();
  }
  for (int i=0; i<objectsToWrite.size(); i++) {
    objectsToWrite[i]->Write();
  }


  fout->Close();


}
