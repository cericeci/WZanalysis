#define DEBUG  false

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"


#include "RooUnfoldResponse.h"

// Replace this with the new tree
#include "WZEvent.h"
#include "WZAnalysis.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>


TH1D * createJetPtHisto(std::string key, std::string title) {

  TH1D * h = new TH1D(key.c_str(),title.c_str(), 10,0., 500.);
  return h;

}


//function declaration goes here:
//
//
//


///////////////////////////////////////////
//
//  Arguments:
//  [ -i <inputFile> |  -l <listOfInputFiles> ]  [-o <outputFile> \]
//

int main(int argc, char **argv)
{
  using namespace std;

  bool debug=DEBUG;

  char * outputName(0);
  char * inputFileName(0);

  char * fileList;
  bool useInputList=false;

  bool gotInput  = false;
  bool gotOutput = false;
  char c;

  while ((c = getopt (argc, argv, "i:o:l:")) != -1)
    switch (c)
      {
      case 'o':
	gotOutput = true;
	outputName = new char[strlen(optarg)+1];
	strcpy(outputName,optarg);
	break;
      case 'i':
	gotInput = true;
	inputFileName = new char[strlen(optarg)+1];
	strcpy(inputFileName,optarg);
	break;
      case 'l':
	useInputList = true;
	fileList = new char[strlen(optarg)+1];
	strcpy(fileList,optarg);
	break;
      default:
	std::cout << "usage: [-k|-g|-l] [-v] [-b <binWidth>]   -i <input> -o <output> \n";
	abort ();
      }

  // OUTPUT ROOT FILE

  TFile * fout;
  if (gotOutput) {
    fout = new TFile(outputName, "RECREATE");    
  } else {
    fout = new TFile("wzJets-test.root", "RECREATE");
  }

  // INPUT TREES

  std::vector<TString>inputName;
  TChain wz("latino");

  if (useInputList) {
    ifstream list(fileList);
    TString name;
    while (list>>name) {
      cout << "adding name: " << name << endl;
      inputName.push_back(name);
    }
  } else   if (gotInput) {
    inputName.push_back(inputFileName);
  } else {   // Some default file in case there's no input given: PROBABLY AN OLDER VERSION OF THE TREE
    inputName.push_back("/users/vuko/phan/wz8tev/latinos/CMSSW_5_3_15/src/WWAnalysis/AnalysisStep/test/step3/latinowz-step3-oldJets.root"); 
  }
  for (int input=0; input< inputName.size(); input++){
    wz.Add(inputName[input]);
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZEvent *cWZ= new WZEvent(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;

  // SETUP ANALYSIS

  float yields[5];
  float genYields[5];
  for (int i=0; i<5; i++) {
    yields[i] = 0;
    genYields[i] = 0;
  }
  float nZYield(0.),nWYield(0.);

  // Set up some event analysis
  WZAnalysis   genAnalysis(cWZ);


  // Setup event lists of selected events for each channel
  std::ofstream eventLists[4];
  for (int i=0; i<4; i++) {
    std::ostringstream fileName;
    fileName << "acceptedEvents_" << i+1;
    std::cout << "file name : "  << fileName.str() << std::endl;
    eventLists[i].open(fileName.str().c_str());
  }

  //  for  (Int_t k = 0; k<events && k<50;k++) {
  for  (Int_t k = 0; k<events; k++) {

    if ( !(k%100000)  ) std::cout << "Processed " << k << " events \n";

    wz_tTree->GetEntry(k);
    cWZ->ReadEvent();

    float pileUpWeight=cWZ->puW;
    int wzGenChannel = cWZ->WZchan;

    bool eventPassed   = cWZ->passesSelection();
    FinalState channel = cWZ->GetFinalState();
    PassedSelectionStep selectionLevel = cWZ->GetSelectionLevel();

    if (cWZ->MZ>71.1876  && cWZ->MZ<111.1876) {
      genYields[wzGenChannel] += pileUpWeight;
    }

    if (eventPassed) {
      yields[channel] += pileUpWeight;

      eventLists[channel-1] << cWZ->run << "\t" 
			    << cWZ->event  << "\t" 
			    << pileUpWeight
			    << std::endl;
    } else {
      if (selectionLevel == passesZSelection) {
	nZYield += pileUpWeight;
      } else  if (selectionLevel == passesWSelection) {
	nWYield += pileUpWeight;
      }
    }

    genAnalysis.EventAnalysis();

  }

  double signalYield = 0;
  double totalGenYield = 0.;
  for (int i=0; i<5; i++) totalGenYield += genYields[i];
  std::cout << "Passes Z Selection : " << nZYield 
	    << "  -  " << cWZ->NumZ() << std::endl;
  std::cout << "Passes W Selection : " << nWYield 
	    << "  -  " << cWZ->NumW() << std::endl;

  for (int i=0; i<5; i++) {

    double Aeff = 1.;
    double Aeff_spanish =1.;
    if (i>0) {
      Aeff = yields[i] / genYields[i-1]; 
      Aeff_spanish = yields[i] / totalGenYield;
    }
    std::cout << "yield  for : " 
	      << i << "\t:\t" << yields[i] 
	      << "\t gen :" << genYields[i-1]
	      << "\t A*eff = " << Aeff
	      << "\t spanish way  = " << Aeff_spanish
	      << std::endl;
    if (i>0) signalYield+= yields[i];
  }
  std::cout << "Total Signal = " << signalYield << std::endl;
  std::cout << "Total Gen Signal = " << totalGenYield << std::endl;

  cWZ->PrintSummary();

  //  fout->cd();
  //  fout->Write();


  genAnalysis.Finish();

  fout->Close();
}
