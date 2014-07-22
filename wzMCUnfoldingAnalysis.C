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
#include "UnfoldingAnalysis.h"

#include "WZAnalysis.h"

#include "constants.h"

//#include "WZ.h"
//#include "WZ2012Data.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>


// TH1D * createJetPtHisto(std::string key, std::string title) {

//   TH1D * h = new TH1D(key.c_str(),title.c_str(), 10,0., 500.);
//   return h;

// }


//function declaration goes here:
//
//
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

  TH1D * hZmassMu1         = new TH1D ("hZmassMu1", "hZmassMu1", 100, 60, 120);  
  TH1D * hZmassEl1         = new TH1D ("hZmassEl1", "hZmassEl1", 100, 60, 120);  

  TH2D * hRecoVsGenChannel     = new TH2D ("hRecoVsGenChannel","Reco vs Gen Channel",
					   5, -0.5, 4.5,
					   5, -0.5, 4.5);

//   TH2D * hXChanelTriggerEff    = new TH2D ("hXChanelTriggerEff",
// 					   "XCH Trig. eff: right vs wrong",
// 					   100, 0.,1.,
// 					   100, 0.,1.);

//   TH2D * hXChanelTriggerEff2    = new TH2D ("hXChanelTriggerEff2",
// 					   "XCH Trig. eff: right vs wrong",
// 					   100, 0.92,1.,
// 					   100, 0.92,1.);
//   TH1D * hXChanelTrigEffDiff    = new TH1D ("hXChanelTrigEffDiff",
// 					   "XCH Trig. eff: right - wrong",
// 					    100, -0.02,0.02);

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
  } else {
    std::cout << "Got no input ROOT file: quit \n";
    return 1;
  }
  for (int input=0; input< inputName.size(); input++){
    wz.Add(inputName[input]);
  }
  TTree *wz_tTree=(TTree*)&wz;
  WZEvent *cWZ= new WZEvent(wz_tTree);
  Int_t events= wz_tTree->GetEntries();
  
  std::cout<<"number of events: "<<events << std::endl;

  float yields[5];
  float genYields[5];
  double totalMCYield=0;
  for (int i=0; i<5; i++) {
    yields[i] = 0;
    genYields[i] = 0;
  }
  float nZYield(0.),nWYield(0.);


  //
  // Creating analyses classes that will be run
  //
  UnfoldingLeadingJetPt unfoldJetPt(cWZ);
  //  unfoldJetPt.Init();
  UnfoldingZPt unfoldZPt(cWZ);
  WZAnalysis   genAnalysis(cWZ);
  genAnalysis.Init();


  // Setup selected event lists 
  std::ofstream eventLists[4];
  for (int i=0; i<4; i++) {
    std::ostringstream fileName;
    fileName << "passedEvents_" << i+1;
    std::cout << "file name : "  << fileName.str() << std::endl;
    eventLists[i].open(fileName.str().c_str());
  }

  //  for  (Int_t k = 0; k<events && k<50;k++) {
  for  (Int_t k = 0; k<events; k++) {

    if ( !(k%100000)  ) std::cout << "Processed " << k << " events \n";

    wz_tTree->GetEntry(k);
    cWZ->ReadEvent();

    if (debug) {
      std::cout << "============  Run: " << cWZ->run << "\t Event: " 
		<< cWZ->event << " ======================= \n";
    }

    float pileUpWeight=cWZ->GetPileupWeight();
    int wzGenChannel = cWZ->WZchan;
    
    //rejecting run 201191



    //    if (cWZ->run==201191) continue;
    //    if (!(cWZ->trigger)) continue;

    bool eventPassed   = cWZ->passesSelection();
    FinalState channel = cWZ->GetFinalState();
    PassedSelectionStep selectionLevel = cWZ->GetSelectionLevel();

    if (debug) {
      std::cout << "Event passed: " << eventPassed
		<< "\t step: " << selectionLevel
		<< "\t channel : " << channel << std::endl;
      
      if (eventPassed) std::cout << "PAAAASSSSS " << pileUpWeight << "\n";
    }

    totalMCYield += pileUpWeight;

    if (cWZ->MZ>71.1876  && cWZ->MZ<111.1876) {
      genYields[wzGenChannel] += pileUpWeight;
    }


    if (eventPassed) {
      yields[channel] += pileUpWeight*cWZ->GetMCWeight();

      cWZ->DumpEvent(eventLists[channel-1],1);

    } else {
      if (selectionLevel == passesZSelection) {
	nZYield += pileUpWeight;
      } else  if (selectionLevel == passesWSelection) {
	nWYield += pileUpWeight;
      }
    }

    genAnalysis.EventAnalysis();

    if (k%2) {    
      unfoldJetPt.FillEvent();
      unfoldZPt.FillEvent();
    } else {
      unfoldJetPt.FillEvent(true);
      unfoldZPt.FillEvent(true);
    }

    if (eventPassed) {
      hRecoVsGenChannel->Fill(wzGenChannel,channel);
    }

    // fill ngen vs nreco jets
    double weight = pileUpWeight; // IMPORTANT: here goes, PU, efficiencies / scale factors, ...

  }
  double signalYield = 0;


  std::cout << "Events in sample: " <<  totalMCYield << std::endl;
  std::cout << "Passes Z Selection : " << nZYield 
	    << "  -  " << cWZ->NumZ() << std::endl;
  std::cout << "Passes W Selection : " << nWYield 
	    << "  -  " << cWZ->NumW() << std::endl;
  double totalGenYield = 0.;

  double lumiNormalizedYields[5];
  double lumiNormalization = LUMINOSITY*WZ3LXSECTION/totalMCYield;


  for (int i=0; i<5; i++) {
    lumiNormalizedYields[i] = yields[i]*lumiNormalization;
    totalGenYield += genYields[i];
  }

  for (int i=0; i<5; i++) {

    double Aeff = 1.;
    double Aeff_spanish =1.;
    if (i>0) {
      Aeff = yields[i] / genYields[i-1]; 
      Aeff_spanish = yields[i] / totalGenYield;
    }
    std::cout << "yield  for : " 
	      << i << "\t:\t" << yields[i] 
	      << "\t scaled to lumi : " << lumiNormalizedYields[i] 
	      << "\t gen :" << genYields[i-1]
	      << "\t A*eff = " << Aeff
	      << "\t spanish way  = " << Aeff_spanish
	      << std::endl;
    if (i>0) signalYield+= yields[i];
  }
  std::cout << "Total Signal = " << signalYield << std::endl;
  std::cout << "Total Gen Signal = " << totalGenYield << std::endl;

  cWZ->PrintSummary();

  fout->cd();
  //  fout->Write();
//   hXChanelTriggerEff->Write();
//   hXChanelTriggerEff2->Write();
//   hXChanelTrigEffDiff->Write();

  hRecoVsGenChannel->Write(); 
  genAnalysis.Finish(fout);
  unfoldJetPt.Finish(fout);
  unfoldZPt.Finish(fout);


  fout->Close();
}
