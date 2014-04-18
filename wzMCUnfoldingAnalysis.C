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

//#include "WZ.h"
//#include "WZ2012Data.h"

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

//   TH1D * hnGenJets            = new TH1D ("hnGenJets", "Nr. gen jets", 6, -0.5, 5.5);
//   TH1D * hleadingGenJetPt     = createJetPtHisto("hleadingGenJetPt", "Leading GenJet Pt");

//   TH1D * hnrecoJets            = new TH1D ("hnrecoJets", "Nr. reco jets", 6, -0.5, 5.5);
//   TH1D * hleadingrecoJetPt     = createJetPtHisto("hleadingRecoJetPt", "Leading reco Jet Pt");

  TH2D * hRecoVsGenChannel     = new TH2D ("hRecoVsGenChannel","Reco vs Gen Channel",
					   5, -0.5, 4.5,
					   5, -0.5, 4.5);

  // Control distributions for unfolding

//   TH1D * genJetPtHistos[4];
//   TH1D * recoJetPtHistos[4];

//   TH1D * controlJetPtHistos[4];
//   TH1D * controlGenJetPtHistos[4];

//   for (int i=0; i<4; i++) {
//     std::ostringstream gjhistoKey;
//     std::ostringstream gjhistoTitle;
//     gjhistoKey << "hGenJetPt_" << i+1;
//     gjhistoTitle << "Leading Gen Jet Pt for channel " << i+1;
//     genJetPtHistos[i] = createJetPtHisto(gjhistoKey.str(), gjhistoTitle.str());

//     std::ostringstream rjhistoKey;
//     std::ostringstream rjhistoTitle;
//     rjhistoKey << "hRecoJetPt_" << i+1;
//     rjhistoTitle << "Leading Reco Jet Pt for channel " << i+1;
//     recoJetPtHistos[i] = createJetPtHisto(rjhistoKey.str(), rjhistoTitle.str());

//     std::ostringstream crjhistoKey;
//     std::ostringstream crjhistoTitle;
//     crjhistoKey << "hControlJetPt_" << i+1;
//     crjhistoTitle << "CR: Leading Reco Jet Pt for channel " << i+1;
//     controlJetPtHistos[i] = createJetPtHisto(crjhistoKey.str(), crjhistoTitle.str());
//     std::ostringstream crgjhistoKey;
//     std::ostringstream crgjhistoTitle;
//     crgjhistoKey << "hControlGenJetPt_" << i+1;
//     crgjhistoTitle << "CR-GEN: Leading Gen Jet Pt for channel " << i+1;
//     controlGenJetPtHistos[i] = createJetPtHisto(crgjhistoKey.str(), crgjhistoTitle.str());

//     // NOW WE NEED TO FILL THESE DISTRIBUTIONS...

//   }
  




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
    //  inputName.push_back("/STORE/lucija/latinosTrees/MC_LooseLooseTypeI/latino_074_WZJetsMad.root");
    inputName.push_back("/users/vuko/phan/wz8tev/latinos/CMSSW_5_3_15/src/WWAnalysis/AnalysisStep/test/step3/latinowz-step3-oldJets.root"); 
    //  inputName.push_back("/users/vuko/phan/wz8tev/latinos/CMSSW_5_3_15/src/WWAnalysis/AnalysisStep/test/step3/latinowz-step3-NoElNoMuNoNuJets.root");
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
  for (int i=0; i<5; i++) {
    yields[i] = 0;
    genYields[i] = 0;
  }
  float nZYield(0.),nWYield(0.);


  // UNFOLDING: SETTING UP

//   RooUnfoldResponse *responseJetPt[4];
//   for (int i=0; i<4; i++) {
//     responseJetPt[i]   = new RooUnfoldResponse(hleadingrecoJetPt,
// 					       hleadingGenJetPt);
//   }

  UnfoldingLeadingJetPt unfoldJetPt(cWZ);
  //  unfoldJetPt.Init();

  UnfoldingZPt unfoldZPt(cWZ);

  WZAnalysis   genAnalysis(cWZ);


  // Setup event lists 
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

    float pileUpWeight=cWZ->puW;
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

  std::cout << "Passes Z Selection : " << nZYield 
	    << "  -  " << cWZ->NumZ() << std::endl;
  std::cout << "Passes W Selection : " << nWYield 
	    << "  -  " << cWZ->NumW() << std::endl;
  double totalGenYield = 0.;


  for (int i=0; i<5; i++) totalGenYield += genYields[i];

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
  unfoldJetPt.Finish(fout);
  unfoldZPt.Finish(fout);


  fout->Close();
}
