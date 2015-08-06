#include "UnfoldingAnalysis.h"

#include "UnfoldingHistogramFactory.h"

#include "constants.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#define USENORMALIZEDWEIGHTS  true
#define NORMALIZETOLUMINOSITY true


UnfoldingAnalysis::UnfoldingAnalysis(std::string k, WZEvent * e) :
    key(k),
    wzevt(e) 
{
//   std::cout << "Calling Unf. constructor \n";
//   CreateBaseHistos();
//   std::cout << "Calling Init method \n";
//   Init();
//   std::cout << "Called Init method \n";

  useNormalizedWeights = USENORMALIZEDWEIGHTS;
  normalizeToLumi      = NORMALIZETOLUMINOSITY;
};



void UnfoldingAnalysis::CreateBaseHistos() {

  for (int i=0; i<5; i++) {
    std::ostringstream genhistoKey;
    std::ostringstream genhistoTitle;
    genhistoKey << "hGen" << key << "_" << i;

    genhistoTitle << "Gen " << key << "  for channel " << i;
    genHistos[i] = createHistogram(genhistoKey.str(), genhistoTitle.str());
    genHistos[i]->Sumw2();


//     std::ostringstream genxshistoKey;
//     std::ostringstream genxshistoTitle;
//     genxshistoKey << "hGenXs" << key << "_" << i;
//     genxshistoTitle << "Truth XS " << key << "  for channel " << i;
//     genXSHistos[i] = createHistogram(genxshistoKey.str(), 
// 				     genxshistoTitle.str());
//     genXSHistos[i]->Sumw2();


    std::ostringstream recohistoKey;
    std::ostringstream recohistoTitle;
    recohistoKey << "hReco" << key << "_" << i;
    recohistoTitle << "Reco " << key << "  for channel " << i;
    recoHistos[i] = createHistogram(recohistoKey.str(), recohistoTitle.str());
    recoHistos[i]->Sumw2();

    std::ostringstream crecohistoKey;
    std::ostringstream crecohistoTitle;
    crecohistoKey << "hControlReco" << key << "_" << i;
    crecohistoTitle << "CR:" << " Reco " << key << " for channel " << i;
    controlRecoHistos[i] = createHistogram(crecohistoKey.str(), crecohistoTitle.str());
    controlRecoHistos[i]->Sumw2();
    std::ostringstream crgenhistoKey;
    std::ostringstream crgenhistoTitle;
    crgenhistoKey << "hControlGen" << key << "_" << i;
    crgenhistoTitle << "CR-GEN:" << key << "  for channel " << i;
    controlGenHistos[i] = createHistogram(crgenhistoKey.str(), crgenhistoTitle.str());
    controlGenHistos[i]->Sumw2();

    // NOW WE NEED TO FILL THESE DISTRIBUTIONS...

    response[i]   = new RooUnfoldResponse(recoHistos[i],
					  genHistos[i]);
    std::ostringstream purityhistoKey;
    std::ostringstream purAllhistoKey;
    std::ostringstream purityhistoTitle;
    std::ostringstream purAllhistoTitle;
    purityhistoKey << "hPurity" << key << "_" << i;
    purityhistoTitle << "Purity for " << key << " for channel " << i;
    purAllhistoKey << "hPurityAll" << key << "_" << i;
    purAllhistoTitle << "Purity All " << key << "  for channel " << i;

    purityPlot[i]   = createHistogram(purityhistoKey.str(),
				      purityhistoTitle.str());
    purityPlotDenominator[i] = createHistogram(purAllhistoKey.str(),
					       purAllhistoTitle.str());

    std::ostringstream stabilityhistoKey;
    std::ostringstream stabAllhistoKey;
    std::ostringstream stabilityhistoTitle;
    std::ostringstream stabAllhistoTitle;
    stabilityhistoKey << "hStability" << key << "_" << i;
    stabilityhistoTitle << "Stability for " << key << " for channel " << i;
    stabAllhistoKey << "hStabilityAll" << key << "_" << i;
    stabAllhistoTitle << "Stability All " << key << "  for channel " << i;

    stabilityPlot[i]   = createHistogram(stabilityhistoKey.str(),
					 stabilityhistoTitle.str());
    stabilityPlotDenominator[i] = createHistogram(stabAllhistoKey.str(),
						  stabAllhistoTitle.str());


  }

  std::string treeName = "resolution";
  treeName+=key;


  // Create tree to study 
  resolutionTree = new TTree(treeName.c_str(),"REsolution");
  resolutionTree->Branch("genValue", trueValue, "genValue/D");
  resolutionTree->Branch("recoValue", recoValue, "genValue/D");

}

void UnfoldingAnalysis::ApplyLuminosityNormalization(double norm){

  if (!normalizeToLumi) return;

  std::cout << "APPLY LUMI NORMALIZATION TO UNF. HISTOS: "
	    << norm << std::endl;

  double BrZToL    = 0.0337;
  double BrWtoE    = 0.1075;
  double BrWtoMu   = 0.1057;
  double BrWZto3L  = BrZToL*BrWtoE;

  for (int i=0; i<5; i++) {
    if (i>0) {
      genHistos[i]->Scale(norm);
      recoHistos[i]->Scale(norm);
      controlRecoHistos[i]->Scale(norm);
      controlGenHistos[i]->Scale(norm);

      // Copy gen histos to diff xs histos
      std::ostringstream genxshistoKey;
      std::ostringstream genxshistoTitle;
      genxshistoKey << "hGenXs" << key << "_" << i;
      genxshistoTitle << "Truth XS " << key << "  for channel " << i;

      genXSHistos[i] = (TH1D*) genHistos[i]->Clone(genxshistoKey.str().c_str());

      // Now loop over bins and rescale each bin by
      //   1/LUMINOSITY * 1/binWidth * 1/BR
      
      double totalXs = 0;

      for (int k=0; k<=genXSHistos[i]->GetNbinsX()+1; k++) {

	double binxs    = genXSHistos[i]->GetBinContent(k);
	double binxs_err= genXSHistos[i]->GetBinError(k);
	double binWidth = genXSHistos[i]->GetBinWidth(k);
	double scaleFactor = 1./ (LUMINOSITY*BrWZto3L*binWidth);
	double bindxs     = binxs*scaleFactor;
	double bindxs_err = binxs_err*scaleFactor;

	totalXs += bindxs;

	genXSHistos[i]->SetBinContent(k,bindxs);
	genXSHistos[i]->SetBinError(k,bindxs_err);

      }
      std::cout << "Total xs in channel : " << i << " = " << totalXs/BrWZto3L
		<< std::endl;


    }
  }
}

void UnfoldingAnalysis::Finish(TFile * fout) {

  std::cout << "Base Finish \n";

  if (!fout) return; // No output file: nothing to do

  fout->cd();

  for (int i=0; i<5; i++) {
    
    if (i>0) {

      genHistos[i]->Write();
      recoHistos[i]->Write();
      controlRecoHistos[i]->Write();
      controlGenHistos[i]->Write();
      genXSHistos[i]->Write();
      std::ostringstream respkey;
      respkey << "response" << key << "_" << i;
      fout->WriteTObject(response[i],respkey.str().c_str());
    }

    purityPlot[i]->Divide(purityPlotDenominator[i]);
    stabilityPlot[i]->Divide(stabilityPlotDenominator[i]);

    purityPlot[i]->Write();
    stabilityPlot[i]->Write();

  }

  if (resolutionTree) resolutionTree->Write();

}

void UnfoldingAnalysis::FillEvent(bool controlSample) {

  EventAnalysis(controlSample);
  FillPurityStability();

}



void UnfoldingAnalysis::FillPurityStability() {

  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->GetPileupWeight();
  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  double weight = pileUpWeight;

  if (eventPassed
      && wzGenChannel >=0 && wzGenChannel <4
      ) {


    int genBin  = stabilityPlot[0]->FindBin(*trueValue);
    int recoBin = stabilityPlot[0]->FindBin(*recoValue);

//     if (*trueValue<0 || *recoValue<0) 
//       std::cout << "PURSTAB: true Value = " << *trueValue
// 		<< "\t genBin = " << genBin
// 		<< "\t reco Value = " << *recoValue 
// 		<< "\t recoBin = " << recoBin
// 		<< std::endl;

    if (*recoValue>=0.) {
      purityPlotDenominator[wzGenChannel+1]->Fill(*recoValue);
      purityPlotDenominator[0]->Fill(*recoValue);
    }
    if (*trueValue>=0) {
      stabilityPlotDenominator[wzGenChannel+1]->Fill(*trueValue);
      stabilityPlotDenominator[0]->Fill(*trueValue);
    }

    if (genBin == recoBin) {
      if (*recoValue>=0.) {
	purityPlot[wzGenChannel+1]->Fill(*recoValue);
	purityPlot[0]->Fill(*recoValue);
      }
      if (*trueValue>=0) {
	stabilityPlot[wzGenChannel+1]->Fill(*trueValue);
	stabilityPlot[0]->Fill(*trueValue);
      }
    }
    //    resolutionTree->Fill();
  }
}

double UnfoldingAnalysis::GetGenWeight() {

  double genWeight     = wzevt->GetPileupWeight();
  return genWeight;

}


double UnfoldingAnalysis::GetRecoWeight() {

  double genWeight     = wzevt->GetPileupWeight();
  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  float mcEfficiency   = wzevt->GetMCWeight();

  double recoWeight = genWeight*mcEfficiency;
  return recoWeight;

}



UnfoldingLeadingJetPt::UnfoldingLeadingJetPt(WZEvent * e) :
  UnfoldingAnalysis("LeadingJetPt", e)
    
{
  std::cout << "Entered constructor of Jet Pt unf. analysis \n";
  trueValue = &leadingGenJetPt;
  recoValue = &leadingRecoJetPt;

  CreateBaseHistos();
  std::cout << "Calling Init method \n";
  std::cout << "Entered Jet Pt Init method \n";
  Init();

};


void UnfoldingLeadingJetPt::Init() {


  //  CreateBaseHistos();
  //  std::cout << "Calling Init method \n";
  //  std::cout << "Entered Jet Pt Init method \n";
  for (int i=0; i<5; i++) {
    std::ostringstream genhistoKey;
    std::ostringstream genhistoTitle;
    genhistoKey << "hnGenJets_" << i;
    genhistoTitle << "Nr. of Gen jets for channel " << i;
    hnGenJets[i] = new TH1D (genhistoKey.str().c_str(), 
			     genhistoTitle.str().c_str(), 
			     6, -0.5, 5.5);

    std::ostringstream recohistoKey;
    std::ostringstream recohistoTitle;
    recohistoKey << "hnRecoJets_" << i;
    recohistoTitle << "Nr. of reco jets for channel " << i;
    hnRecoJets[i] =  new TH1D (recohistoKey.str().c_str(), 
			       recohistoTitle.str().c_str(), 
			       6, -0.5, 5.5);
  }

  resolutionTree->Branch("recoJetPhi", &leadingRecoJetPhi, "recoJetPhi/D");
  resolutionTree->Branch("genJetPhi", &leadingGenJetPhi, "genJetPhi/D");
  resolutionTree->Branch("recoJetEta", &leadingRecoJetEta, "recoJetEta/D");
  resolutionTree->Branch("genJetEta", &leadingGenJetEta, "genJetEta/D");
  resolutionTree->Branch("recoJetDRZl", &leadingRecoJetDRZl, "recoJetDRZl/D");
  resolutionTree->Branch("recoJetDRWl", &leadingRecoJetDRWl, "recoJetDRWl/D");


  // Create tree to study resolution




TTree* wz = new TTree();
float  WZgen_ptZ = 0.;
int  WZpass_3e = 0.;
float  WZweight_total = 0.;
float  WZweight_reweighted = 0.;
wz->Branch("WZgen_ptZ", &WZgen_ptZ, "WZgen_ptZ/F");
wz->Branch("WZpass_3e", &WZpass_3e, "WZpass_3e/I");
wz->Branch("WZweight_total", &WZweight_total, "WZweight_total/F");
wz->Branch("WZweight_reweighted", &WZweight_reweighted,
"WZweight_reweighted/F");



}

void UnfoldingLeadingJetPt::EventAnalysis(bool controlSample) {


  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->GetPileupWeight();
  float mcEfficiency=wzevt->GetMCWeight();

  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  double genWeight     = pileUpWeight;
  double recoWeight = pileUpWeight*mcEfficiency;
  //  double recoWeight = genWeight;


  // 
  // Count gen jets
  // 
  
  int nGenJets          = 0;
  leadingGenJetPt       = -9999.;
  int leadingGenJet     = -1;
  
  for (int i=0; i < wzevt->genJets.size(); i++) {

    // should be away from the gen leptons
    
    bool closeToLepton = false;
    for (int igl=0; igl<wzevt->genLeptons.size() ; igl++) {
      double dR = wzevt->genJets[i].DeltaR(wzevt->genLeptons[igl]);
      if (dR < 0.5) {
	//	  std::cout << "genJet close to lepton: " << dR << std::endl;
	closeToLepton = true;
      }
    }
    
    if (closeToLepton) continue;
    
    double jetPt  = wzevt->genJets[i].Pt();
    double jetEta = wzevt->genJets[i].Eta();
    //      cout << "Jet " << i << " : Pt = " << jetPt << "\t eta = " << jetEta << std::endl;
    if (jetPt > 30. && fabs(jetEta) < 2.5) {
      nGenJets++;
      if (jetPt > leadingGenJetPt) {
	leadingGenJetPt = jetPt;
	leadingGenJet = i;
      }
    }
  }
  if (leadingGenJet>=0) {
    leadingGenJetPhi = wzevt->genJets[leadingGenJet].Phi();
    leadingGenJetEta = wzevt->genJets[leadingGenJet].Eta();
  }

  if (wzGenChannel >=0 && wzGenChannel <4) {
    hnGenJets[wzGenChannel+1]->Fill(nGenJets,GetGenWeight());
  }


  /////////////////////////////////
  // Look at RECO jets
  

  int nRecoJets = 0;
  leadingRecoJetPt = -9999.;
  leadingRecoJetPhi = -9999.;
  leadingRecoJetEta = -9999.;
  leadingRecoJetDRZl = -9999.;
  leadingRecoJetDRWl = -9999.;

  int leadingRecojet = -1;
  
  if (eventPassed) {
    for (int i=0; i<wzevt->recoJets.size(); i++) {
      
      if (wzevt->recoJets[i].Pt() > 30 && fabs(wzevt->recoJets[i].Eta()) < 2.5) {
	
	// Is this jet close to a reco lepton: skip if yes
	bool closeToLepton = false;
	float drMin = 3.;
	for (int il=0; il<wzevt->leptons.size(); il++) {
	  if (wzevt->recoJets[i].DeltaR(wzevt->leptons[il])<0.5) {
	    closeToLepton = true;
	  }
	}
	if (closeToLepton) continue;
	
	nRecoJets++;
	if (wzevt->recoJets[i].Pt() > leadingRecoJetPt) {
	  leadingRecoJetPt = wzevt->recoJets[i].Pt();
	  leadingRecojet = i;
	}
      }
    }
    if (leadingRecoJetPt >=0) {
	  leadingRecoJetPhi = wzevt->recoJets[leadingRecojet].Phi();
	  leadingRecoJetEta = wzevt->recoJets[leadingRecojet].Eta();
	  double drzl1 = wzevt->ZLepton(0)->DeltaR(wzevt->recoJets[leadingRecojet]);
	  double drzl2 = wzevt->ZLepton(1)->DeltaR(wzevt->recoJets[leadingRecojet]);
	  double drwl  = wzevt->WLepton()->DeltaR(wzevt->recoJets[leadingRecojet]);
	  leadingRecoJetDRZl = TMath::Min(drzl1,drzl2);
	  leadingRecoJetDRWl = drwl;

    }

    if (wzGenChannel >=0 && wzGenChannel <4) {
      hnRecoJets[wzGenChannel+1]->Fill(nRecoJets); 
    }
  }
  
  // Fill Unfolding matrix for Jet Pt spectrum

  // fill ngen vs nreco jets

  // Check that it is the MC channel I want to look at
  if ( !controlSample ) {
    //    if (wzevt->MZ>71. && wzevt->MZ<111.) {
    if (wzevt->PassesGenCuts()) {
      
      //      hnGenJets->Fill(nGenJets);
      if (nGenJets>0) {
	if (wzGenChannel >=0 && wzGenChannel <4) {
	  (genHistos[wzGenChannel+1])->Fill(leadingGenJetPt,GetGenWeight());
	}
      }
      
      //      hnrecoJets->Fill(nRecoJets); 
      if (nRecoJets > 0 ) {
	if (wzGenChannel >=0 && wzGenChannel <4) {
	  (recoHistos[wzGenChannel+1])->Fill(leadingRecoJetPt,GetRecoWeight());
	}
      }
      
      for (int ch=1; ch<5; ch++) {
	
	if (ch != wzGenChannel+1 ) continue;
	
	if (eventPassed) {
	  if (ch!= recoChannel) {
	    std::cout << "Mismatch in reco & Gen Channel: Gen " << ch
		      << "\t Reco: " << recoChannel << std::endl;
	    std::cout<<GetRecoWeight()<<std::endl;
	  }
	  
	}
	// Fake
	//	if (nGenJets<1 && nRecoJets>0) responseJetPt[recoChannel-1]->Fake(leadingRecoJetPt, weight);
	if (nGenJets<1 && nRecoJets>0) response[ch]->Fake(leadingRecoJetPt, 
							  GetRecoWeight());
	
	// Miss
	if (nGenJets>0 && nRecoJets<1) response[ch]->Miss(leadingGenJetPt, 
							  GetGenWeight());
	
	
	// Fill
	if (nGenJets>0 && nRecoJets>0) {
	  //	  response[ch]->Fill(leadingRecoJetPt,leadingGenJetPt, recoWeight);
	  response[ch]->Fill(leadingRecoJetPt,leadingGenJetPt, 
			     GetRecoWeight());
	  if (useNormalizedWeights) {
	    response[ch]->Miss(leadingGenJetPt,
			       GetGenWeight()*(1-mcEfficiency));
	  }
	  if (!eventPassed) std::cout << "ALARM: filling response matrix for not passed event \n";
	}
	
      }
    }
  } else { // FIll CONTROL SAMPLE PLOTS
    //    if (wzevt->MZ>71. && wzevt->MZ<111.) {

    if (wzevt->PassesGenCuts()) {
      if (wzGenChannel >=0 && wzGenChannel <4) {
	if (nGenJets>0) {
	  controlGenHistos[wzGenChannel+1]->Fill(leadingGenJetPt, GetGenWeight());
	}
	if ( nRecoJets>0) {
	  controlRecoHistos[wzGenChannel+1]->Fill(leadingRecoJetPt, GetRecoWeight());
	}
      }
    }
  }

}

void UnfoldingLeadingJetPt::Finish(TFile * fout) {

  std::cout << "Jet Pt Finish \n";

  //END OF TESTTING
  UnfoldingAnalysis::Finish(fout);

  fout->cd();
  for (int i=0; i<4; i++) {
    hnGenJets[i]->Write();
    hnRecoJets[i]->Write();
  }


}


TH1D * UnfoldingLeadingJetPt::createHistogram(std::string s, 
					      std::string title) {


  TH1D * h = UnfoldingHistogramFactory::createLeadingJetHistogram(s,title);

  return h;

//   TH1D * h;

//   if (false) {
//     h = new TH1D(s.c_str(),title.c_str(), 10,0., 500.);

//   } else { // Variable Bin Size

//     std::vector<double> binLimits;
//     double value = 1.;
//     double binSize = 25.;
//     while (value<500.) {
//       binLimits.push_back(value);
//       if (value>200) binSize = 40.;
//       if (value>300) binSize = 50.;
//       value += binSize;
//     }
//     int nBins = binLimits.size() - 1;
//     double *  bins = new double[nBins+1];
//     for (int i=0; i< binLimits.size(); i++) {
//       bins[i] = binLimits[i];
//     }
//     h = new TH1D(key.c_str(),title.c_str(), nBins, bins);

//   }

//   return h;


}

////////////////////////////////////////////////////////////////


UnfoldingZPt::UnfoldingZPt(WZEvent * e) :
  UnfoldingAnalysis("Zpt", e)
    
{
  std::cout << "Entered constructor of Z Pt unf. analysis \n";
  trueValue = &genZPt;
  recoValue = &recoZPt;

  CreateBaseHistos();
  std::cout << "Calling Init method \n";
  std::cout << "Entered Jet Pt Init method \n";
  Init();

};


void UnfoldingZPt::Init () {

  if (!resolutionTree) return;

  resolutionTree->Branch("recoZPhi", &recoZPhi, "recoZPhi/D");
  resolutionTree->Branch("genZPhi",  &genZPhi,  "genZPhi/D");
  resolutionTree->Branch("recoZEta", &recoZEta, "recoZEta/D");
  resolutionTree->Branch("genZEta",  &genZEta,  "genZEta/D");

}



TH1D * UnfoldingZPt::createHistogram(std::string s, 
				     std::string title) {

  TH1D * h = UnfoldingHistogramFactory::createZPtHistogram(s, title);


  //   TH1D * h = new TH1D(s.c_str(),title.c_str(), 20,0., 400.);
  return h;


}

void UnfoldingZPt::EventAnalysis(bool isControlSample) {

  genZPt = wzevt->PtZ;

  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->GetPileupWeight();
  float mcEfficiency = wzevt->GetMCWeight();
  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  //  double weight = pileUpWeight;
  double genWeight     = pileUpWeight;
  double recoWeight = pileUpWeight*mcEfficiency;


  if (eventPassed) {
    recoZPt  = wzevt->SelectedZPt();
    TLorentzVector zp4 = wzevt->SelectedZP4();
    recoZPhi = zp4.Phi();
    recoZEta = zp4.Eta();
  }


  if (wzevt->PassesGenCuts()) {

  if (wzGenChannel >=0 && wzGenChannel <4) {
    if (isControlSample) {
      controlGenHistos[wzGenChannel+1]->Fill(wzevt->PtZ,GetGenWeight());      
    } else {
      genHistos[wzGenChannel+1]->Fill(wzevt->PtZ, GetGenWeight());      
    }
  }
  if (eventPassed) {
    if (recoZPt>=0.) {
      if (wzGenChannel >=0 && wzGenChannel <4) {
	if (isControlSample) {
	  controlRecoHistos[wzGenChannel+1]->Fill(recoZPt,GetRecoWeight());
	} else {
	  recoHistos[wzGenChannel+1]->Fill(recoZPt, GetRecoWeight());
	}
      }
    } else {
      std::cout << "Z Pt negative: " << recoZPt << std::endl;

    }
  }

  if (!isControlSample 
      && wzevt->PassesGenCuts() 
      && wzGenChannel >=0 && wzGenChannel <4) {
    
    // Fake: passes reco but not gen: DOES NOT APPLY HERE
    //    if (nGenJets<1 && nRecoJets>0) response[ch]->Fake(leadingRecoJetPt, weight);
    
    // Miss
    if (!eventPassed ) response[wzGenChannel+1]->Miss(genZPt, GetGenWeight());

    // Fill
    if (eventPassed  ) {
      response[recoChannel]->Fill(recoZPt,genZPt, 
				  GetRecoWeight());
      //				  genWeight*mcEfficiency);
      if (useNormalizedWeights) {
	response[wzGenChannel+1]->Miss(genZPt, 
				       GetGenWeight()*(1-mcEfficiency));
      }
      if (!eventPassed) std::cout << "ALARM: filling response matrix for not passed event \n";
    }
  }

  }

}

UnfoldingNjets::UnfoldingNjets(WZEvent * e) :
  UnfoldingAnalysis("Njets", e)
    
{
  std::cout << "Entered constructor of Z Pt unf. analysis \n";
  //  trueValue=0;
  //recoValue=0;
  trueValue = &(nGenJets);
  recoValue = &(nRecoJets);

  CreateBaseHistos();
  //std::cout << "Calling Init method \n";
  //std::cout << "Entered Jet Pt Init method \n";
  //Init();

};

TH1D * UnfoldingNjets::createHistogram(std::string s, 
				     std::string title) {

  TH1D * h = UnfoldingHistogramFactory::createNjetsHistogram(s, title);


  return h;

}


void UnfoldingNjets::EventAnalysis(bool controlSample) {


  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->GetPileupWeight();
  float mcEfficiency=wzevt->GetMCWeight();

  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  double genWeight     = pileUpWeight;
  double recoWeight = pileUpWeight*mcEfficiency;
  //  double recoWeight = genWeight;


  // 
  // Count gen jets
  // 
  
  //  int nGenJets          = 0;
  nGenJets          = 0;

  for (int i=0; i < wzevt->genJets.size(); i++) {

    // should be away from the gen leptons
    
    bool closeToLepton = false;
    for (int igl=0; igl<wzevt->genLeptons.size() ; igl++) {
      double dR = wzevt->genJets[i].DeltaR(wzevt->genLeptons[igl]);
      if (dR < 0.5) {
	closeToLepton = true;
      }
    }
    
    if (closeToLepton) continue;
    
    double jetPt  = wzevt->genJets[i].Pt();
    double jetEta = wzevt->genJets[i].Eta();

    if (jetPt > 30. && fabs(jetEta) < 2.5) {
      nGenJets++;
    }
  }
  /*
  if (wzGenChannel >=0 && wzGenChannel <4) {
    hnGenJets[wzGenChannel+1]->Fill(nGenJets,GetGenWeight());
  }
  */

  /////////////////////////////////
  // Look at RECO jets
  

  //int nRecoJets = 0;
  nRecoJets = 0;

  if (eventPassed) {
    for (int i=0; i<wzevt->recoJets.size(); i++) {
      
      if (wzevt->recoJets[i].Pt() > 30 && fabs(wzevt->recoJets[i].Eta()) < 2.5) {
	
	// Is this jet close to a reco lepton: skip if yes
	bool closeToLepton = false;
	float drMin = 3.;
	for (int il=0; il<wzevt->leptons.size(); il++) {
	  if (wzevt->recoJets[i].DeltaR(wzevt->leptons[il])<0.5) {
	    closeToLepton = true;
	  }
	}
	if (closeToLepton) continue;
	
	nRecoJets++;
      }
    }
    /*        
    if (wzGenChannel >=0 && wzGenChannel <4) {
      hnRecoJets[wzGenChannel+1]->Fill(nRecoJets, GetRecoWeight()); 
    }
    */
  }
  //  trueValue= &nGjets;
  //  recoValue= &nRjets;
  
  // Fill Unfolding matrix for N of jets spectrum
  if ( !controlSample ) {
    if (wzevt->PassesGenCuts()) {
      
      //if (nGenJets>0) {
	if (wzGenChannel >=0 && wzGenChannel <4) {
	  (genHistos[wzGenChannel+1])->Fill(nGenJets,GetGenWeight());
	}
	//}
      
      //      hnrecoJets->Fill(nRecoJets); 
      if (wzGenChannel >=0 && wzGenChannel <4) {
	  (recoHistos[wzGenChannel+1])->Fill(nRecoJets,GetRecoWeight());
	}
      }

      for (int ch=1; ch<5; ch++) {
	
	if (ch != wzGenChannel+1 ) continue;
	
	if (eventPassed) {
	  if (ch!= recoChannel) {
	    std::cout << "Mismatch in reco & Gen Channel: Gen " << ch
		      << "\t Reco: " << recoChannel << std::endl;
	  }
	}
      }
      if (!controlSample 
	  && wzevt->PassesGenCuts() 
	  && wzGenChannel >=0 && wzGenChannel <4) {
	
	// Fake: passes reco but not gen: DOES NOT APPLY HERE
	//    if (nGenJets<1 && nRecoJets>0) response[ch]->Fake(leadingRecoJetPt, weight);
	
	// Miss
	if (!eventPassed ) response[wzGenChannel+1]->Miss(nGenJets, GetGenWeight());
	
	// Fill
	if (eventPassed  ) {
	  response[recoChannel]->Fill(nRecoJets, nGenJets, 
				  GetRecoWeight());
	  //				  genWeight*mcEfficiency);
	  if (useNormalizedWeights) {
	    response[wzGenChannel+1]->Miss(nGenJets, 
					   GetGenWeight()*(1-mcEfficiency));
	  }
	  if (!eventPassed) std::cout << "ALARM: filling response matrix for not passed event \n";
	}
      }	
  } else { //FILL CONTROL SAMPLE PLOTS
    if (wzevt->PassesGenCuts()) {
      if (wzGenChannel >=0 && wzGenChannel <4) {
	if (nGenJets>=0) {
	  controlGenHistos[wzGenChannel+1]->Fill(nGenJets, GetGenWeight());
	}
	if ( nRecoJets>=0) {
	  controlRecoHistos[wzGenChannel+1]->Fill(nRecoJets, GetRecoWeight());
	}
      }
    }
  }
}

  


			 


