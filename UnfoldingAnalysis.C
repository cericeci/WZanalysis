#include "UnfoldingAnalysis.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


UnfoldingAnalysis::UnfoldingAnalysis(std::string k, WZEvent * e) :
    key(k),
    wzevt(e) 
{
//   std::cout << "Calling Unf. constructor \n";
//   CreateBaseHistos();
//   std::cout << "Calling Init method \n";
//   Init();
//   std::cout << "Called Init method \n";
};



void UnfoldingAnalysis::CreateBaseHistos() {

  for (int i=0; i<5; i++) {
    std::ostringstream genhistoKey;
    std::ostringstream genhistoTitle;
    genhistoKey << "hGen" << key << "_" << i;
    genhistoTitle << "Gen " << key << "  for channel " << i;
    genHistos[i] = createHistogram(genhistoKey.str(), genhistoTitle.str());

    std::ostringstream recohistoKey;
    std::ostringstream recohistoTitle;
    recohistoKey << "hReco" << key << "_" << i;
    recohistoTitle << "Reco " << key << "  for channel " << i;
    recoHistos[i] = createHistogram(recohistoKey.str(), recohistoTitle.str());

    std::ostringstream crecohistoKey;
    std::ostringstream crecohistoTitle;
    crecohistoKey << "hControlReco" << key << "_" << i;
    crecohistoTitle << "CR:" << " Reco " << key << " for channel " << i;
    controlRecoHistos[i] = createHistogram(crecohistoKey.str(), crecohistoTitle.str());
    std::ostringstream crgenhistoKey;
    std::ostringstream crgenhistoTitle;
    crgenhistoKey << "hControlGen" << key << "_" << i;
    crgenhistoTitle << "CR-GEN:" << key << "  for channel " << i;
    controlGenHistos[i] = createHistogram(crgenhistoKey.str(), crgenhistoTitle.str());

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
      std::ostringstream respkey;
      respkey << "response" << key << "_" << i;
      fout->WriteTObject(response[i],respkey.str().c_str());
    }

    purityPlot[i]->Divide(purityPlotDenominator[i]);
    stabilityPlot[i]->Divide(stabilityPlotDenominator[i]);

    purityPlot[i]->Write();
    stabilityPlot[i]->Write();

  }

}

void UnfoldingAnalysis::FillEvent(bool controlSample) {

  EventAnalysis(controlSample);
  FillPurityStability();

}



void UnfoldingAnalysis::FillPurityStability() {

  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->puW;
  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  double weight = pileUpWeight;

  if (eventPassed
      && wzGenChannel >=0 && wzGenChannel <4) {

    int genBin  = stabilityPlot[0]->FindBin(*trueValue);
    int recoBin = stabilityPlot[0]->FindBin(*recoValue);

    purityPlotDenominator[wzGenChannel+1]->Fill(*recoValue);
    stabilityPlotDenominator[wzGenChannel+1]->Fill(*trueValue);

    purityPlotDenominator[0]->Fill(*recoValue);
    stabilityPlotDenominator[0]->Fill(*trueValue);

    if (genBin == recoBin) {
      purityPlot[wzGenChannel+1]->Fill(*recoValue);
      stabilityPlot[wzGenChannel+1]->Fill(*trueValue);

      purityPlot[0]->Fill(*recoValue);
      stabilityPlot[0]->Fill(*trueValue);
    }
  }
}


UnfoldingLeadingJetPt::UnfoldingLeadingJetPt(WZEvent * e) :
  UnfoldingAnalysis("LeadJetPt", e)
    
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

}

void UnfoldingLeadingJetPt::EventAnalysis(bool controlSample) {


  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->puW;
  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  double weight = pileUpWeight;


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

  if (wzGenChannel >=0 && wzGenChannel <4) {
    hnGenJets[wzGenChannel+1]->Fill(nGenJets,weight);
  }


  /////////////////////////////////
  // Look at RECO jets
  

  int nRecoJets = 0;
  leadingRecoJetPt = -9999.;
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

    if (wzGenChannel >=0 && wzGenChannel <4) {
      hnRecoJets[wzGenChannel+1]->Fill(nRecoJets); 
    }


  }
  
  // Fill Unfolding matrix for Jet Pt spectrum
  //    std::cout << "DO THE UNFOLDING \n";


  // fill ngen vs nreco jets

  // Check that it is the MC channel I want to look at
  if ( !controlSample ) {
    if (wzevt->MZ>71. && wzevt->MZ<111.) {
      
      //      hnGenJets->Fill(nGenJets);
      if (nGenJets>0) {
	if (wzGenChannel >=0 && wzGenChannel <4) {
	  (genHistos[wzGenChannel+1])->Fill(leadingGenJetPt,weight);
	}
      }
      
      //      hnrecoJets->Fill(nRecoJets); 
      if (nRecoJets > 0 ) {
	if (wzGenChannel >=0 && wzGenChannel <4) {
	  (recoHistos[wzGenChannel+1])->Fill(leadingRecoJetPt,weight);
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
	
	// Fake
	//	if (nGenJets<1 && nRecoJets>0) responseJetPt[recoChannel-1]->Fake(leadingRecoJetPt, weight);
	if (nGenJets<1 && nRecoJets>0) response[ch]->Fake(leadingRecoJetPt, weight);
	
	// Miss
	if (nGenJets>0 && nRecoJets<1) response[ch]->Miss(leadingGenJetPt, weight);
	
	// Fill
	if (nGenJets>0 && nRecoJets>0) {
	  response[ch]->Fill(leadingRecoJetPt,leadingGenJetPt, weight);
	  if (!eventPassed) std::cout << "ALARM: filling response matrix for not passed event \n";
	}
	
      }
    }
  } else { // FIll CONTROL SAMPLE PLOTS

    if (wzevt->MZ>71. && wzevt->MZ<111.) {
      if (wzGenChannel >=0 && wzGenChannel <4) {
	if (nGenJets>0) {
	  controlGenHistos[wzGenChannel+1]->Fill(leadingGenJetPt, weight);
	}
	if ( nRecoJets>0) {
	  controlRecoHistos[wzGenChannel+1]->Fill(leadingRecoJetPt, weight);
	}
      }
    }
  }

}

void UnfoldingLeadingJetPt::Finish(TFile * fout) {

  std::cout << "Jet Pt Finish \n";

  UnfoldingAnalysis::Finish(fout);

  fout->cd();
  for (int i=0; i<4; i++) {
    hnGenJets[i]->Write();
    hnRecoJets[i]->Write();
  }

}


TH1D * UnfoldingLeadingJetPt::createHistogram(std::string s, 
					      std::string title) {

  TH1D * h;

  if (false) {
    h = new TH1D(s.c_str(),title.c_str(), 10,0., 500.);

  } else { // Variable Bin Size

    std::vector<double> binLimits;
    double value = 1.;
    double binSize = 25.;
    while (value<500.) {
      binLimits.push_back(value);
      if (value>200) binSize = 40.;
      if (value>300) binSize = 50.;
      value += binSize;
    }
    int nBins = binLimits.size() - 1;
    double *  bins = new double[nBins+1];
    for (int i=0; i< binLimits.size(); i++) {
      bins[i] = binLimits[i];
    }
    h = new TH1D(key.c_str(),title.c_str(), nBins, bins);

  }

  return h;


}

////////////////////////////////////////////////////////////////


UnfoldingZPt::UnfoldingZPt(WZEvent * e) :
  UnfoldingAnalysis("ZPt", e)
    
{
  std::cout << "Entered constructor of Jet Pt unf. analysis \n";
  trueValue = &genZPt;
  recoValue = &recoZPt;

  CreateBaseHistos();
  std::cout << "Calling Init method \n";
  std::cout << "Entered Jet Pt Init method \n";
  Init();

};

TH1D * UnfoldingZPt::createHistogram(std::string s, 
				     std::string title) {

  TH1D * h = new TH1D(s.c_str(),title.c_str(), 20,0., 400.);
  return h;


}

void UnfoldingZPt::EventAnalysis(bool isControlSample) {

  genZPt = wzevt->PtZ;

  bool eventPassed = (wzevt->GetSelectionLevel() == passesFullSelection);
  FinalState recoChannel = wzevt->GetFinalState();
  int wzGenChannel = wzevt->WZchan;
  float pileUpWeight=wzevt->puW;
  // IMPORTANT: here goes, PU, efficiencies / scale factors, ...
  double weight = pileUpWeight;

  if (eventPassed) {
    recoZPt = wzevt->SelectedZPt();
  }


  if (wzevt->PassesGenCuts()) {

  if (wzGenChannel >=0 && wzGenChannel <4) {
    if (isControlSample) {
      controlGenHistos[wzGenChannel+1]->Fill(wzevt->PtZ);      
    } else {
      genHistos[wzGenChannel+1]->Fill(wzevt->PtZ);      
    }
  }
  if (eventPassed) {
    if (recoZPt>=0.) {
      if (isControlSample) {
	controlRecoHistos[wzGenChannel+1]->Fill(recoZPt);
      } else {
	recoHistos[wzGenChannel+1]->Fill(recoZPt);
      }
    } else {
      std::cout << "Z Pt negative: " << recoZPt << std::endl;

    }
  }

  if (wzevt->PassesGenCuts() && 
      wzGenChannel >=0 && wzGenChannel <4) {
    
    // Fake: passes reco but not gen: DOES APPLY HERE
    //    if (nGenJets<1 && nRecoJets>0) response[ch]->Fake(leadingRecoJetPt, weight);
    
    // Miss
    if (!eventPassed ) response[wzGenChannel+1]->Miss(genZPt, weight);

    // Fill
    if (eventPassed  ) {
      response[recoChannel]->Fill(recoZPt,genZPt, weight);
      if (!eventPassed) std::cout << "ALARM: filling response matrix for not passed event \n";
    }
  }

  }

}

