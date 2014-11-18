#include "WZAnalysis.h"

#include <ios>
#include <iostream>
#include <fstream>


WZAnalysis::WZAnalysis(WZEvent * e) {

  wzevt = e;

  totalNrEvents=0;
  for (int i=0; i<NRWDECAYTYPES; i++) {
    wDecaysByType[i] = 0;

  }
  for (int i=0; i<NRZDECAYTYPES; i++) {
    zDecaysByType[i] = 0;
  }
  weirdEventsList = new std::ofstream("weirdEvents.list"); // ,ios::out);
}


void WZAnalysis::Init() {

  // Create tree to study 
  jetResolutionTree = new TTree("jetres","REsolution");

  jetResolutionTree->Branch("genJetPt",  &_genJetPt,  "genJetPt/F");
  jetResolutionTree->Branch("genJetPhi", &_genJetPhi, "genJetPhi/F");
  jetResolutionTree->Branch("genJetEta", &_genJetEta, "genJetEta/F");

  jetResolutionTree->Branch("recoJetPt",  &_recoJetPt,  "recoJetPt/F");
  jetResolutionTree->Branch("recoJetPhi", &_recoJetPhi, "recoJetPhi/F");
  jetResolutionTree->Branch("recoJetEta", &_recoJetEta, "recoJetEta/F");

  jetResolutionTree->Branch("dR", &_drRecoGenJet, "dR/F");



}

void WZAnalysis::EventAnalysis() {


  if (wzevt->MZ<71.1876 || wzevt->MZ>111.1876) return;

  totalNrEvents++;

  int iwlepton=-1;
  int izleptons[2] = {-1,-1};
  int nwleptons = 0;
  int nzleptons = 0;

  for (int igl=0; igl<wzevt->genLeptons.size() ; igl++) {
    int bosonId = wzevt->genLeptons[igl].MotherBoson();
    if (abs(bosonId) == 24) {
      iwlepton = igl;
      nwleptons++;
    }
    if (abs(bosonId) == 23) {
      if (nzleptons<2) {
	izleptons[nzleptons] = igl;
	nzleptons++;
      }
    }
  }
  if (nwleptons>3) {
    std::cout << "Found more than one W lepton : " << nwleptons
	      << "\t id = " << wzevt->genLeptons[iwlepton].Id()
	      << "\t boson id = " << wzevt->genLeptons[iwlepton].MotherBoson()
	      << "\t from tau = " << wzevt->genLeptons[iwlepton].ComesFromTau()
	      << std::endl;
  }

  if (nzleptons>2) {
    std::cout << "Found more than two Z leptons : " << nzleptons
	      << std::endl;
  }

  if (nwleptons>1 || nzleptons>2) {
    //    std::cout << "PROBLEM: # W leptons: " << nwleptons
    //	      << "\t # Z leptons : " << nzleptons
    //	      << std::endl;
    *weirdEventsList << wzevt->run << ":" 
		     << wzevt->lumi << ":"<< wzevt->event 
		     << "\t # W l: " << nwleptons << "\t # Z l:" << nzleptons
		     << "\t mZ = " << wzevt->MZ
		     << std::endl;
    return; // Something wrong: IGNORE IT FOR NOW
  }

  // Look at W decay type
  if (iwlepton>=0) {
    int lid = abs(wzevt->genLeptons[iwlepton].Id());
    if (wzevt->genLeptons[iwlepton].ComesFromTau()) {
      if (lid == 11) { 
	  wDecaysByType[wtaue]++;
	} else if (lid == 13) { 
	  wDecaysByType[wtaumu]++;
	} else {
	  std::cout << "THIS SHOULD NOT BE HAPPENING: stable lepton neither e nor mu \n";
	}
      } else {
	if (lid == 11) { 
	  wDecaysByType[we]++;
	} else if (lid == 13) { 
	  wDecaysByType[wu]++;
	} else {
	  std::cout << "THIS SHOULD NOT BE HAPPENING: stable lepton neither e nor mu \n";
	}
      }
  } else {
    wDecaysByType[wtauh]++;
  }


  // Look at Z decay type

  ZDecayType decay=undefinedZDecay;
  if (nzleptons>=0) {

    int zlid1 = abs(wzevt->genLeptons[izleptons[0]].Id());

    bool isTauDecay1 = false;
    if (wzevt->genLeptons[izleptons[0]].ComesFromTau()) {
      isTauDecay1 = true;
    }
    if (!isTauDecay1) {
      if (zlid1 == 11) {
	decay = zee;
      } else if (zlid1 == 13) { 
	decay = zmumu;
      }
    } else {
      if (zlid1 == 11) {
	decay = ztteh;
      } else if (zlid1 == 13) { 
	decay = zttmuh;
      } else {
	decay = ztthh;
      }
    }

    if (nzleptons==2) {
      int zlid2 = abs(wzevt->genLeptons[izleptons[1]].Id());
      bool isTauDecay2 = wzevt->genLeptons[izleptons[1]].ComesFromTau();
      // Check that the 2 sides are consistent
      if (isTauDecay1 != isTauDecay2) {
	std::cout << "INCONSISTENT TAU DECAYS: " 
		  << isTauDecay1 << isTauDecay2  << "\t"
		  << wzevt->genLeptons[izleptons[0]].ComesFromTau()
		  << wzevt->genLeptons[izleptons[1]].ComesFromTau()
		  << std::endl;
	return;
      }
      if (!isTauDecay2) {
	if (zlid2 != zlid1) {
	  std::cout << "INCONSITENT Z LEPTON FLAVOURS: "
		    << zlid1 << "\t" << zlid2 
		    << "\t Tau decays: " 
		    << isTauDecay1 << isTauDecay2 << "\t"
		    << std::endl;
	  return;
	}
      } else { // Both tau decays: which type
	if (zlid1 == zlid2) {
	  if (zlid1 == 11) {
	    decay = zttee;
	  } else if (zlid1 == 13) {
	    decay = zttmumu;
	  }
	}  else {
	  decay = zttemu;
	}
      }
    }
  } else {
    decay = ztthh;
  }
  zDecaysByType[decay]++;


  // 
  // JET RESOLUTION ANALYSIS
  // 



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

    // Now find matching reco jet

    int matchedJet = -1;
    float drmax = 10.;
    for (int irj=0; irj<wzevt->recoJets.size(); irj++) {
      double dRgr = wzevt->genJets[i].DeltaR(wzevt->recoJets[irj]);
      if (dRgr < drmax) {
	matchedJet = irj;
	drmax = dRgr;
      }
    }
    if (matchedJet>=0 && drmax < 0.3) {
      _genJetPt = wzevt->genJets[i].Pt();
      _genJetPhi = wzevt->genJets[i].Phi();
      _genJetEta = wzevt->genJets[i].Eta();
      _recoJetPt = wzevt->recoJets[matchedJet].Pt();
      _recoJetPhi = wzevt->recoJets[matchedJet].Phi();
      _recoJetEta = wzevt->recoJets[matchedJet].Eta();
      _drRecoGenJet = drmax;

      jetResolutionTree->Fill();

    }

  }
}

void WZAnalysis::Finish(TFile * fout) {




  // PDG Values
  float BrWtoE    = 0.1075;
  float BrWtoMu   = 0.1057;
  float BrWtoTau  = 0.1125;
  float BrZToL    = 0.0337;

  float BrTauToE  = 0.1783;
  float BrTauToMu = 0.1741;


  float BrTauToH = 1.- BrTauToE - BrTauToMu;


  std::cout << " =========== W DECAYS ============= \n \n";

  int totWDecays = 0;
  for (int i=0; i<NRWDECAYTYPES; i++) {
    totWDecays += wDecaysByType[i];
  }
  // 
  for (int i=0; i<NRWDECAYTYPES; i++) {
    double fraction = (double) wDecaysByType[i] / (double) totWDecays;
    double error = sqrt(fraction*(1-fraction) / (double) totWDecays );
    std::cout << "W Type " << i << "\t events: " 
	      << wDecaysByType[i] << "\t BR = " << fraction 
	      << " +/- " << error
	      << std::endl;
  }
  std::cout << "Total: " << totWDecays
	    << "\t ot of " <<   totalNrEvents 
	    << " = " << 100.* (float) totWDecays / (float) totalNrEvents 
	    << " % " << std::endl << std::endl;




  std::cout << "EXPECTED USING PDG BR values: \n";
  std::cout << "W -> e       : " << BrWtoE / (BrWtoE + BrWtoMu + BrWtoTau) << std::endl;
  std::cout << "W -> mu      : " << BrWtoMu/ (BrWtoE + BrWtoMu + BrWtoTau) << std::endl;
  std::cout << "W -> tau ->e : " 
	    << BrWtoTau*BrTauToE / (BrWtoE + BrWtoMu + BrWtoTau) << std::endl;
  std::cout << "W -> tau ->mu : " 
	    << BrWtoTau*BrTauToMu / (BrWtoE + BrWtoMu + BrWtoTau)
	    << std::endl;
  std::cout << "W -> tau ->h  : " 
	    << BrWtoTau*BrTauToH / (BrWtoE + BrWtoMu + BrWtoTau)
	    << std::endl;



  std::cout << " =========== Z DECAYS ============= \n \n";

  int totZDecays = 0;
  for (int i=0; i<NRZDECAYTYPES; i++) {
    totZDecays += zDecaysByType[i];
  }
  // 
  for (int i=0; i<NRZDECAYTYPES; i++) {
    double fraction = (double) zDecaysByType[i] / (double) totZDecays;
    double error = sqrt(fraction*(1-fraction) / (double) totZDecays );
    std::cout << "Z Type " << i << "\t events: " 
	      << zDecaysByType[i] << "\t BR = " << fraction 
	      << " +/- " << error
	      << std::endl;
  }
  std::cout << "Total: " << totZDecays
	    << "\t ot of " <<   totalNrEvents 
	    << " = " << 100.* (float) totZDecays / (float) totalNrEvents 
	    << " % " << std::endl << std::endl;

  std::cout << "EXPECTED USING PDG BR values: \n";
  std::cout << "Z->ee / mumu / tautau : " 
	    << 1/3. << std::endl;
  std::cout << "Z->ttee               : " 
	    << BrTauToE*BrTauToE / 3 << std::endl;
  std::cout << "Z->ttmumu             : " 
	    << BrZToL*BrTauToMu*BrTauToMu / (3.*BrZToL) << std::endl;
  std::cout << "Z->ttemu               : " 
	    << 2.*BrTauToE*BrTauToMu / 3 << std::endl;
  std::cout << "Z->tteh               : " 
	    << 2.*BrTauToE*BrTauToH / 3 << std::endl;
  std::cout << "Z->ttmuh               : " 
	    << 2.*BrTauToMu*BrTauToH / 3 << std::endl;
  std::cout << "Z->tthh               : " 
	    << BrTauToH*BrTauToH / 3 << std::endl;


  //

  weirdEventsList->close();

  // 
  if (fout) {
    jetResolutionTree->Write();
  }

}
