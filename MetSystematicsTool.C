#include "MetSystematicsTool.h"

#include "metsys.h"

#include <iostream>



// Global static pointer used to ensure a single instance of the class.

MetSystematicsTool * MetSystematicsTool::_instance = NULL; 



MetSystematicsTool::MetSystematicsTool() {

  metType=-1;
  met_fileName="";


}

MetSystematicsTool * MetSystematicsTool::GetInstance(){

  if (! _instance) {  // Only allow one instance of class to be generated.

    _instance = new MetSystematicsTool();

  }

  return _instance;

}


void MetSystematicsTool::LoadMETMap(int type)
{

  METEventMap * newmap = new METEventMap();

  TChain metch("MET");
  metch.Add(met_fileName.c_str());

  TTree *met_tTree  =(TTree*)&metch;
  metsys *cMET      = new metsys(met_tTree);
  Int_t events= met_tTree->GetEntries();

  //  TLeaf * leaf = met_tTree->GetLeaf("patPFMet");
  //  float * pleaf;
  //  leaf->SetAddress(pleaf);

  //  std::cout << "Leaf address: " << leaf << std::endl;

  for  (Int_t k = 0; k<events; k++) {

    met_tTree->GetEntry(k);

    float metPt = -1;
    float metPhi = -999;
    if (type == 1) {
      metPt  = cMET->patType1CorrectedPFMet;
      metPhi = cMET->patType1CorrectedPFMetPhi;
    } else if (type == 2) {
      metPt  =    cMET->patType1CorrectedPFMetElectronEnDown;
      metPhi =    cMET->patType1CorrectedPFMetElectronEnDownPhi;
    } else if (type == 3) {
      metPt  =    cMET->patType1CorrectedPFMetElectronEnUp;
      metPhi =    cMET->patType1CorrectedPFMetElectronEnUpPhi;
    } else if (type == 4) {
      metPt  =    cMET->patType1CorrectedPFMetJetEnDown;
      metPhi =    cMET->patType1CorrectedPFMetJetEnDownPhi;
    } else if (type == 5) {
      metPt =    cMET->patType1CorrectedPFMetJetEnUp;
      metPhi =    cMET->patType1CorrectedPFMetJetEnUpPhi;
    } else if (type == 6) {
      metPt  =    cMET->patType1CorrectedPFMetJetResDown;
      metPhi =    cMET->patType1CorrectedPFMetJetResDownPhi;
    } else if (type == 7) {
      metPt  =    cMET->patType1CorrectedPFMetJetResUp;
      metPhi =    cMET->patType1CorrectedPFMetJetResUpPhi;
    } else if (type == 8) {
      metPt  =    cMET->patType1CorrectedPFMetMuonEnDown;
      metPhi =    cMET->patType1CorrectedPFMetMuonEnDownPhi;
    } else if (type == 9) {
      metPt  =    cMET->patType1CorrectedPFMetMuonEnUp;
      metPhi =    cMET->patType1CorrectedPFMetMuonEnUpPhi;
    } else if (type == 10) {
      metPt  =    cMET->patType1CorrectedPFMetTauEnDown;
      metPhi =    cMET->patType1CorrectedPFMetTauEnDownPhi;
    } else if (type == 11) {
      metPt  =    cMET->patType1CorrectedPFMetTauEnUp;
      metPhi =    cMET->patType1CorrectedPFMetTauEnUpPhi;
    } else if (type == 12) {
      metPt  =    cMET->patType1CorrectedPFMetUnclusteredEnDown;
      metPhi =    cMET->patType1CorrectedPFMetUnclusteredEnDownPhi;
    } else if (type == 13) {
      metPt  =    cMET->patType1CorrectedPFMetUnclusteredEnUp;
      metPhi =    cMET->patType1CorrectedPFMetUnclusteredEnUpPhi;
    } else if (type == 20) {
      metPt  = cMET->patPFMet;
      metPhi = cMET->patPFMetPhi;
    } 


    //    metValues.insert(make_pair(make_pair(cMET->run,cMET->event), 
    //			       metPt));

    TLorentzVector metp4;
    metp4.SetPtEtaPhiM(metPt,0,metPhi,0);

    newmap->insert(make_pair(make_pair(cMET->run,cMET->event), 
			       metp4));


  }

  metMaps.insert(make_pair(type,newmap));

}

void MetSystematicsTool::SetInputFile(string s)
{
  met_fileName = s;
}


float MetSystematicsTool::GetMETValue(int run, int event, int type){

  TLorentzVector pmet = GetMETVector(run, event, type);
  return pmet.Pt();

}

TLorentzVector MetSystematicsTool::GetMETVector(int run, int event, int type)
{


  TLorentzVector pmet;

  //  std::cout << "Entering MET MAP: " << type
  //    << "\t # of maps : " << metMaps.size() << std::endl;

  // Load the map if it isn't loaded already

  std::map< int, METEventMap *>::iterator iter;

  iter = metMaps.find(type);

  if (iter==metMaps.end()) {
    LoadMETMap(type);
    iter = metMaps.find(type);
  }

  std::map< pair<int,int>, TLorentzVector>::iterator iter2;  

  iter2 = iter->second->find(make_pair(run,event));
  float met_pt;
  if (iter2 != iter->second->end() ) {
    pmet = iter2->second;
    //    met_pt = iter2->second.Pt();
  }
  return pmet;

}
