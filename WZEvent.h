#ifndef WZEvent_h
#define WZEvent_h

#define WZBASECLASS WZGenEvent

#include "TLorentzVector.h"

#include <vector>


// #include "WZGenNew.h"
#include "WZGenEvent.h"

enum FinalState { undefined,
		  eee,
		  eem,
		  mme,
		  mmm};

enum PassedSelectionStep { failsSelection,
			   passesZSelection,
			   passesWSelection,
			   passesFullSelection};
		  


class RecoLepton : public TLorentzVector {

public:

  RecoLepton(double pt, double eta, double phi,
	     float ch, float id ) {
    SetPtEtaPhiM(pt, eta, phi, 0.);
    charge = ch;
    pdgid  = id;
  };

  float GetScaleFactor();

protected:
  float pdgid;
  float charge;

};


class GenS1Lepton : public TLorentzVector {

public:

  GenS1Lepton(double pt, double eta, double phi,
	      float id=0, float oVpid=-9999., 
	      float imTau = -9999. ) {

    mass = 0;
    pdgid  = id;    
    if (abs(id)==11) {
      mass = 0.000511;
    } else if (abs(id) == 13) {
      mass = 0.1057;
    }
    SetPtEtaPhiM(pt, eta, phi, mass);

    motherBoson     = oVpid;
    isTauDescendent = imTau;

  };

  int Id() { return pdgid; };
  int MotherBoson() { return motherBoson; };
  int ComesFromTau() { return isTauDescendent; };
  

protected:
  int   pdgid;
  int   motherBoson;
  int   isTauDescendent;
  float mass;
};



class WZEvent : public WZBASECLASS
{
public:
  WZEvent(TTree *tree);
  bool passesSelection();

  void ReadEvent();

  void PrintSummary();

  bool PassesGenCuts();


  float LeptonPt(int i);
  float LeptonBDT(int i);
  float LeptonCharge(int i);
  float LeptonEta(int i);
  float LeptonPhi(int i);
  int   LeptonPass2012ICHEP(int i);
  float LeptonIso(int i);
  float LeptonIsomva(int i);


  FinalState GetFinalState(){return final_state;}
  PassedSelectionStep GetSelectionLevel() { return selection_level;}
  

//     float pt[leptonNumber]={cWZ->pt1, cWZ->pt2, cWZ->pt3, cWZ->pt4};
//     float bdt[leptonNumber]={cWZ->bdt1, cWZ->bdt2, cWZ->bdt3, cWZ->bdt4};
//     float ch[leptonNumber]={cWZ->ch1, cWZ->ch2, cWZ->ch3, cWZ->ch4};
//     float eta[leptonNumber]={cWZ->eta1, cWZ->eta2, cWZ->eta3, cWZ->eta4};
//     float phi[leptonNumber]={cWZ->phi1, cWZ->phi2, cWZ->phi3, cWZ->phi4};
//     int pass2012ICHEP[leptonNumber]={cWZ->pass2012ICHEP1, cWZ->pass2012ICHEP2, cWZ->pass2012ICHEP3, cWZ->pass2012ICHEP4};
//     float iso[leptonNumber]={cWZ->iso1, cWZ->iso2, cWZ->iso3, cWZ->iso4};
//     float isomva[leptonNumber]={cWZ->isomva1, cWZ->isomva2, cWZ->isomva3, cWZ->isomva4};


//  std::vector<TLorentzVector> genLeptons;
  std::vector<GenS1Lepton> genLeptons;
  std::vector<int> genLeptonsIds;
  std::vector<TLorentzVector> genJets;
  std::vector<TLorentzVector> recoJets;

  std::vector<RecoLepton>  leptons;

  float NumZ() { return numZ;}
  float NumW() { return numW;}

  double SelectedZPt() {
    return selectedZPt;
  }


protected:

  std::vector<float*> pt;
  std::vector<float*> bdt;
  std::vector<float*> pdgid;
  std::vector<float*> ch;
  std::vector<float*> eta;
  std::vector<float*> phi;
  std::vector<int*> pass2012ICHEP;
  std::vector<float*> iso;
  std::vector<float*> isomva;



  FinalState          final_state;
  PassedSelectionStep selection_level;
  
  double selectedZPt;


  float numZ;
  float numW;
  float numMET;
  float num3e;
  float num2e1mu;
  float num1e2mu;
  float num3mu;
  float numMET3e;
  float numMET2e1mu;
  float numMET1e2mu;
  float numMET3mu;


};


#endif // #ifdef WZ_cxx
