#ifndef WZAnalysis_h
#define WZAnalysis_h

#include "WZEvent.h"
#include <fstream>

#define NRWDECAYTYPES 6
#define NRZDECAYTYPES 9


enum WDecayType { undefinedWDecay,
		  we,
		  wu,
		  wtaue,
		  wtaumu,
		  wtauh };

enum ZDecayType { undefinedZDecay,
		  zee,
		  zmumu,
		  zttee,
		  zttmumu,
		  zttemu,
		  ztteh,
		  zttmuh,
		  ztthh };

class WZAnalysis {

public:

  WZAnalysis(WZEvent * e);

  virtual void Init();

  void CreateBaseHistos();

  virtual void EventAnalysis();

  virtual void Finish(TFile * fout =0);

protected:

  WZEvent * wzevt;

  int wDecaysByType[NRWDECAYTYPES];
  int zDecaysByType[NRZDECAYTYPES];

  int totalNrEvents;

  std::ofstream * weirdEventsList;

  // For Jet Resolution analysis

  float _genJetPt;
  float _genJetPhi;
  float _genJetEta;

  float _recoJetPt;
  float _recoJetPhi;
  float _recoJetEta;

  float _drRecoGenJet;

  TTree * jetResolutionTree;

  int   _channel;
  float _trigEff;

  TTree * eventSummaryTree;

};

#endif
