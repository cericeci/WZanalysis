#ifndef UnfoldingAnalysis_h
#define UnfoldingAnalysis_h

#include "RooUnfoldResponse.h"

#include "WZEvent.h"

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"


#include <string>
#include <fstream>


class UnfoldingAnalysis {


public:

  UnfoldingAnalysis(std::string k, WZEvent * e);

  virtual void Init() {};

  void CreateBaseHistos();

  void FillEvent(bool controlSample=false);

  virtual void EventAnalysis(bool controlSample=false)=0;

  virtual void Finish(TFile * fout=0);

  void FillPurityStability();
  void ApplyLuminosityNormalization(double norm);


  void UseModifiedShape(bool use=true);

protected:

  virtual TH1D * createHistogram(std::string s, std::string title)=0;

  double GetRecoWeight();
  double GetGenWeight();

  virtual double GetShapeWeight(double var); // { return 1.;}

  
  WZEvent * wzevt;

  TH1D * genHistos[5];
  TH1D * recoHistos[5];

  TH1D * genXSHistos[5];

  TH1D * controlRecoHistos[5];
  TH1D * controlGenHistos[5];

  //  TH1D * stabilityHistos[4];
  //  TH1D * purityHistos[4];

  TH1D * purityPlot[5];
  TH1D * purityPlotDenominator[5];

  TH1D * stabilityPlot[5];
  TH1D * stabilityPlotDenominator[5];

  RooUnfoldResponse *response[5];

  std::string key;

  double * trueValue;
  double * recoValue;

  bool    normalizeToLumi;
  bool    useNormalizedWeights;
  bool    useModifiedShape;

  TTree * resolutionTree;

  double * shapeWeights;

  double missYield;
  double fakeYield;
  double fillYield;



};

class  UnfoldingLeadingJetPt : public UnfoldingAnalysis {

public:
  
  UnfoldingLeadingJetPt(WZEvent * e);
//     UnfoldingAnalysis("LeadJetPt", e)
    
//   {
//     trueValue = &leadingGenJetPt;
//     recoValue = &leadingRecoJetPt;
//   };

  void EventAnalysis(bool controlSample=false);


  void Finish(TFile * fout=0);

  void Init();

protected:

  //  double GetShapeWeight(double var);

  double leadingRecoJetPt;
  double leadingGenJetPt;

  double leadingRecoJetPhi;
  double leadingGenJetPhi;;

  double leadingRecoJetEta;
  double leadingGenJetEta;

  double leadingRecoJetDRZl;
  double leadingRecoJetDRWl;

  TH1D * createHistogram(std::string s, std::string title);

  TH1D * hnGenJets[5];
  TH1D * hnRecoJets[5];
  

};


class  UnfoldingNjets : public UnfoldingAnalysis {

public:
  
  UnfoldingNjets(WZEvent * e);

  void EventAnalysis(bool controlSample=false);


  // void Finish(TFile * fout=0);

  void Init();

protected:
  /*
  double leadingRecoJetPt;
  double leadingGenJetPt;

  double leadingRecoJetPhi;
  double leadingGenJetPhi;;

  double leadingRecoJetEta;
  double leadingGenJetEta;

  double leadingRecoJetDRZl;
  double leadingRecoJetDRWl;
  */
  double nGenJets;
  double nRecoJets;

  TH1D * createHistogram(std::string s, std::string title);

  TH1D * hnGenJets[5];
  TH1D * hnRecoJets[5];
  
};


class  UnfoldingZPt : public UnfoldingAnalysis {

public:
  
  UnfoldingZPt(WZEvent * e);

  void EventAnalysis(bool controlSample=false);
  //  void Finish(TFile * fout=0);

  void Init();

protected:
  
  double recoZPt;
  double genZPt;
  double recoZPhi;
  double genZPhi;
  double recoZEta;
  double genZEta;


  TH1D * createHistogram(std::string s, std::string title);

};


#endif
