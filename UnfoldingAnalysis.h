#ifndef UnfoldingAnalysis_h
#define UnfoldingAnalysis_h

#include "RooUnfoldResponse.h"

#include "WZEvent.h"

#include "TFile.h"
#include "TH1D.h"


#include <string>



class UnfoldingAnalysis {


public:

  UnfoldingAnalysis(std::string k, WZEvent * e);

  virtual void Init() {};

  void CreateBaseHistos();

  void FillEvent(bool controlSample=false);

  virtual void EventAnalysis(bool controlSample=false)=0;

  virtual void Finish(TFile * fout=0);

  void FillPurityStability();

protected:

  virtual TH1D * createHistogram(std::string s, std::string title)=0;
  
  WZEvent * wzevt;

  TH1D * genHistos[5];
  TH1D * recoHistos[5];

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

  double leadingRecoJetPt;
  double leadingGenJetPt;

  TH1D * createHistogram(std::string s, std::string title);

  TH1D * hnGenJets[5];
  TH1D * hnRecoJets[5];

};


class  UnfoldingZPt : public UnfoldingAnalysis {

public:
  
  UnfoldingZPt(WZEvent * e);

  void EventAnalysis(bool controlSample=false);
  //  void Finish(TFile * fout=0);

  //  void Init();

protected:

  double recoZPt;
  double genZPt;

  TH1D * createHistogram(std::string s, std::string title);

};


#endif
