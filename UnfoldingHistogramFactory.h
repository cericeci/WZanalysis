#ifndef UnfoldingHistoFactory_h
#define UnfoldingHistoFactory_h


#include "TH1D.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

// Define pure static class with histo definitions


class UnfoldingHistogramFactory {


public:
  static UnfoldingHistogramFactory * GetInstance();
  void SetBinning(std::string file);

  TH1D * createHistogram(std::string key, std::string name, std::string title);
  static TH1D * createZPtHistogram(std::string name, std::string title);
  static TH1D * createLeadingJetHistogram(std::string name, std::string title);


private:

  UnfoldingHistogramFactory() {};
  UnfoldingHistogramFactory(UnfoldingHistogramFactory const & ) {};


  void createDefaultBinning(std::string key);

  static UnfoldingHistogramFactory * _instance;

  std::map<string,vector<double> > binLimitsMap;


}; 


#endif
