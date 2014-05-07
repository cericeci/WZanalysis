#ifndef UnfoldingHistoFactory_h
#define UnfoldingHistoFactory_h


#include "TH1D.h"

#include <string>


// Define pure static class with histo definitions

class UnfoldingHistogramFactory {
public:
  static TH1D * createZPtHistogram(std::string s, std::string title);
  static TH1D * createLeadingJetHistogram(std::string s, std::string title);

}; 


#endif
