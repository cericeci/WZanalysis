#include "UnfoldingHistogramFactory.h"

#include <vector>


TH1D * UnfoldingHistogramFactory::createLeadingJetHistogram(std::string s, std::string title) {
  TH1D * h;

  if (true) {
    h = new TH1D(s.c_str(),title.c_str(), 10,30., 530.);

  } else { // Variable Bin Size

    std::vector<double> binLimits;
    double value = 1.;
    double binSize = 20;
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
    h = new TH1D(s.c_str(),title.c_str(), nBins, bins);

  }

  return h;

}


TH1D *  UnfoldingHistogramFactory::createZPtHistogram(std::string s, std::string title) {

  TH1D * h = new TH1D(s.c_str(),title.c_str(), 10,0., 500.);
  return h;

}
