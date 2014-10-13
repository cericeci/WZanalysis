#include "UnfoldingHistogramFactory.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

// Global static pointer used to ensure a single instance of the class.

UnfoldingHistogramFactory * UnfoldingHistogramFactory::_instance = NULL; 


UnfoldingHistogramFactory * UnfoldingHistogramFactory::GetInstance(){

  if (! _instance) {  // Only allow one instance of class to be generated.

    _instance = new UnfoldingHistogramFactory;

  }

  return _instance;

}


TH1D *  UnfoldingHistogramFactory::createHistogram(std::string key,
						   std::string hname,
						   std::string htitle) {

  TH1D * h;

  //  auto it = binLimitsMap.find(key);
  map<string,vector<double> >::iterator it = binLimitsMap.find(key);
  if (it == binLimitsMap.end())  {
    createDefaultBinning(key);
  }

  int nBins = binLimitsMap[key].size() - 1;
  double *  bins = new double[nBins+1];
  for (int i=0; i< binLimitsMap[key].size(); i++) {
    bins[i] = (binLimitsMap[key])[i];
  }

  h = new TH1D(hname.c_str(),htitle.c_str(), nBins, bins);

  delete [] bins;
      
  return h;


}

void UnfoldingHistogramFactory::createDefaultBinning(std::string key) {

  std::cout << "CREATE DEFAULT BINNING FOR : " << key << std::endl;

  std::vector<double> binLimits;
  if (key == "ZPt") {

    std::vector<double> binLimits;
    double value = 30.;
    double binSize = 40;
    while (value<500.) {
      binLimits.push_back(value);
      if (value>100) binSize = 25.;
      if (value>300) binSize = 40.;
      if (value>400) binSize = 50.;
      value += binSize;
    }
  } else if (key == "LeadingJetPt") {

    double value = 30.;
    double binSize = 40;
    while (value<500.) {
      binLimits.push_back(value);
      if (value>150) binSize = 50.;
      if (value>200) binSize = 60.;
      if (value>300) binSize = 80.;
      if (value>400) binSize = 120.;
      value += binSize;
    }
  }

  binLimitsMap[key] = binLimits;

}



TH1D * UnfoldingHistogramFactory::createLeadingJetHistogram(std::string hname, 
							    std::string htitle) {

  UnfoldingHistogramFactory * instance = UnfoldingHistogramFactory::GetInstance();

  return instance->createHistogram("LeadingJetPt",hname,htitle);

}




TH1D *  UnfoldingHistogramFactory::createZPtHistogram(std::string hname, 
						      std::string htitle) {


  UnfoldingHistogramFactory * instance = UnfoldingHistogramFactory::GetInstance();

  return instance->createHistogram("ZPt",hname,htitle);

}


void UnfoldingHistogramFactory::SetBinning(string file) {


  std::cout << "SET BINNING BY HAND " << std::endl;

  ifstream infile(file.c_str());

  string line;
  while (std::getline(infile, line))
    {
      std::istringstream iss(line);
      string key;
      double x;
      if (!(iss >> key)) { break; } // error

      vector<double> limits;
      while (iss >> x) {
	limits.push_back(x);
      }
      if (limits.size()>1) {
	binLimitsMap.insert(pair<string,vector<double> >(key,limits));
      }
    }
}
  



