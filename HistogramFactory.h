#ifndef HistogramFactory_h
#define HistogramFactory_h

#include "TH1F.h"

#include <string>

class HistogramFactory {
public:
  static TH1F * createZmassBigHisto(std::string s, std::string title);
  static TH1F * createZmassHisto(std::string s, std::string title);
  static TH1F * createMETbigHisto(std::string s, std::string title);
  static TH1F * createMETHisto(std::string s, std::string title);
  static TH1F * createZptHisto(std::string s, std::string title);
  static TH1F * createWcandptHisto(std::string s, std::string title);
  static TH1F * createLeadingJetptHisto(std::string s, std::string title);
  static TH1F * createNjetsHisto(std::string s, std::string title);
  static TH1F * createDeltaPhi(std::string s, std::string title);
  static TH1F * createMTW(std::string s, std::string title);
  //
}; 


#endif
