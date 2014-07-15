#include "HistogramFactory.h"

#include <vector>

TH1F *  HistogramFactory::createZmassBigHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 100,60, 120);
  return h;

}

TH1F *  HistogramFactory::createZmassHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 40,70, 110);
  return h;

}

TH1F *  HistogramFactory::createMETbigHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 150,0., 300);
  return h;

}

TH1F *  HistogramFactory::createMETHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 150,0., 150);
  return h;

}

TH1F *  HistogramFactory::createZptHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 40,0., 400);
  return h;

}

TH1F *  HistogramFactory::createWcandptHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 40,0., 400);
  return h;

}

TH1F *  HistogramFactory::createLeadingJetptHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(), 40,0., 200);
  return h;

}

TH1F *  HistogramFactory::createNjetsHisto(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(),4, 0, 4);
  return h;

}

TH1F *  HistogramFactory::createDeltaPhi(std::string s, std::string title) {

  TH1F * h = new TH1F(s.c_str(),title.c_str(),32, 0, 3.2);
  return h;

}


