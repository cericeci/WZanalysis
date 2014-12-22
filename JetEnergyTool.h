#ifndef JetEnergyTool_h
#define JetEnergyTool_h


#include "TH1D.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

// Define pure static class with histo definitions

class JES_map {

public:

  JES_map(std::string);

  float etaMin, etaMax;
  std::vector<float> pt_bins;
  std::vector<float> scales;

};


class JetEnergyTool {


public:
  static JetEnergyTool * GetInstance();
  void LoadEnergyScaleMap();
  float JetEnergyScale(float jetpt, float jeteta);
  float GetJetEnergyScale(float jetpt, float jeteta);

  void SetJESFile(string s);

  void FillMaps();

private:

  JetEnergyTool();

  static JetEnergyTool * _instance;

  std::string jes_fileName;


  std::vector<JES_map *> jesmaps;

}; 


#endif
