#ifndef JetEnergyTool_h
#define JetEnergyTool_h


#include "TH1D.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

// Define pure static class with histo definitions


class JetEnergyTool {


public:
  static JetEnergyTool * GetInstance();
  void LoadEnergyScaleMap();
  float JetEnergyScale(float jetpt, float jeteta);

  void SetJESFile(string s) { jes_fileName = s;}

private:

  JetEnergyTool();

  static JetEnergyTool * _instance;

  std::string jes_fileName;


}; 


#endif
