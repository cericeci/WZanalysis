#ifndef MetSystematicsTool_h
#define MetSystematicsTool_h


#include "TH1D.h"
#include "TLorentzVector.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

// Define pure static class with histo definitions


typedef std::map< pair<int,int>, TLorentzVector> METEventMap;

class MetSystematicsTool {


public:
  static MetSystematicsTool * GetInstance();

  //  void LoadMETValues();
  void LoadMETMap(int type);
  void SetInputFile(string s);
  float GetMETValue(int run, int event, int type);
  TLorentzVector GetMETVector(int run, int event, int type);


private:

  MetSystematicsTool();
  static MetSystematicsTool * _instance;

  int metType;
  std::string met_fileName;


  //  std::map< pair<int,int>, float> metValues;
  //  std::map< pair<int,int>, TLorentzVector> metVectors;
  //  METEventMap metVectors;

  std::map<int,METEventMap * > metMaps;

  

}; 


#endif
