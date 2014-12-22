#include "JetEnergyTool.h"



#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

// Global static pointer used to ensure a single instance of the class.

JetEnergyTool * JetEnergyTool::_instance = NULL; 



JES_map::JES_map(std::string line) {

  pt_bins.clear();
  scales.clear();

  std::istringstream iss(line);
  string word;
  double x,etamin,etamax,pt,s;
  if (!(iss >> etamin >> etamax >> word)) { return; } // error

  etaMin = etamin;
  etaMax = etamax;
  
  while ( iss >> pt >> s >> x) {
    pt_bins.push_back(pt);
    scales.push_back(s);	

  }

}

JetEnergyTool::JetEnergyTool() {

  jes_fileName = "";

}

JetEnergyTool * JetEnergyTool::GetInstance(){

  if (! _instance) {  // Only allow one instance of class to be generated.

    _instance = new JetEnergyTool();

  }

  return _instance;

}


void JetEnergyTool::SetJESFile(string s) {

   jes_fileName = s;
   FillMaps();
}


float JetEnergyTool::JetEnergyScale(float jetpt, float jeteta) {

  if (jes_fileName == "") {
    std::cout << "ERROR: CANOT APPLY JES CORRECTIONS - NO FILE GIVEN \n";
    return -999;
  }

  ifstream infile(jes_fileName.c_str());

  if (! infile.is_open()) {
    std::cout << "ERROR FILE WITH JES CORRECTIONS UNREADABLE !!! \n";
    return -999;
  }
  float scale=-1;

  string line;
  while (std::getline(infile, line))
    {

      if (line[0] == '#'
	  || line[0] == '{' ) continue;

      std::istringstream iss(line);
      string word;
      double x,etamin,etamax,pt,s;
      if (!(iss >> etamin >> etamax >> word)) { break; } // error

      if (jeteta<etamin  || jeteta > etamax) continue;
      std::vector<float> pt_bins;
      std::vector<float> scales;
      while ( iss >> pt >> s >> x) {
	pt_bins.push_back(pt);
	scales.push_back(s);

      }
      for (int i=0; i<pt_bins.size(); i++) {
	if (jetpt>=pt_bins[i] ) scale = scales[i];
	if (jetpt<pt_bins[i] ) {
	  if (i==0) scale = scales[i];
	  break;
	}
      }
    }

  return scale;
}

void JetEnergyTool::FillMaps() {

  if (jes_fileName == "") {
    std::cout << "ERROR: CANOT APPLY JES CORRECTIONS - NO FILE GIVEN \n";
  }

  ifstream infile(jes_fileName.c_str());

  if (! infile.is_open()) {
    std::cout << "ERROR FILE WITH JES CORRECTIONS UNREADABLE !!! \n";
  }
  float scale=-1;

  string line;
  while (std::getline(infile, line))
    {

      if (line[0] == '#'
	  || line[0] == '{' ) continue;

      JES_map * jm = new JES_map(line);
      jesmaps.push_back(jm);

    }

}

float JetEnergyTool::GetJetEnergyScale(float jetpt, float jeteta) {

  float scale=-999;

  for (int im=0; im<jesmaps.size(); im++) {

    if (jeteta<jesmaps[im]->etaMin  
	|| jeteta > jesmaps[im]->etaMax) continue;    

    for (int ipt=0; ipt<jesmaps[im]->pt_bins.size(); ipt++) {
      if (jetpt>=jesmaps[im]->pt_bins[ipt] ) 
	scale = jesmaps[im]->scales[ipt];
	if (jetpt < jesmaps[im]->pt_bins[ipt] ) {
	  if (ipt==0) scale = jesmaps[im]->scales[ipt];
	  break;
	}
      }
  }

  return scale;


}
