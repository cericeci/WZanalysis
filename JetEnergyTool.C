#include "JetEnergyTool.h"



#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

// Global static pointer used to ensure a single instance of the class.

JetEnergyTool * JetEnergyTool::_instance = NULL; 


JetEnergyTool::JetEnergyTool() {

  jes_fileName = "";

}

JetEnergyTool * JetEnergyTool::GetInstance(){

  if (! _instance) {  // Only allow one instance of class to be generated.

    _instance = new JetEnergyTool();

  }

  return _instance;

}


float JetEnergyTool::JetEnergyScale(float jetpt, float jeteta) {

  if (jes_fileName == "") {
    std::cout << "ERROR: CANOT APPLY JES CORRECTIONS - NO FILE GIVEN \n";
    return -999;
  }

  ifstream infile(jes_fileName.c_str());

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

