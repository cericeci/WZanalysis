#include "SystematicsManager.h"



#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

// Global static pointer used to ensure a single instance of the class.

SystematicsManager * SystematicsManager::_instance = NULL; 


SystematicsManager::SystematicsManager() {

  SetDefaultValues();

}

void SystematicsManager::SetDefaultValues() {

  values["ele_SF"]     = 0.;
  values["mu_SF"]      = 0.;
  values["ele_scale"]  = 0.;
  values["mu_scale"]   = 0.;
  values["PU"]         = 0.;
  values["JER"]        = 0.;
  values["JES"]        = 0.;
  

}


void SystematicsManager::Setup(string filename) {

  std::cout << "Systematics SETUP \n";

  ifstream infile(filename.c_str());


  string line;
  while (std::getline(infile, line))
    {

      std::cout << "\t Processing line: " << line << std::endl;
      if (line[0] == '#')  continue;

      std::istringstream iss(line);
      string key;
      float value;

      if ((iss >> key >> value)) { 

	map<string,float >::iterator it = values.find(key);
	if (it == values.end())  {
	  std::cout << "DEFINING UNKNOWN SYSTEMATICS: " 
		    << key << std::endl;
	  continue;
	}
	values[key] = value;
	std::cout << "DEF SYST: " 
		  << key << "\t:\t" << value << std::endl;
      }
    }
}

SystematicsManager * SystematicsManager::GetInstance(){

  if (! _instance) {  // Only allow one instance of class to be generated.

    _instance = new SystematicsManager();

  }

  return _instance;

}


float SystematicsManager::GetValue(string key) {

  map<string,float >::iterator it = values.find(key);
  if (it == values.end())  {
    std::cout << "UNKNOWN SYSTEMATICS: " << key << std::endl;
    return -999;
  }
  return values[key];

}

