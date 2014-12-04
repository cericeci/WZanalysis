#ifndef SystematicsManager_h
#define SystematicsManager_h


#include "TH1D.h"

#include <string>
#include <map>
#include <vector>

using namespace std;

// Define pure static class with histo definitions


class SystematicsManager {


public:
  static SystematicsManager * GetInstance();

  float GetValue(string key);
  void Setup(string filename);

private:

  SystematicsManager();
  void SetDefaultValues();

  static SystematicsManager * _instance;

  std::map<string,float> values;


}; 


#endif
