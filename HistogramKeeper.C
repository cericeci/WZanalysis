#include "HistogramKeeper.h"


#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;

// Global static pointer used to ensure a single instance of the class.

HistogramKeeper * HistogramKeeper::_instance = NULL; 


HistogramKeeper::HistogramKeeper() {

}



HistogramKeeper * HistogramKeeper::GetInstance(){

  if (! _instance) {  // Only allow one instance of class to be generated.

    _instance = new HistogramKeeper();

  }

  return _instance;

}


void HistogramKeeper::AddObject(TObject * obj, std::string key) {

  map<string, TObject *  >::iterator it = objects.find(key);
  if (it == objects.end())  {
    objects.insert(pair<string, TObject *>(key, obj));
  } else {
    std::cout << "ERROR-HistogramKeeper: adding object with already used key : "
	      << key << endl;
  }
}


void HistogramKeeper::WriteAll(TFile * fout){

  if (!fout) return;

  fout->cd();

  map<string, TObject * >::iterator it;

  for (it=objects.begin(); it!= objects.end(); it++) {
    it->second->Write(it->first.c_str());
  }




}
