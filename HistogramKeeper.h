#ifndef HistogramKeeper_h
#define HistogramKeeper_h


#include "TH1D.h"
#include "TFile.h"

#include <string>
#include <map>
#include <vector>


// Singleton class


class HistogramKeeper {

public:

  static HistogramKeeper * GetInstance();

  void AddObject(TObject *, std::string);
  void WriteAll(TFile *);

private:

  HistogramKeeper();

  static HistogramKeeper * _instance;

  std::map<std::string,TObject *> objects;


};




#endif
