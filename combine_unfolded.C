#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>

#include "UnfoldingHistogramFactory.h"

#define NR_BINS     10
#define NR_CHANNELS 4



// Gives index as a function of bin and channel
// All indices assume to start from zero


int index_from_bin_and_channel(int ibin, int ichan) {

  if ( ibin<0 || ibin>=NR_BINS 
       || ichan < 0 || ichan > NR_CHANNELS) {

    std::cout << "index_vs_bin_and_channel : INVALID INPUT " << ibin << "\t" << ichan << std::endl;
    return -1;
  }
  return ibin*NR_CHANNELS + ichan;
}

void bin_and_channel_from_index(int index, int & bin, int & channel) {

  if (index<0 || index > NR_BINS * NR_CHANNELS) {
    std::cout << "bin_and_channel_from_index: INVALID INPUT : " << index << std::endl;
    bin     = -1;
    channel = -1;
    return;
  }
  bin      = index/NR_CHANNELS;
  channel  = index % NR_CHANNELS;

}

using namespace std;

int main(int argc, char **argv)
{
  //definirati file-ove iz kojih citam i u koje pisem
  char * variableName(0);
  char * inputFile(0);
  bool gotVarName  = false;
  char c;

  while ((c = getopt (argc, argv, "v:i:")) != -1)
    switch (c)
	{
	case 'v':
	  gotVarName = true;
	  variableName = new char[strlen(optarg)+1];
	  strcpy(variableName,optarg);
	  break;
	case 'i':
	  inputFile = new char[strlen(optarg)+1];
	  strcpy(inputFile,optarg);
	  break;
	default:
	  std::cout << "usage: -r responseFile [-d <dataFile>]   \n";
	  abort ();
	}
  


  string variable = variableName;

  // Read matrices from input


  TFile * fin= new TFile(inputFile,"READ");

  fin->ls();

  int nr_bins     = 9;
  int nr_channels = 4;

  int dimension = nr_bins*nr_channels;

  TVectorD * measuredValues[4];


  vector<TMatrixD*>  unfoldingCovariance;   // [nr_channels]
  vector<TMatrixD*>  channelCovariance;     // [nr_bins]

  for (int ibin=1; ibin<=nr_bins; ibin++) {

    std::ostringstream matName;
    matName << "covMatrix_" << variable << "_bin" << ibin;

    TObject * matrix =  fin->Get(matName.str().c_str());
    std::cout << "matrix pointer: " << matrix << endl;
    if (!matrix) {
      std::cout << "Matrix missing in input file, cannot proceed : " << matName.str() << std::endl;
      return -1;
    }
    channelCovariance.push_back((TMatrixD*) matrix->Clone(matName.str().c_str()));

  }

  for (int ich=0; ich<nr_channels; ich++) {

    std::ostringstream matName;
    matName << "covarianceMatrix_unfolding_" << variable << "_ch" << ich;
    TMatrixD * matrix = (TMatrixD*) fin->Get(matName.str().c_str());
    if (!matrix) {
      std::cout << "Matrix missing in input file, cannot proceed : " << matName.str() << std::endl;
      return -1;
    }
    unfoldingCovariance.push_back((TMatrixD*) matrix->Clone(matName.str().c_str()));    
  }








  TMatrixD  fullCovariance(dimension,dimension);

  // Build full covariance matrix

  for (int ibin1=0; ibin1<nr_bins; ibin1++) {
    for (int ibin2=0; ibin2<nr_bins; ibin2++) {

      for (int ichan1=0; ichan1<nr_channels; ichan1++) {
	for (int ichan2=0; ichan2<nr_channels; ichan2++) {

	  double cov;
	  if (ibin1 == ibin2 && ichan1 == ichan2) {
	    std::cout << "WHAT DO WE DO HERE??? \n";
	    cov = (*channelCovariance[ichan1])(ibin1,ibin2);
	  }
	  else if (ibin1 == ibin2) {
	    cov = (*channelCovariance[ibin1])(ichan1,ichan2);
	  } else if (ichan1 == ichan2) {
	    cov = (*channelCovariance[ichan1])(ibin1,ibin2);
	  } else {
	    cov = 0.;
	  }
	  int index1 = index_from_bin_and_channel(ibin1,ichan1);
	  int index2 = index_from_bin_and_channel(ibin2,ichan2);
	  fullCovariance(index1,index2) = cov;
	}  
      }    
    }
  }

  // Build full measurements vector

  int nall = nr_bins*nr_channels;

  TVectorD  a(nall);
  for (int ibin=0; ibin<nr_bins; ibin++) {
    for (int ichan=0; ichan<nr_channels; ichan++) {
      double val = (*measuredValues[ichan])(ibin);
      int index = index_from_bin_and_channel(ibin,ichan);
      a(index) = val;
    }
  }


  // Invert full covariance matrix and auxiliary vector
  // 
  // M = Cov-1
  // 
  // f = M*a

  TMatrixD M = fullCovariance.Invert();
  TVectorD f = M*a;


  // Now builds Matrix D: nbins*nbins


  TMatrix D(nr_bins,nr_bins);
  for (int ibin1=0; ibin1<nr_bins; ibin1++) {
    for (int ibin2=0; ibin2<nr_bins; ibin2++) {

      double val=0;
      for (int ichan1=0; ichan1<nr_channels; ichan1++) {
	for (int ichan2=0; ichan2<nr_channels; ichan2++) {
	  int index1 = index_from_bin_and_channel(ibin1,ichan1);
	  int index2 = index_from_bin_and_channel(ibin2,ichan2);
	  val += M(index1,index2);

	}
      }
    }
  }


  // Build vector of n_bins component, each summing all channels of that bin

  TVectorD g(nr_bins);
  for (int ibin1=0; ibin1<nr_bins; ibin1++) {
    double val =0;
    for (int ichan1=0; ichan1<nr_channels; ichan1++) {
      int index1 = index_from_bin_and_channel(ibin1,ichan1);      
      val += f(index1);
    }
    g(ibin1) = val;
  }


  // And now get result and final covariance matrix

  TVectorD result(nr_bins);
  TMatrixD finalCovariance(nr_bins,nr_bins);

  TMatrixD Dinv=D.Invert();
  
  result = Dinv * g;


  finalCovariance = Dinv;
 
}
