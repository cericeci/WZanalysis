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

void Combine_unfolded(
		      std::map<int, TH1D* >      & h_dsigma,
		      std::map<int, TH1D* >      & h_dsigmadx,
		      std::map<int, TMatrixD *>  & unfoldingCovariance,
		      std::map<int, TMatrixD *>  & binCovariance
		      )
{


  // Read measurements from input

  
  std::map<int, TH1D* >::iterator ithist, ithist2;

  std::map<int, TH1D* >   differential_xsections;

  for (ithist=h_dsigma.begin(); ithist!=h_dsigma.end(); ithist++) {

    int channel = ithist->first;

    std::cout << "CHANNEL = " << channel << std::endl
	      << "============== \n";

    TH1D * h = ithist->second;
    int nbins = h->GetNbinsX();

    // ithist2=h_dsigmadx.find(ithist->first);
    // TH1D * h2(0);
    // if (ithist2 != h_dsigmadx.end()) {
    //   int nbins2 = ithist2->second->GetNbinsX();
    //   if (nbins2 != nbins) {
    // 	std::cout << "MISMATCH IN # of bins!!! \n";
    //   } else {
    // 	h2 = ithist2->second;
    //   }
    // }
       
    std::ostringstream diffXsName;
    diffXsName << "h_diffXS_" << channel;
    TH1D * h_diffxs = (TH1D*) h->Clone(diffXsName.str().c_str());

    for (int ibin=0; ibin<=nbins+1; ibin++) {
      double val = h->GetBinContent(ibin);
      double wid = h->GetBinWidth(ibin);
      double err = h->GetBinError(ibin);

      h_diffxs->SetBinContent(ibin,val/wid);
      h_diffxs->SetBinError(ibin,err/wid);

      std::cout << "Bin : " << ibin
		<< "  Val = " << val << "  Error = " << err 
		<< "  Val/dX = " << val/wid
		<< "  Rel. =  " << err/val ;
      // if (h2) {
      // 	double val2 = h2->GetBinContent(ibin);
      // 	double err2 = h2->GetBinError(ibin);
      // 	std::cout << "\t DSDX  Val = " << val2 << "  Error = " << err2 
      // 		  << "  Rel. =  " << err2/val2 
      // 		  << "\t R = " << (err/val) / (err2/val2);
      // }
      std::cout << std::endl;

    }
    differential_xsections.insert(std::pair<int,TH1D *> (channel, h_diffxs));
  }

  // Read matrices from input

  int nr_bins     = binCovariance.size();
  int nr_channels = unfoldingCovariance.size();

  std::map<int, TMatrixD *>::iterator it;

  std::cout << "FULL BLUE UNFOLDING " 
	    << "\t # bins = " << nr_bins
	    << "\t # channels = " << nr_channels << endl;


  // BUILD CORRELATION AND UNFOLDING COVARINACE MATRIX
  ///   DOUBT: 
  /// we have nbins measured bins
  /// BUT: unfolding covariance matrix is (nbins+2)*(nbins+2)
  // Contains under- and overflow bins: for bin=0 and bin=NBins
  //   OUR CHOICE:

  // NOTE NOTE NOTE:
  // We DO *NOT* TAKE UNDER- AND OVERFLOW BINS

  std::map<int, TMatrixD *>  unfoldingCorrelation;
  std::map<int, TMatrixD *>  unfoldingErrorMatrix;

  for (it=unfoldingCovariance.begin(); it!=unfoldingCovariance.end(); it++) {
    std::cout << "Channel " << it->first << std::endl;
    int channel = it->first;
    //    it->second->Print();
    int dim=it->second->GetNrows();
    // remove under- and overflow bins
    dim -= 2;

    TMatrixD * correlationMatrix    = new TMatrixD(dim,dim);
    TMatrixD * errorMatrix          = new TMatrixD(dim,dim);
    for (int irow=0; irow<correlationMatrix->GetNrows(); irow++) {
      for (int icol=irow; icol<correlationMatrix->GetNcols(); icol++) {
	double cov   = 	(*it->second)(irow+1,icol+1);
	double sig_x = sqrt( (*it->second)(irow+1,irow+1)); 
	double sig_y = sqrt( (*it->second)(icol+1,icol+1));
	double corr_xy = cov/(sig_x*sig_y);

	//	std::cout << "row : " << irow << "  col : " << icol
	//		  << " cov = " << cov << "\t sx = " << sig_x << "\t xy = " << sig_y << std::endl;
	(*correlationMatrix)(irow,icol) = corr_xy;
	(*correlationMatrix)(icol,irow) = corr_xy;
	//  And now produce real unfolding covariance
	// Correctly normalized to cross section

	
	double unf_error_x = differential_xsections[channel]->GetBinError(irow+1);
	double unf_error_y = differential_xsections[channel]->GetBinError(icol+1);
	double cov_xy   = corr_xy*unf_error_x*unf_error_y;

	(*errorMatrix)(irow,icol) = cov_xy;
	(*errorMatrix)(icol,irow) = cov_xy;
      }
    }

    unfoldingCorrelation.insert(std::pair<int,TMatrixD*> ( it->first, correlationMatrix));
    unfoldingErrorMatrix.insert(std::pair<int,TMatrixD*> ( it->first, errorMatrix));
    std::cout << "Correlation Matrix \n";
    correlationMatrix->Print();
    std::cout << "Error Matrix \n";
    errorMatrix->Print();

  }

  //

  int dimension = nr_bins*nr_channels;

  //  vector<TMatrixD*>  unfoldingCovariance;   // [nr_channels]
  //  vector<TMatrixD*>  channelCovariance;     // [nr_bins]

  TMatrixD  fullCovariance(dimension,dimension);

  // Build full covariance matrix

  for (int ibin1=0; ibin1<nr_bins; ibin1++) {
    for (int ibin2=0; ibin2<nr_bins; ibin2++) {

      // SHIFT TO JUMP OVER UNDERFLOW BIN
      int binIndex = ibin1+1;

      for (int ichan1=0; ichan1<nr_channels; ichan1++) {
	for (int ichan2=0; ichan2<nr_channels; ichan2++) {

	  double cov(0);
	  if (ibin1 == ibin2 && ichan1 == ichan2) {
	    cov = (*binCovariance[binIndex])(ichan1,ichan2);
	  }
	  else if (ibin1 == ibin2) {
	    cov = (*binCovariance[binIndex])(ichan1,ichan2);

	  } else if (ichan1 == ichan2) {
	    cov = (*unfoldingErrorMatrix[ichan1])(ibin1,ibin2);
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

  // std::cout << "FULL BLOW COVARIANCE : \n";
  // std::cout << "====================== \n \n";
  // fullCovariance.Print();
  std::cout << "Full Covariance Determinant = " << fullCovariance.Determinant() << std::endl;

  // Build full measurements vector

  //  int nall = nr_bins*nr_channels;

  TVectorD  a(dimension);
  for (int ibin=0; ibin<nr_bins; ibin++) {
    for (int ichan=0; ichan<nr_channels; ichan++) {
      double val = differential_xsections[ichan]->GetBinContent(ibin+1);
      int index = index_from_bin_and_channel(ibin,ichan);
      a(index) = val;
    }
  }

  std::cout << "MEASUREMENT VECTOR: \n";
  std::cout << "====================== \n \n";
  //  a.Print();

  // Invert full covariance matrix and auxiliary vector
  // 
  // M = Cov-1
  // 
  // f = M*a

  TMatrixD M = fullCovariance.Invert();
  TVectorD f = M*a;
  //  TVectorD f(a);
  //  f *= M;

  std::cout << "First inversion done \n";

  //  M.Print();

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
      D(ibin1,ibin2) = val;
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

  //  D.Print();

  TVectorD result(nr_bins);
  TMatrixD finalCovariance(nr_bins,nr_bins);

  TMatrixD Dinv=D.Invert();


  std::cout << "Second inversion done \n";
  
  result = Dinv * g;

  finalCovariance = Dinv;

  Dinv.Print();

  std::cout << "FINAL RESULTS: \n \n";
  for (int ibin=0; ibin<nr_bins ;ibin++) {
    std::cout << "Bin " << ibin 
	      << "\t Value = " << result(ibin)
      //	      << " +/- " << sqrt(finalCovariance(ibin,ibin)) << std::endl;
	      << " +/- " << sqrt(D(ibin,ibin)) << std::endl;
  }
 
}
