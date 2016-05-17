#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompLU.h"
#include "TLorentzVector.h"

#include "Math/SMatrix.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <set>


#define NR_BINS     10
#define NR_CHANNELS 4



// Gives index as a function of bin and channel
// All indices assume to start from zero

typedef  ROOT::Math::SMatrix< double, 16, 16, ROOT::Math::MatRepStd<double,16,16 > >  AlgebraicMatrix16x16;
typedef  ROOT::Math::SMatrix< double, 16, 16, ROOT::Math::MatRepSym<double,16> >      AlgebraicSymMatrix16x16;


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


TMatrixD *  GetFinalCovariance(int nbins, int nchannels, 
			       TMatrixD M ) {
  // M is the inverse of the full Covariance Matrix


  //  TMatrix M=fullCovariance.Invert();

  TMatrixD D(nbins,nbins);

  for (int ibin1=0; ibin1<nbins; ibin1++) {
    for (int ibin2=0; ibin2<nbins; ibin2++) {

      double val=0;
      for (int ichan1=0; ichan1<nchannels; ichan1++) {
	for (int ichan2=0; ichan2<nchannels; ichan2++) {
	  int index1 = index_from_bin_and_channel(ibin1,ichan1);
	  int index2 = index_from_bin_and_channel(ibin2,ichan2);
	  val += M(index1,index2);
	}
      }
      D(ibin1,ibin2) = val;
    }
  }

  TMatrixD Dcopy(D);

  TMatrixD * finalCovariance  = new TMatrixD(D.Invert());

  TMatrixD testUnity = Dcopy*(*finalCovariance);

  std::cout << "======= SANITY CHECK =============== \n";
  testUnity.Print();
  std::cout << "======= END SANITY CHECK =============== \n";

  return finalCovariance;



}

void DumpMatrixToFile(TMatrixD matrix,std::string name="matrix.dat") {

  std::ofstream outf(name.c_str(),std::ofstream::out);

  outf << matrix.GetNrows() << "\t" << matrix.GetNcols() << std::endl;
  for (int irow=0; irow<matrix.GetNrows(); irow++) {
    for (int icol=0; icol<matrix.GetNcols(); icol++) {
      outf << "\t" << matrix(irow,icol);
    }
    outf << endl;
  }
  outf.close();

}

void CheckMatrixInvertible(TMatrixD a)
{

  TDecompLU lu(a);
  TMatrixD b;
  if (!lu.Decompose()) {
    cout << "Decomposition failed, matrix singular ?" << endl;
    cout << "condition number = " << lu.GetCondition() << endl;
  } else {
    lu.Invert(b);
    cout << "Invertible \n";
  }

    for (int irow=1; irow < a.GetNrows(); irow++) {
      int irowref = 4*(irow/4);
      std::cout << "Row Ratio: " << irowref << "  /  " << irow << "   :   ";
      for (int icol=0; icol<a.GetNcols(); icol++) {
	std::cout << "\t" << a(irowref,icol)/a(irow,icol);
      }
      std::cout << endl;
    }


}

int main(int argc, char **argv)
{
  //definirati file-ove iz kojih citam i u koje pisem
  char * variableName(0);
  char * inputFile(0);
  bool gotVarName  = false;
  bool useFullCorrelations=false;
  bool showPlots = false;
  char c;

  while ((c = getopt (argc, argv, "v:i:Fp")) != -1)
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
	case 'F':
	  useFullCorrelations=true;
	  break;
	case 'p':
	  showPlots = true;
	  break;
	default:
	  std::cout << "usage: -r responseFile [-d <dataFile>]   \n";
	  abort ();
	}
  



  std::ostringstream outfilename;
  outfilename<<"unfoldingFinalResults/blue_combination_"<< variableName <<".root";
  TFile * fout = new TFile(outfilename.str().c_str(),"RECREATE");


  string variable = variableName;

  //void Combine_unfolded(
  std::map<int, TH1D* >       h_dsigma;
  std::map<int, TH1D* >       h_dsigmadx;
  std::map<int, TMatrixD *>   unfoldingCovariance;
  std::map<int, TMatrixD *>   binCovariance;

  TFile * fin= new TFile(inputFile,"READ");

  //  fin->ls();


  //   READ RESULTS AND ERROR MATRICES FROM INPUT
  //
  // 

  int nr_of_bins = -1;

  int ichan=0;

  for (int ichan=0; ichan<4; ichan++) {
    std::ostringstream name,nameCopy, nameDiff, nameDiffcopy;
    name      << "h_Dsigma_" << ichan;
    nameCopy  << "h_crossSection_inclusive" << ichan;
    nameDiff  << "h_crossSection_incl_diff_" << ichan;
    TH1D* h = (TH1D*) fin->Get(name.str().c_str());
    if (h) { 
      //      name << "_copy";
      TH1D* hc    = (TH1D*) h->Clone(nameCopy.str().c_str());
      TH1D* hdiff = (TH1D*) h->Clone(nameDiff.str().c_str());
      hc->SetDirectory(0);
      nr_of_bins = hc->GetNbinsX();
      h_dsigma.insert(pair<int,TH1D*>(ichan,hc));
      //      h_dsigmadx.insert(pair<int,TH1D*>(ichan,hdiff));
    } else {
      std::cout << "DID NOT FIND HISTO: " << name.str() << std::endl;
      //      break;
    }
    // 
    std::ostringstream unfMatrixName;
    unfMatrixName << "unfoldingCov_" << ichan;
    TMatrixD * matrix = (TMatrixD*) fin->Get(unfMatrixName.str().c_str());
    if (matrix) {
      TMatrixD * mc = new TMatrixD(*matrix);
      //      mc->SetDirectory(0);
      unfoldingCovariance.insert(pair<int, TMatrixD*>(ichan,mc));
    } else {
      std::cout << "DID NOT FIND UNF Matrix: " << unfMatrixName.str() << std::endl;
    }
  }


  for (int ibin=0; ibin<nr_of_bins; ibin++) {
    std::ostringstream covMatrixName;
    covMatrixName << "errorCov_" << ibin+1;
    TMatrixD * matrix = (TMatrixD*) fin->Get(covMatrixName.str().c_str());    
    if (matrix) {
      TMatrixD * mc = new TMatrixD(*matrix);
      //      mc->SetDirectory(0);
      binCovariance.insert(pair<int, TMatrixD*>(ibin+1,mc));
    } else {
      std::cout << "DID NOT FIND ERR. Matrix: " << covMatrixName.str() << std::endl;
    }
  }
    

  //  fin->Close();

  //  return 1;





  int verbosityLevel=0;

  // Read measurements from input

  
  std::map<int, TH1D* >::iterator ithist, ithist2;

  std::map<int, TH1D* >   differential_xsections;

  for (ithist=h_dsigma.begin(); ithist!=h_dsigma.end(); ithist++) {

    int channel = ithist->first;
    if (verbosityLevel>0) {
    std::cout << "CHANNEL = " << channel << std::endl
	      << "============== \n";
    }
    TH1D * h = ithist->second;
    int nbins = h->GetNbinsX();

    std::ostringstream diffXsName;
    //    diffXsName << "h_diffXS_" << channel;
    diffXsName << "h_crossSection_incl_diff_" << channel;
    TH1D * h_diffxs = (TH1D*) h->Clone(diffXsName.str().c_str());

    for (int ibin=0; ibin<=nbins+1; ibin++) {
      double val = h->GetBinContent(ibin);
      double wid = h->GetBinWidth(ibin);
      double err = h->GetBinError(ibin);

      h_diffxs->SetBinContent(ibin,val/wid);
      h_diffxs->SetBinError(ibin,err/wid);

      if (verbosityLevel>0) {
      std::cout << "Bin : " << ibin
		<< "  Val = " << val << "  Error = " << err 
		<< "  Val/dX = " << val/wid
		<< "  Rel. =  " << err/val ;
      std::cout << std::endl;
      }

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

	(*correlationMatrix)(irow,icol) = corr_xy;
	(*correlationMatrix)(icol,irow) = corr_xy;
	//  And now produce real unfolding covariance
	// Correctly normalized to cross section

	double unf_error_x = h_dsigma[channel]->GetBinError(irow+1);
	double unf_error_y = h_dsigma[channel]->GetBinError(icol+1);

	double cov_xy   = corr_xy*unf_error_x*unf_error_y;

	(*errorMatrix)(irow,icol) = cov_xy;
	(*errorMatrix)(icol,irow) = cov_xy;
      }
    }

    unfoldingCorrelation.insert(std::pair<int,TMatrixD*> ( it->first, correlationMatrix));
    unfoldingErrorMatrix.insert(std::pair<int,TMatrixD*> ( it->first, errorMatrix));

    if (verbosityLevel>1) {
      std::cout << "Correlation Matrix \n";
      correlationMatrix->Print();
      std::cout << "Error Matrix \n";
      errorMatrix->Print();
    }
  }

  //

  int dimension = nr_bins*nr_channels;

  //  vector<TMatrixD*>  unfoldingCovariance;   // [nr_channels]
  //  vector<TMatrixD*>  channelCovariance;     // [nr_bins]

  TMatrixD  fullCovariance(dimension,dimension);
  //  TMatrixD  lumiCovariance(dimension,dimension);
  TMatrixDSym  lumiCovariance(dimension);
  //  TMatrixDSym  miniLumiCov(4);

  //  ROOT::Math::SMatrix< double, 16, 16, ROOT::Math::MatRepStd<double,16,16 > > lumiCovarianceSM;
  //  AlgebraicSymMatrix16x16 lumiCovarianceSM;

  //  TMatrixD  lumiCovariance(dimension);
  TMatrixD  statCovariance(dimension,dimension);

  // Build full covariance matrix

  for (int ibin1=0; ibin1<nr_bins; ibin1++) {
    for (int ibin2=0; ibin2<nr_bins; ibin2++) {

      // SHIFT TO JUMP OVER UNDERFLOW BIN
      int binIndex = ibin1+1;

      for (int ichan1=0; ichan1<nr_channels; ichan1++) {
	for (int ichan2=0; ichan2<nr_channels; ichan2++) {


	  double sig1 = h_dsigma[ichan1]->GetBinContent(ibin1+1);
	  double sig2 = h_dsigma[ichan2]->GetBinContent(ibin2+1);
	  double cov(0);
	  double lumiCov(0);
	  double statCov(0);

	  // std::cout << "ibin1: " << ibin1
	  // 	    << "\t ibin2: " << ibin2
	  // 	    << "\t ichan1: " << ichan1
	  // 	    << "\t ichan2: " << ichan2
	  // 	    << "\t index: " << binIndex << std::endl;

	  if (ibin1 == ibin2 && ichan1 == ichan2) {
	    cov     = (*binCovariance[binIndex])(ichan1,ichan2);
	    lumiCov = 0.026*0.026*sig1*sig2;
	    statCov = (*unfoldingErrorMatrix[ichan1])(ibin1,ibin2);
	  }
	  else if (ibin1 == ibin2) {
	    cov = (*binCovariance[binIndex])(ichan1,ichan2);
	    lumiCov = 0.026*0.026*sig1*sig2;

	  } else if (ichan1 == ichan2) {
	    if (useFullCorrelations) {
	      cov     = (*unfoldingErrorMatrix[ichan1])(ibin1,ibin2);
	      statCov = cov;
	      lumiCov = 0.026*0.026*sig1*sig2;
	    }
	  } else {
	    cov = 0.;
	    //	    lumiCov = 0.026*0.026*sig1*sig2;
	  }
	  int index1 = index_from_bin_and_channel(ibin1,ichan1);
	  int index2 = index_from_bin_and_channel(ibin2,ichan2);
	  fullCovariance(index1,index2) = cov;
	  //if (index1<=index2) 
	  lumiCovariance(index1,index2) = lumiCov;
	  //	  lumiCovarianceSM(index1,index2) = lumiCov;
	  statCovariance(index1,index2) = statCov;

	  // 
	  // if (index1<4 && index2<4) {
	  //   miniLumiCov(index1,index2) = lumiCov;
	  // }

	}  
      }    
    }
  }

  // std::cout << "MINI LUMI : Det = " << miniLumiCov.Determinant() << std::endl;
  // double minidet;
  // miniLumiCov.Print();
  // miniLumiCov.Invert(&minidet);
  // miniLumiCov.Print();
  // //  CheckMatrixInvertible(miniLumiCov);
  // std::cout << "END OF MINI : Det = " << minidet << std::endl;


  // std::cout << "FULL BLOW COVARIANCE : \n";
  // std::cout << "====================== \n \n";
  // fullCovariance.Print();
  //  CheckMatrixInvertible(lumiCovariance);

  // std::cout << "FULL BLOWN LUMI COVARIANCE : Det -= " <<   lumiCovariance.Determinant() << std::endl;
  // std::cout << "========================= \n \n";
  // lumiCovariance.Print();
  // std::cout << "Full Covariance Determinant = " << fullCovariance.Determinant() << std::endl;

  //  std::cout << "S-matrix rep. : " << lumiCovarianceSM << std::endl;

  // int ifail=0;
  // AlgebraicSymMatrix16x16 lumiInverseSM = lumiCovarianceSM.Inverse(ifail);
  // AlgebraicMatrix16x16 product  = lumiCovarianceSM * lumiInverseSM;


  // std::cout << "Inverted S-matrix rep. : " << ifail << std::endl << lumiInverseSM 
  // 	    << "\n Product : \n: " << product << std::endl;
  // Build full measurements vector

  //  int nall = nr_bins*nr_channels;

  TVectorD  a(dimension);
  for (int ibin=0; ibin<nr_bins; ibin++) {
    for (int ichan=0; ichan<nr_channels; ichan++) {
      double val    = h_dsigma[ichan]->GetBinContent(ibin+1);
      double valdx  = differential_xsections[ichan]->GetBinContent(ibin+1);
      int index = index_from_bin_and_channel(ibin,ichan);
      a(index) = val;
    }
  }

  if (verbosityLevel>1) {
    std::cout << "MEASUREMENT VECTOR: \n";
    std::cout << "====================== \n \n";
    a.Print();
  }

  // Print out measurements with errors
  for (int ibin=0; ibin<nr_bins && verbosityLevel>0; ibin++) {
    double binWidth = h_dsigma[0]->GetBinWidth(ibin+1);
    for (int ichan=0; ichan<nr_channels; ichan++) {
      int index = index_from_bin_and_channel(ibin,ichan);
      double val = a(index);
      double err = fullCovariance(index,index);
      double stat_err = statCovariance(index,index);
      double lumi_err = lumiCovariance(index,index);
      std::cout << "Measurement in bin " << ibin << "\t channel: " << ichan
		<< "\t:\t" << val << " +/- " << sqrt(err)/binWidth
		<< "\t " << sqrt(stat_err)/binWidth << " (stat.)"
		<< "\t " << sqrt(lumi_err)/binWidth << " (lumi.)" 
		<< std::endl;
    }
    std::cout << "========================================================== \n";
  }

  // Invert full covariance matrix and auxiliary vector
  // 
  // M = Cov-1
  // 
  // f = M*a

  //  std::cout << "LUMI COVARIANCE (AGAIN) \n";
  //  lumiCovariance.Print();
  //  DumpMatrixToFile(lumiCovariance, "lumicov.txt");
  TMatrixDSym lumiCovCopy(lumiCovariance);
  TMatrixDSym  sima =  lumiCovariance.Invert();
  std::cout << "inverting lumi Covariance \n";
  TMatrixD * Mlumi = new TMatrixD(lumiCovariance);
  //  lumiCovariance.Print();
  //  Mlumi->Invert();
  //  Mlumi->Print();
  // Sanity check: multiply with the inverse: do we get unity???
  //  TMatrixD testUnity = lumiCovCopy*lumiCovariance; // lumiCovCopy;
  //  std::cout << "Should be unity: M*Minv = \n";
  //  testUnity.Print();

  std::cout << "inverting full Covariance \n";
  TMatrixD M = fullCovariance.Invert();

  //  lumiCovariance.Print();
  std::cout << "inverting stat Covariance \n";
  TMatrixD Mstat = statCovariance.Invert();
  TVectorD f = M*a;

  // std::cout << "LUMI COV INVERSE : \n";
  // Mlumi.Print();
  // 

  TMatrixD * finalCov     = GetFinalCovariance(nr_bins,nr_channels,M);
  TMatrixD * finalCovLumi = GetFinalCovariance(nr_bins,nr_channels,*Mlumi);
  TMatrixD * finalCovStat = GetFinalCovariance(nr_bins,nr_channels,Mstat);

  std::cout << "First inversion done \n";

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

  /*
  TMatrixD Dinv=D.Invert();
  */

  std::cout << "Second inversion done \n";
  
  //  result = Dinv * g;
  result = (*finalCov) * g;

  //  finalCovariance = Dinv;
  finalCovariance = *finalCov;

  //  Dinv.Print();

  if (verbosityLevel>0) {
    finalCov->Print();
    std::cout << " Lumi. Covariance: \n";
    finalCovLumi->Print();
    std::cout << " Stat. Covariance: \n";
    finalCovStat->Print();
  }
  // Get bin widths to correct for it...
  double * binWidths(0);

  //  std::map<int, TH1D* >::iterator ithist;
  ithist = h_dsigma.find(0);
  if (ithist != h_dsigma.end() ) {
    int nbins = ithist->second->GetNbinsX();
    binWidths = new double[nbins];
    for (int ibin=0; ibin<nbins; ibin++) {
      binWidths[ibin] = ithist->second->GetBinWidth(ibin+1);
    }
  } else {
    std::cout << "DID NOT FIND HISTOGRAM TO GET WIDTHS... \n";
  }


  TH1D* h_combined_ds    = (TH1D*) h_dsigma[0]->Clone("h_xs_comb");
  TH1D* h_combined_dsdx  = (TH1D*) h_dsigma[0]->Clone("h_xs_comb_diff");
  h_combined_ds->Reset();
  h_combined_dsdx->Reset();
  h_combined_dsdx->Print();
  //  h_combined_dsdx->SetNameTitle("hComb_diff","Combined diff. xs");
  //  h_combined_dsdx->SetName("hComb_diff");

  ofstream outCombined;
  std::ostringstream outCombName;
  outCombName << "outCombination-" << variable << ".txt";
  outCombined.open(outCombName.str().c_str());
  std::cout << "FINAL RESULTS: \n \n";
  for (int ibin=0; ibin<nr_bins ;ibin++) {
    double val      = result(ibin);
    double err      = sqrt(finalCovariance(ibin,ibin));

    h_combined_ds->SetBinContent(ibin+1,val);
    h_combined_ds->SetBinError(ibin+1,err);


    double err_stat = sqrt( (*finalCovStat)(ibin,ibin));
    //    double err_lumi = sqrt( (*finalCovLumi)(ibin,ibin));
    double err_lumi = 0.026*val;
    double err_syst  = sqrt(err*err - err_stat*err_stat - err_lumi*err_lumi);
    if (binWidths) {
      val      *= (1./binWidths[ibin]);
      err      *= (1./binWidths[ibin]);
      err_stat *= (1./binWidths[ibin]);
      err_lumi *= (1./binWidths[ibin]);
      err_syst *= (1./binWidths[ibin]);
    }
    std::cout << "Bin " << ibin 
	      << "\t Value = " << val
      //	      << " +/- " << sqrt(finalCovariance(ibin,ibin)) << std::endl;
	      << " +/- " << err
	      << "\t: " << err_stat << "(stat)"
	      << "\t" << err_syst << "(syst)"
	      << "\t" << err_lumi << "(lumi)"
	      << std::endl;
    outCombined << ibin+1 << "\t" << val
		<< "\t" << err_stat 
		<< "\t" << err_syst 
		<< "\t" << err_lumi 
		<< std::endl;

    h_combined_dsdx->SetBinContent(ibin+1,val);
    h_combined_dsdx->SetBinError(ibin+1,err);

  }
  outCombined.close();

  h_combined_dsdx->Print();
  fout->cd();
  h_combined_ds->Write();
  h_combined_dsdx->Write();

  for (ithist=h_dsigma.begin(); ithist!=h_dsigma.end(); ithist++) {
    ithist->second->Write();
    //    h_dsigmadx[ithist->first]->Write();
    differential_xsections[ithist->first]->Write();
  }
  fout->Close();

  return 1;


  TApplication *theApp = new TApplication("app", &argc, argv); //<======

  // TStyle * m_gStyle;
  // m_gStyle = new TStyle();
  // m_gStyle->SetOptStat("0000");
  TCanvas * can = new TCanvas();
  can->Divide(2,2);
  for (int ichan=0; ichan<4; ichan++) {
    can->cd(ichan+1);
    unfoldingCorrelation[ichan]->Draw("colztext");
  }

  std::ostringstream correlationNamePdf, correlationNameRoot;
  correlationNamePdf  << "unfoldingCorrelation-" << variable << ".pdf";
  correlationNameRoot << "unfoldingCorrelation-" << variable << ".root";

  can->SaveAs(correlationNamePdf.str().c_str());
  can->SaveAs(correlationNameRoot.str().c_str());
  //  can->SaveAs("unfoldingCorrelation.root");

  if (showPlots)   theApp->Run();

 
}
