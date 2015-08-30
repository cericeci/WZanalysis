#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMatrixD.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <set>

#include "UnfoldingHistogramFactory.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////
//    BLUE Combination between the 4 WZ channels bin by bin
//       for unfolded distributions
//
//     For each channel independently:
//
//      Build Covariance Matrix ErrMat = E
//
//    xs(comb) = sum_i alpha_i xs(ch. i)
//
//   Var(xs(comb) ) = alpha_t * E * alpha
//
//         alpha = 1/N * ( E-1 * U)
// where:
//      N = U_t * E-1 * U
//      U = (1,1,1,1)
//
// Constructing E:
//
//   diagonal elements = sum_ syst_i^2 + stat^2 + lumi^2
//
//   off-diagonal elemnents (k.l): sum_(syst. sources i) sigma(
//   ( this is assuming full correlation: cov(x,y) = sigma_x*sigma_y )
//

int main(int argc, char **argv)
{
  //read histograms form .root files
  //definirati file-ove iz kojih citam i u koje pisem
  bool gotHistoBinning(false);
  bool latexOutput(true);  
  char * binningFileName(0);
  char * variableName(0);
  bool gotVarName  = false;
  char c;
  ofstream outMain, outError1, outError2, outError3;
  outMain.open("outMain.txt");
  outError1.open("outError1.txt");
  outError2.open("outError2.txt");
  outError3.open("outError3.txt");

  while ((c = getopt (argc, argv, "v:H:")) != -1)
    switch (c)
      {
      case 'v':
	gotVarName = true;
	variableName = new char[strlen(optarg)+1];
	strcpy(variableName,optarg);
	break;
      case 'H':
	gotHistoBinning = true;
	binningFileName = new char[strlen(optarg)+1];
	strcpy(binningFileName,optarg);
	break;
      default:
	std::cout << "usage: -r responseFile [-d <dataFile>]   \n";
	abort ();
      }
  
  string variable;
  
  if (gotVarName) {
    variable = variableName;
  } else {
    variable = "LeadJetPt";
  }

  if (variable != "Njets" && variable != "LeadingJetPt" && variable !="Zpt"){
    std::cout<<"UNKNOWN VARIABLE!!!"<<std::endl;
    return -1;
  }

  bool printBLUEmatrix(true);



  std::ostringstream outfilename, infilename;
  outfilename<<"unfoldingFinalResults/combination_"<<variable<<".root";
  infilename<<"sysResults/systematics_"<<variable<<".root";
  TFile * fout= new TFile(outfilename.str().c_str(), "RECREATE");  
  TFile * finput= TFile::Open(infilename.str().c_str());
  if (gotHistoBinning) {

    UnfoldingHistogramFactory * histoFac = UnfoldingHistogramFactory::GetInstance();
    histoFac->SetBinning(binningFileName);

  }

  std::vector<string> sysSources;
  std::vector<string> allHistos;

  sysSources.push_back("qcdScale");
  sysSources.push_back("PDFsys");
  //  sysSources.push_back("leptTrgEff_el");
  //  sysSources.push_back("leptTrgEff_mu");
  sysSources.push_back("ele_SF");
  sysSources.push_back("mu_SF");
  sysSources.push_back("Etsys");
  sysSources.push_back("muMomScale");
  sysSources.push_back("elEnScale");
  sysSources.push_back("pileupSys");
  sysSources.push_back("ZZxs");
  sysSources.push_back("Zgammaxs");
  //  sysSources.push_back("dataDrivensys");
  sysSources.push_back("dataDrivenElsys");
  sysSources.push_back("dataDrivenMusys");
  sysSources.push_back("bckgSys");
  sysSources.push_back("unfSyst");


  if (variable == "LeadingJetPt" || variable =="Njets") {
    sysSources.push_back("JESsys");
    sysSources.push_back("JERsys");
  }

  // copy all of it in the 2nd vector
  allHistos.insert(allHistos.begin(), sysSources.begin(), sysSources.end());

  //  for (int isyst=0; isyst<sysSources.size(); isyst++) {
  //    allHistos.push_back(sysSources[isyst]);
  //  }

  allHistos.push_back("lumi");
  allHistos.push_back("crossSection");
  allHistos.push_back("crossSection_inclusive");
  allHistos.push_back("crossSection_incl_diff");
  allHistos.push_back("totalSyst");
  allHistos.push_back("totalSyst_diff");
  allHistos.push_back("totalStat");
  allHistos.push_back("totalStat_diff");


  for (int ivar=0; ivar<allHistos.size(); ivar++) {
    std:: cout << "Variable : " << allHistos[ivar] << std::endl;
  }
  
  const int nChannels(4);

  std::map<string, std::vector<TH1D *> > inputHistos;

  TH1D * h_crossSection_combination;
  TH1D * h_crossSection_comb_diff;
  TH1D * h_combStat;
  TH1D * h_combSyst;
  TH1D * h_combLumi;
  TH1D * h_combStat_diff;
  TH1D * h_combSyst_diff;

  TH1D * h_totalSyst[nChannels];
  TH1D * h_totalSyst_diff[nChannels];
  TH1D * h_totalStat[nChannels];
  TH1D * h_totalStat_diff[nChannels];



  double crossSection[4]={0,0,0,0};
  
  //ovo treba deklarirat!!!!!!!!!
  for (int fill=0; fill<nChannels; fill++) {
    std::ostringstream totalSysName, totalSysDiffName, totalStatName, totalStatDiffName;
    totalSysName<<"h_totalSyst_"<<fill;
    totalSysDiffName<<"h_totalSyst_diff_"<<fill;
    totalStatName<<"h_totalStat_"<<fill;
    totalStatDiffName<<"h_totalStat_diff_"<<fill;

    h_totalSyst[fill]=UnfoldingHistogramFactory::createHistogramForVar(variable, 
								       totalSysName.str().c_str(), 
								       totalSysName.str().c_str());
    h_totalSyst_diff[fill]=UnfoldingHistogramFactory::createHistogramForVar(variable,
									    totalSysDiffName.str().c_str(),
									    totalSysDiffName.str().c_str());
    h_totalStat[fill]=UnfoldingHistogramFactory::createHistogramForVar(variable,
								       totalStatName.str().c_str(), 
								       totalStatName.str().c_str());
    h_totalStat_diff[fill]=UnfoldingHistogramFactory::createHistogramForVar(variable,
									    totalStatDiffName.str().c_str(),
									    totalStatDiffName.str().c_str());

  }


  h_crossSection_combination= UnfoldingHistogramFactory::createHistogramForVar(variable,
									       "h_xs_comb", "h_xs_comb");
  h_crossSection_comb_diff= UnfoldingHistogramFactory::createHistogramForVar(variable,
									     "h_xs_comb_diff", "h_xs_comb_diff");
  h_combStat = UnfoldingHistogramFactory::createHistogramForVar(variable,
								"h_combStat", "h_combStat"); 
  h_combSyst = UnfoldingHistogramFactory::createHistogramForVar(variable,
								"h_combSyst", "h_combSyst");
  h_combLumi = UnfoldingHistogramFactory::createHistogramForVar(variable,
								"h_combLumi", "h_combLumi");
  h_combStat_diff = UnfoldingHistogramFactory::createHistogramForVar(variable,
								     "h_combStat_diff", "h_combStat_diff"); 
  h_combSyst_diff = UnfoldingHistogramFactory::createHistogramForVar(variable,
								     "h_combSyst_diff", "h_combSyst_diff");


  for (int hist=0; hist<nChannels; hist++){
    for (int ivar=0; ivar<allHistos.size(); ivar++) {

      std::ostringstream histoKey;
      std::string        inputKey;
      if (allHistos[ivar] != "ele_SF" && allHistos[ivar] != "mu_SF" ) {
	histoKey << "h_";
      }
      histoKey << allHistos[ivar];
      if (allHistos[ivar] != "crossSection_inclusive") {
	histoKey << "_";
      }
      histoKey << hist;
      if (allHistos[ivar] != "crossSection_incl_diff") {
	inputKey = histoKey.str();
      } else {
	std::ostringstream histoKey_input;
	histoKey_input << "h_crossSection_inclusive" << hist;
	inputKey = histoKey_input.str();
      }

      TH1D * h = (TH1D*) finput ->Get(inputKey.c_str());
      if ( h!= 0) {
	//	  inputHistos[allHistos][hist] = (TH1D*) h->Clone(histoKey.str().c_str());
	inputHistos[allHistos[ivar]].push_back( (TH1D*) h->Clone(histoKey.str().c_str()));
      } else {
	std::cout << "NONEXISTENT HISTOGRAM: " << inputKey << "\t var = " << allHistos[ivar]
		  << "\t exiting... "<< std::endl;
	exit;
      } 
    }
  }

  //loop over each bin 

  int nBins = (inputHistos[allHistos[0]][0])->GetNbinsX();
  std::map<std::string,TMatrixD *> covarianceMatrices;


  for (int bin=1; bin< nBins +1; bin++) {
      double elements[16];
      for (int el=0; el<16; el++) elements[el]=0;
      double statisticError[nChannels];
      double systematicError2[nChannels];
      double systematicError[nChannels];
      double lumiError[nChannels];
      //statistic and systematic errors:

      for (int nCh=0; nCh<nChannels; nCh++) {

	//// DOUBLE CHECK THE NEXT COUPLE OF LINES

	//      statisticError[nCh]=h_crossSection[nCh]->GetBinError(bin);
	statisticError[nCh]=   inputHistos["crossSection_inclusive"][nCh]->GetBinError(bin);
	//h_crossSection_final[nCh]->GetBinError(bin);
	double sys_err=0;
	for (int isys=0; isys<sysSources.size(); isys++ ) {

	  sys_err += pow( (inputHistos[sysSources[isys]][nCh])->GetBinContent(bin),2);

	}
	std::cout << "FIRST PRINTOUT \n";
	systematicError2[nCh]=sys_err;

	h_totalStat[nCh]->SetBinContent(bin, statisticError[nCh]);
	systematicError[nCh]     = (sqrt(systematicError2[nCh]))
	  *((inputHistos["crossSection_inclusive"][nCh])->GetBinContent(bin));
	lumiError[nCh]           = (inputHistos["lumi"][nCh])->GetBinContent(bin)
	  *( (inputHistos["crossSection_inclusive"][nCh])->GetBinContent(bin));

	h_totalSyst[nCh]->SetBinContent(bin, (systematicError[nCh]));
	std::cout<<bin<<" , "<<systematicError[nCh]<<" , "
		 << (inputHistos["crossSection_inclusive"][nCh])->GetBinContent(bin)<<std::endl;
      }
      //common elements
      double commonSys[4][4];
      double lumiSys[4][4];

      std::vector<string> commonSources;
      commonSources.push_back("qcdScale");
      commonSources.push_back("PDFsys");
      commonSources.push_back("Etsys");
      commonSources.push_back("pileupSys");
      commonSources.push_back("ZZxs");
      commonSources.push_back("Zgammaxs");
      commonSources.push_back("bckgSys");
      commonSources.push_back("unfSyst");
      commonSources.push_back("lumi");
      if (variable!="Zpt") {
	commonSources.push_back("JESsys");
      }


      for (int cha=0; cha<4; cha++){
	for (int chb=0; chb<4; chb++){
	  double common_sys=0;
	  for (int isys=0; isys<commonSources.size(); isys++) {
	    common_sys += (inputHistos[commonSources[isys]][cha])->GetBinContent(bin)
	      * (inputHistos[commonSources[isys]][chb])->GetBinContent(bin);
	  }
	  commonSys[cha][chb] = common_sys;

	  lumiSys[cha][chb]= 
	    (inputHistos["lumi"][cha])->GetBinContent(bin) *
	    (inputHistos["lumi"][chb])->GetBinContent(bin) *
	    (inputHistos["crossSection_inclusive"][cha])->GetBinContent(bin) *
	    (inputHistos["crossSection_inclusive"][chb])->GetBinContent(bin);

	  std::cout << "common sys [" << cha << " , " << chb << " ] = " << 
	    commonSys[cha][chb] << std::endl;


	}
      }
    
      //diagonal elements
      double wo_lumi_1=pow(systematicError[0],2) + pow(statisticError[0],2);
      double wo_lumi_2=pow(systematicError[1],2) + pow(statisticError[1],2);
      double wo_lumi_3=pow(systematicError[2],2) + pow(statisticError[2],2);
      double wo_lumi_4=pow(systematicError[3],2) + pow(statisticError[3],2);

      elements[0]=pow(systematicError[0],2) + pow(statisticError[0],2) +pow(lumiError[0],2);
      elements[5]=pow(systematicError[1],2) + pow(statisticError[1],2) +pow(lumiError[1],2);
      elements[10]=pow(systematicError[2],2) + pow(statisticError[2],2) +pow(lumiError[2],2);
      elements[15]=pow(systematicError[3],2) + pow(statisticError[3],2) +pow(lumiError[3],2);

      inputHistos["crossSection_inclusive"][0]->SetBinError(bin, sqrt(wo_lumi_1));
      inputHistos["crossSection_inclusive"][1]->SetBinError(bin, sqrt(wo_lumi_2));
      inputHistos["crossSection_inclusive"][2]->SetBinError(bin, sqrt(wo_lumi_3));
      inputHistos["crossSection_inclusive"][3]->SetBinError(bin, sqrt(wo_lumi_4));

      double corr_matrix_dd[4][4]={
	{1, 0, 1, 0},
	{0, 1, 0, 1}, 
	{1, 0, 1, 0},
	{0, 1, 0, 1}
      };
      //matrix is symetric
      //channels 0 and 1

      vector<string> uncorrelatedSources;

      uncorrelatedSources.push_back("elEnScale");
      uncorrelatedSources.push_back("muMomScale");
      uncorrelatedSources.push_back("ele_SF");
      uncorrelatedSources.push_back("mu_SF");
      uncorrelatedSources.push_back("dataDrivenElsys");
      uncorrelatedSources.push_back("dataDrivenMusys");

      TMatrixD  auxiliaryMatrix(4,4);

      for (int irow=0; irow<4; irow++) {


	for (int icol=irow; icol<4; icol++) {
	  if (irow==icol) continue;  // Diagonal elements treated specially earlier

	  double x = commonSys[irow][icol];

	  for (int isys=0; isys<uncorrelatedSources.size(); isys++) {
	    std::string source=uncorrelatedSources[isys];
	    x += (inputHistos[source][irow]->GetBinContent(bin))
	      *(inputHistos[source][icol]->GetBinContent(bin));
	    if (irow==0 && icol==2 && false) {
	      std::cout << source << " = " << (inputHistos[source][irow]->GetBinContent(bin))
			<< " " << (inputHistos[source][icol]->GetBinContent(bin)) << std::endl;
	    }
	  }
	  x *= (inputHistos["crossSection_inclusive"][irow]->GetBinContent(bin))
	    * (inputHistos["crossSection_inclusive"][icol]->GetBinContent(bin));

	  // Matrix is symmetric
	  auxiliaryMatrix(irow,icol) = x;
	  auxiliaryMatrix(icol,irow) = x;	
	}
      }
      //  Matrix ordering for 1D indices
      //     0   4   8  12    
      //     1   5   9  13
      //     2   6  10  14
      //     3   7  11  15
      elements[4]  = elements[1] = auxiliaryMatrix[0][1];
      elements[8]  = elements[2] = auxiliaryMatrix[0][2];
      elements[12] = elements[3] = auxiliaryMatrix[0][3];
      elements[9]  = elements[6] = auxiliaryMatrix[1][2];
      elements[13] = elements[7] = auxiliaryMatrix[1][3];
      elements[11] = elements[14]= auxiliaryMatrix[2][3];
    
      if (printBLUEmatrix) {
	std::cout<<"Error matrix finished!!"<<std::endl;
	for (int emi=0;emi<4;emi++) {
	  for (int emj=0;emj<4;emj++) {
	    std:: cout << "       " << elements[emj*4+emi] << "       ";
	  }
	  std::cout << endl;
	}
	std::cout << endl; 
      }
      TMatrixD errMat(4,4,elements);
      TMatrixD errMatInv(4,4);
      TMatrixD errMatCopy(errMat);

      std::ostringstream matrixName;
      matrixName << "covMatrix_" << variable << "_" << bin;

      covarianceMatrices.insert(std::pair<std::string,TMatrixD*>
				(matrixName.str(), new TMatrixD(errMat)));


      //invert !
      // NOTE: after a=b.Invert()
      // Inverted matrix is contained in both a and b !!!
      // The original matrix is lost (but was copied in errMatCopy above)
      errMatInv = errMat.Invert();
      Double_t *mRefTest = errMat.GetMatrixArray();
      if (printBLUEmatrix) {
	std::cout << endl << " Matrix Inverse:" << endl;
	for (int emi=0;emi<4;emi++) {
	  for (int emj=0;emj<4;emj++) {
	    std::cout << "       " << mRefTest[emj*4+emi] << "       ";
	  }
	  std::cout << endl;
	}
	std::cout << endl;
      }
      //get norm, and alpha factors for each channel
      Double_t *mRef= errMat.GetMatrixArray();
      Double_t norm=0.;
      Double_t alphaCH[4]={0.,0.,0.,0.};
    
      for (int imatrix=0;imatrix<16;imatrix++) norm+=mRef[imatrix];
    
      for (size_t im=0;im<4;im++) {
	for (size_t jm=0;jm<4;jm++) {
	  alphaCH[im]+=mRef[im*4+jm];
	}
	alphaCH[im]/=norm;
      }
      if (printBLUEmatrix){
	std::cout << "al0 " << alphaCH[0] << " al1 " << alphaCH[1] << " al2 " << alphaCH[2] << " al3 " << alphaCH[3] << endl;
	std::cout << "consistency check:" << alphaCH[0]+alphaCH[1]+alphaCH[2]+alphaCH[3] <<endl;
	std::cout << endl;
      }
      double final_Xsec = 0;
      for (int ich=0; ich<4; ich++) {
	final_Xsec += alphaCH[ich]*(inputHistos["crossSection_inclusive"][ich]->GetBinContent(bin));
      }
      // double final_Xsec = alphaCH[0]*(h_crossSection_final[0]->GetBinContent(bin)) + alphaCH[1]*(h_crossSection_final[1]->GetBinContent(bin)) 
      //   + alphaCH[2]*(h_crossSection_final[2]->GetBinContent(bin)) + alphaCH[3]*(h_crossSection_final[3]->GetBinContent(bin));
    
      Double_t combined_error=0;
      Double_t *copyRef = errMatCopy.GetMatrixArray();
      for (int ier=0;ier<4;ier++)
	for (int jer=0;jer<4;jer++) {
	  combined_error+=alphaCH[ier]*alphaCH[jer]*copyRef[ier*4+jer];
	}
      //lumi error
      double lumi_error2(0);
      for (int ier=0; ier<4; ier++){
	for (int jer=0; jer<4; jer++){
	  lumi_error2+=alphaCH[ier]*alphaCH[jer]*lumiSys[ier][jer];
	}
      }

      Double_t stat_err_tot2 =  pow(alphaCH[0]*statisticError[0],2) 
	+ pow(alphaCH[1]*statisticError[1],2)
	+ pow(alphaCH[2]*statisticError[2],2) 
	+ pow(alphaCH[3]*statisticError[3],2);
    
      Double_t stat_err_tot = sqrt(stat_err_tot2);
      Double_t syst_err_tot = sqrt (combined_error - stat_err_tot2- lumi_error2);
    
      combined_error=sqrt(combined_error);
    
      if(printBLUEmatrix){
	std::cout << "combined sigma(WZ) = " << final_Xsec << " +- " << stat_err_tot << "(stat.) +- " << syst_err_tot << "(syst) +- "
		  << sqrt(lumi_error2) << " (lumi) pb " << endl;
	std::cout << endl;
      }
      //double total_error=sqrt(stat_err_tot*stat_err_tot + syst_err_tot*syst_err_tot + sqrt(lumi_error2));
      double lumi_error=sqrt(lumi_error2);
      /*
	std::cout<<"COMBINED ERROR:"<<combined_error<<std::endl;
	std::cout<<"STAT ERROR: "<<stat_err_tot<<std::endl;
	std::cout<<"SYST ERROR: "<<syst_err_tot<<std::endl;
	std::cout<<"TOTAL ERROR:"<<total_error<<std::endl;
	std::cout<<"BIN:   "<<bin<<std::endl;
      */
      h_crossSection_combination->SetBinContent(bin, final_Xsec);
      //    h_crossSection_combination->SetBinError(bin, total_error);
      h_crossSection_combination->SetBinError(bin, combined_error);
      h_combStat->SetBinContent(bin, stat_err_tot);
      h_combSyst->SetBinContent(bin, syst_err_tot);
      h_combLumi->SetBinContent(bin, lumi_error);
    }   // End of loop over bins


	 for (int i=1; i<=h_crossSection_combination->GetNbinsX(); i++) {
	   double value = h_crossSection_combination->GetBinContent(i);
	   double error = h_crossSection_combination->GetBinError(i);
	   double errorStat = h_combStat->GetBinContent(i);
	   double errorSyst = h_combSyst->GetBinContent(i);
	   double errorLumi = h_combLumi->GetBinContent(i);
	   double width = h_crossSection_combination->GetBinWidth(i);
	   double dsdx = value/width;
	   double dsdx_err = dsdx*error/value;
	   double dsdx_errStat = dsdx*errorStat/value;
	   double dsdx_errSyst = dsdx*errorSyst/value;
	   double dsdx_errLumi = dsdx*errorLumi/value;  
	   h_crossSection_comb_diff->SetBinContent(i,dsdx);
	   h_crossSection_comb_diff->SetBinError(i,dsdx_err);
	   h_combStat->SetBinContent(i,dsdx_errStat);
	   h_combSyst->SetBinContent(i,dsdx_errSyst);
	   h_combLumi->SetBinContent(i,dsdx_errLumi);
	 }

	 for (int channels=0; channels<4; channels++){
	   for (int i=1; i<=(inputHistos["crossSection_inclusive"])[channels]->GetNbinsX(); i++) {
	     double value2 =(inputHistos["crossSection_inclusive"][channels])->GetBinContent(i);
	     double error2 = (inputHistos["crossSection_inclusive"][channels])->GetBinError(i);
	     double errorStat2 = h_totalStat[channels]->GetBinContent(i);
	     double errorSyst2 = h_totalSyst[channels]->GetBinContent(i);
	     double width2 = (inputHistos["crossSection_inclusive"][channels])->GetBinWidth(i);
	     double dsdx2 = value2/width2;
	     double dsdx_err2 = dsdx2*error2/value2;
	     double dsdx_errorStat= dsdx2*errorStat2/value2;
	     double dsdx_errorSyst= dsdx2*errorSyst2/value2;
	     (inputHistos["crossSection_incl_diff"][channels])->SetBinContent(i,dsdx2);
	     (inputHistos["crossSection_incl_diff"][channels])->SetBinError(i,dsdx_err2);
	     h_totalSyst_diff[channels]->SetBinContent(i, dsdx_errorSyst);
	     h_totalStat_diff[channels]->SetBinContent(i, dsdx_errorStat);
	   }
	 }
	 

	 //LATEX OUTPUT
	 if (latexOutput){
	   // TString ranges[9]={"","","","","","","","","",""};
	   // TString rangesZpt[9]={"0-20 GeV", "20-40 GeV", "40-60 GeV", "60-80 GeV", "80-100 GeV", "100-120 GeV", "120-140 GeV", "140-200 GeV", "200-300 GeV"};
	   // TString rangesLeadingJetPt[9]={"30-60 GeV", "60-100 GeV", "100-150 GeV", "150-250 GeV"};
	   // TString rangesNjets[9]={"0 jets", "1 jet", "2 jets", "3 jets", "4 jets"};

	   std::vector<std::string>  variableRanges;


	   if (variable == "Njets") {
	     variableRanges.push_back("0 jets");
	     variableRanges.push_back("1 jet");
	     variableRanges.push_back("2 jets");
	     variableRanges.push_back("3 jets");
	     variableRanges.push_back("4 jets");
	   } else if (variable == "LeadingJetPt") {
	     variableRanges.push_back("30-60 GeV");
	     variableRanges.push_back("60-100 GeV");
	     variableRanges.push_back("100-150 GeV");
	     variableRanges.push_back("150-250 GeV");
	   } else if ( variable =="Zpt"){
	     variableRanges.push_back("0-20 GeV");
	     variableRanges.push_back("20-40 GeV");
	     variableRanges.push_back(	"40-60 GeV");
	     variableRanges.push_back(	"60-80 GeV");
	     variableRanges.push_back( "80-100 GeV");
	     variableRanges.push_back( "100-120 GeV");
	     variableRanges.push_back( "120-140 GeV");
	     variableRanges.push_back( "140-200 GeV");
	     variableRanges.push_back( "200-300 GeV");
	   }

	   map<std::string,int> variablePrecision;
	   variablePrecision["Zpt"]          = 4;
	   variablePrecision["LeadingJetPt"] = 4;
	   variablePrecision["Njets"]        = 3;

	   std::cout<<"--------------------------------------------------------"<<std::endl;
	   std::cout<<"Latex output: "<<std::endl;
	   std::cout<<"bin & 3e & 2e1mu & 1e2mu & 3mu & combination \\\\"<<std::endl;
	   std::cout<<"\\hline"<<std::endl;

	   for (int output=1; output<=h_crossSection_combination->GetNbinsX(); output++) {

	     outMain<<output<<" "
		    << (inputHistos["crossSection_incl_diff"][0])->GetBinContent(output) << " "
		    << (inputHistos["crossSection_incl_diff"][1])->GetBinContent(output) << " "
		    << (inputHistos["crossSection_incl_diff"][2])->GetBinContent(output) << " "
		    << (inputHistos["crossSection_incl_diff"][3])->GetBinContent(output) << " "
		    <<h_crossSection_comb_diff->GetBinContent(output)<<std::endl; 

	     outError1 << output<<" "
		       << h_totalStat_diff[0]->GetBinContent(output)<<" "
		       << h_totalStat_diff[1]->GetBinContent(output)<<" "
		       << h_totalStat_diff[2]->GetBinContent(output)<<" "
		       << h_totalStat_diff[3]->GetBinContent(output)<<" "
		       << h_combStat->GetBinContent(output)<<std::endl;

	     outError2 << output<<" "
		       << h_totalSyst_diff[0]->GetBinContent(output)<<" "
		       << h_totalSyst_diff[1]->GetBinContent(output)<<" "
		       << h_totalSyst_diff[2]->GetBinContent(output)<<" "
		       << h_totalSyst_diff[3]->GetBinContent(output)<<" "
		       << h_combSyst->GetBinContent(output)<<std::endl; 

	     outError3 << output<<" "
		       << 0.026*(inputHistos["crossSection_incl_diff"][0])->GetBinContent(output) << " "
		       << 0.026*(inputHistos["crossSection_incl_diff"][1])->GetBinContent(output) << " "
		       << 0.026*(inputHistos["crossSection_incl_diff"][2])->GetBinContent(output) << " "
		       << 0.026*(inputHistos["crossSection_incl_diff"][3])->GetBinContent(output) << " "
		       << h_combLumi->GetBinContent(output)<<std::endl; 

	     if (variable == "Zpt") {
	       std::cout<<scientific ; 
	     } else {
	       std::cout<<fixed;
	     }
	     std::cout << setprecision(variablePrecision[variable])
		       << variableRanges[output-1] <<" & ";
	     for (unsigned int ich=0; ich<4; ich++) {

	       std::cout << (inputHistos["crossSection_incl_diff"][ich])->GetBinContent(output) <<" $\\pm$ "
			 << h_totalStat_diff[ich]->GetBinContent(output)<<" $\\pm$ "
			 << h_totalSyst_diff[ich]->GetBinContent(output)<<" $\\pm$ "
			 << 0.026*(inputHistos["crossSection_incl_diff"][ich])->GetBinContent(output) <<" & ";

	     }
	     std::cout << h_crossSection_comb_diff->GetBinContent(output)<<" $\\pm$ "
		       << h_combStat->GetBinContent(output)<<" $\\pm$ "
		       << h_combSyst->GetBinContent(output)<<"$\\pm$ "
		       << h_combLumi->GetBinContent(output)<<"\\\\"<<std::endl;

	   }
  
	   outMain.close();

	 }
	 //END OF LATEX OUTPUT
	 fout->cd();
       for (unsigned int ich=0; ich<4; ich++) {
	 (inputHistos["crossSection_inclusive"][ich])->Write();
	 (inputHistos["crossSection_incl_diff"][ich])->Write();
	 
       }
	 
	 
	 h_crossSection_combination->Write();
       h_crossSection_comb_diff->Write();
       h_combStat->Write();
       h_combSyst->Write();
       

       fout->Close();

       TFile * fmat = new TFile("matrices.root","recreate");
       fmat->cd();
       for (std::map<std::string,TMatrixD*>::iterator it = covarianceMatrices.begin();
	    it != covarianceMatrices.end(); it++ ) {
	 it->second->Write(it->first.c_str()); //
       }
    fmat->Close();
       
       }
