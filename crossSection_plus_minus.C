/*
#include "forXS/numbers.h"
#include "forXS/numMC.h"
#include "forXS/numMM_met.h"
#include "forXS/numGEN_met.h"
#include "forXS/numData_met.h"
*/

#include "forXS/numbers.h"
#include "forXS/numMC.h"
#include "forXS/numMM_met.h"
#include "forXS/numGEN_met.h"
#include "forXS/numData_met.h"
#include "forXS/syst.h"

#include "forXS/syst.h"
#include "VVV.h"    //for systematics
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

void crossSection_plus_minus()
{
  //  bool latexOutput(false);
  bool latexOutput(true);
    //0: eee, 1:eem, 2:emm, 3:mmm
  ofstream fileNum;

  double AxeZ_minus[4]={dAxe3eZ_2_minus, dAxe2e1muZ_2_minus, dAxe1e2muZ_2_minus, dAxe3muZ_2_minus};
  double AxeZ_plus[4]={dAxe3eZ_2_plus, dAxe2e1muZ_2_plus, dAxe1e2muZ_2_plus, dAxe3muZ_2_plus};
  double sAxeZ_minus[4]={dsAxe3eZ_minus, dsAxe2e1muZ_minus, dsAxe1e2muZ_minus, dsAxe3muZ_minus};
  double sAxeZ_plus[4]={dsAxe3eZ_plus, dsAxe2e1muZ_plus, dsAxe1e2muZ_plus, dsAxe3muZ_plus};
  double N_good_minus[4]={dN_good3e_minus, dN_good2e1mu_minus, dN_good1e2mu_minus, dN_good3mu_minus};
  double sN_good_minus[4]={dsN_good3e_minus, dsN_good2e1mu_minus, dsN_good1e2mu_minus, dsN_good3mu_minus};
  double N_good_plus[4]={dN_good3e_plus, dN_good2e1mu_plus, dN_good1e2mu_plus, dN_good3mu_plus};
  double sN_good_plus[4]={dsN_good3e_plus, dsN_good2e1mu_plus, dsN_good1e2mu_plus, dsN_good3mu_plus};
  //  double N_fake[4]={dN_fake3e, dN_fake2e1mu, dN_fake1e2mu, dN_fake3mu};
  //double sN_fake[4]={dsN_fake3e, dsN_fake2e1mu, dsN_fake1e2mu, dsN_fake3mu};
  //  double N_data[4]={dN_data3e, dN_data2e1mu, dN_data1e2mu, dN_data3mu};
  double NZZ_minus[4]={dNZZ_3e_minus, dNZZ_2e1mu_minus, dNZZ_1e2mu_minus, dNZZ_3mu_minus};
  double NZZ_plus[4]={dNZZ_3e_plus, dNZZ_2e1mu_plus, dNZZ_1e2mu_plus, dNZZ_3mu_plus};
  double sNZZ_minus[4]={dsNZZ_3e_minus, dsNZZ_2e1mu_minus, dsNZZ_1e2mu_minus, dsNZZ_3mu_minus};
  double sNZZ_plus[4]={dsNZZ_3e_plus, dsNZZ_2e1mu_plus, dsNZZ_1e2mu_plus, dsNZZ_3mu_plus};
  double NZgamma_minus[4]={dNZgamma_3e_minus, dNZgamma_2e1mu_minus, dNZgamma_1e2mu_minus, dNZgamma_3mu_minus};
  double NZgamma_plus[4]={dNZgamma_3e_plus, dNZgamma_2e1mu_plus, dNZgamma_1e2mu_plus, dNZgamma_3mu_plus};
  double sNZgamma_minus[4]={dsNZgamma_3e_minus, dsNZgamma_2e1mu_minus, dsNZgamma_1e2mu_minus, dsNZgamma_3mu_minus};
  double sNZgamma_plus[4]={dsNZgamma_3e_plus, dsNZgamma_2e1mu_plus, dsNZgamma_1e2mu_plus, dsNZgamma_3mu_plus};
  double NWV_minus[4]={dNWV_3e_minus, dNWV_2e1mu_minus, dNWV_1e2mu_minus, dNWV_3mu_minus};
  double NWV_plus[4]={dNWV_3e_plus, dNWV_2e1mu_plus, dNWV_1e2mu_plus, dNWV_3mu_plus};

  double sNWV_minus[4]={dsNWV_3e_minus, dsNWV_2e1mu_minus, dsNWV_1e2mu_minus, dsNWV_3mu_minus};
  double sNWV_plus[4]={dsNWV_3e_plus, dsNWV_2e1mu_plus, dsNWV_1e2mu_plus, dsNWV_3mu_plus};
  double NVVV_minus[4]={dNVVV_3e_minus, dNVVV_2e1mu_minus, dNVVV_1e2mu_minus, dNVVV_3mu_minus};
  double NVVV_plus[4]={dNVVV_3e_plus, dNVVV_2e1mu_plus, dNVVV_1e2mu_plus, dNVVV_3mu_plus};

  double sNVVV_minus[4]={dsNVVV_3e_minus, dsNVVV_2e1mu_minus, dsNVVV_1e2mu_minus, dsNVVV_3mu_minus};
  double sNVVV_plus[4]={dsNVVV_3e_plus, dsNVVV_2e1mu_plus, dsNVVV_1e2mu_plus, dsNVVV_3mu_plus};

  double tauFactor_minus[4]={dtauFactor3e_minus, dtauFactor2e1mu_minus, dtauFactor1e2mu_minus, dtauFactor3mu_minus};
  double tauFactor_plus[4]={dtauFactor3e_plus, dtauFactor2e1mu_plus, dtauFactor1e2mu_plus, dtauFactor3mu_plus};
  double stauFactor_minus[4]={dstauFactor3e_minus, dstauFactor2e1mu_minus, dstauFactor1e2mu_minus, dstauFactor3mu_minus};
  double stauFactor_plus[4]={dstauFactor3e_plus, dstauFactor2e1mu_plus, dstauFactor1e2mu_plus, dstauFactor3mu_plus};
  double crossSectionZagreb_minus[4], crossSectionZagreb_plus[4], errorCsZagreb_minus[4], errorCsZagreb_plus[4];
  double ratio_minus_plus[4], ratio_plus_minus[4];
  double csZgNum_minus[4], csZgNum_plus[4], csZgDen_minus[4], csZgDen_plus[4], errorCsZgDen_minus[4], errorCsZgDen_plus[4];

  double lumiSyst[4]={lumiSyst3e/100, lumiSyst2e1mu/100, lumiSyst1e2mu/100, lumiSyst3mu/100};
  double qcdScale[4]={qcdScale3e/100, qcdScale2e1mu/100, qcdScale1e2mu/100, qcdScale3mu/100};
  double PDFsys[4]={PDF3e/100, PDF2e1mu/100, PDF1e2mu/100, PDF3mu/100};
  //  double leptTrgEff[4] ={LepTrgEff3e/100, LepTrgEff2e1mu/100, LepTrgEff1e2mu/100, LepTrgEff3mu/100};
  double leptTrgEff_el[4] ={LepTrgEff3e_el/100, LepTrgEff2e1mu_el/100, LepTrgEff1e2mu_el/100, LepTrgEff3mu_el/100};
  double leptTrgEff_mu[4] ={LepTrgEff3e_mu/100, LepTrgEff2e1mu_mu/100, LepTrgEff1e2mu_mu/100, LepTrgEff3mu_mu/100};
  double Etsys[4] = {Etmiss3e/100, Etmiss2e1mu/100, Etmiss1e2mu/100, Etmiss3mu/100};
  double muMomScale[4] = {muMomScale3e/100, muMomScale2e1mu/100, muMomScale1e2mu/100, muMomScale3mu/100};
  double elEnScale[4]= {elEnScale3e/100, elEnScale2e1mu/100, elEnScale1e2mu/100, elEnScale3mu/100};
  double pileUpsys[4]= {pileUp3e/100, pileUp2e1mu/100, pileUp1e2mu/100, pileUp3mu/100};
  double ZZcs[4]= {ZZcs3e/100, ZZcs2e1mu/100, ZZcs1e2mu/100, ZZcs3mu/100};
  double Zgammacs[4]= {Zgammacs3e/100, Zgammacs2e1mu/100, Zgammacs1e2mu/100, Zgammacs3mu/100};
  double bckgSys[4] = {back3e/100, back2e1mu/100, back1e2mu/100, back3mu/100};
  double dataDrivensys[4]={dataDriven3e/100, dataDriven2e1mu/100, dataDriven1e2mu/100, dataDriven3mu/100};
  double error_qcd[4], error_pdf[4], error_elScale[4], error_muScale[4], error_mc[4], error_leptTrgEff_el[4], error_leptTrgEff_mu[4], error_MET[4], 
    error_pileup[4], error_DD[4], error_lumi[4],  error_stat_relative[4], error_syst[4], error_syst_minus_plus[4], error_syst_plus_minus[4],
    error_stat_relative_plus[4], error_stat_relative_minus[4], error_stat_plus_minus[4], error_stat_minus_plus[4];
  double Z2ll(dZ2ll), W2e(dW2e), W2mu(dW2mu), W2tau(dW2tau);
  // double WZbr[4]={dWZ23e, dWZ22e1mu, dWZ21e2mu, dWZ23mu};
  double WZbr[4]={0.03363*0.1075,0.03363*0.1057,0.03366*0.1075,0.03366*0.1057};
  double Nsig_minus[4], Nsig_plus[4];
  double  systematicErrorZagreb_minus[4]={0,0,0,0};
  double  systematicErrorZagreb_plus[4]={0,0,0,0};
  double sys[4]={totalSys3e/100, totalSys2e1mu/100, totalSys1e2mu/100, totalSys3mu/100};    double luminosity (LUMI);

  double sNZZ_all[4], sNZgamma_all[4], sNWV_all[4], sNVVV_all[4], sN_fake_all[4];
  double sysNZZ[4]={0.16, 0.15, 0.16, 0.15};
  double sysNZgamma[4]={0.15, 0.15, 0.18, 0};
  double sysNWV[4]={0,0,0,0.21};
  double sysNVVV[4]={0.33, 0.32, 0.32, 0.34};
  double sysN_fake[4]={0.607, 0.64, 0.59, 0.6};
  //adding systematic errors                                                                                                                                    
  double error_stat_relative_minus[4]={0, 0, 0, 0};
  double error_stat_relative_plue[4]={0, 0, 0, 0};

  double WZ23lnu=3.*0.033658*(0.1125+0.1075+0.1057);
  TString names[4]={"3e", "2e1mu","1e2mu", "3mu"};

  for (int i=0; i<4; i++){
    Nsig_minus[i]=(N_good_minus[i]-NZgamma_minus[i]-NWV_minus[i]-NVVV_minus[i]-NZZ_minus[i]);
    Nsig_plus[i]=(N_good_plus[i]-NZgamma_plus[i]-NWV_plus[i]-NVVV_plus[i]-NZZ_plus[i]);

    csZgNum_minus[i]=(1-tauFactor_minus[i])*Nsig_minus[i];
    csZgNum_plus[i]=(1-tauFactor_plus[i])*Nsig_plus[i];
    csZgDen_minus[i]=(AxeZ_minus[i]*luminosity*WZbr[i]);
    csZgDen_plus[i]=(AxeZ_plus[i]*luminosity*WZbr[i]);
    crossSectionZagreb_minus[i] =csZgNum_minus[i]/csZgDen_minus[i];
    crossSectionZagreb_plus[i] =csZgNum_plus[i]/csZgDen_plus[i];


    errorCsZagreb_minus[i]=sqrt(pow(((sN_good_minus[i]*(1-tauFactor_minus[i]))/csZgDen_minus[i]),2)+pow((((1-tauFactor_minus[i])*Nsig_minus[i]*luminosity*WZbr[i]*sAxeZ_minus[i])/(csZgDen_minus[i]*csZgDen_minus[i])), 2));
    errorCsZagreb_plus[i]=sqrt(pow(((sN_good_plus[i]*(1-tauFactor_plus[i]))/csZgDen_plus[i]),2)+pow((((1-tauFactor_plus[i])*Nsig_plus[i]*luminosity*WZbr[i]*sAxeZ_plus[i])/(csZgDen_plus[i]*csZgDen_plus[i])), 2));

    
    //systematics
    systematicErrorZagreb_minus[i]=sys[i]*crossSectionZagreb_minus[i];
    systematicErrorZagreb_plus[i]=sys[i]*crossSectionZagreb_plus[i];
    

    error_stat_relative_minus[i]=errorCsZagreb_minus[i]/crossSectionZagreb_minus[i];
    error_stat_relative_plus[i]= errorCsZagreb_plus[i]/crossSectionZagreb_plus[i];

    ratio_minus_plus[i]=crossSectionZagreb_minus[i]/crossSectionZagreb_plus[i];
    ratio_plus_minus[i]=crossSectionZagreb_plus[i]/crossSectionZagreb_minus[i];


    //ratio systematics


    //correlated: qcd, pdf,el. en. scale, mu.mom.scale, MC
    //this doesn't contribute
    error_qcd[i]=0;
    error_pdf[i]=0;
    error_elScale[i]=0;
    error_muScale[i]=0;
    error_mc[i]=0;

    //uncorrelated lept&trigg eff., MET, pileup, data driven, statistics
    error_leptTrgEff_el[i]=sqrt(2*leptTrgEff_el[i]*leptTrgEff_el[i]);
    error_leptTrgEff_mu[i]=sqrt(2*leptTrgEff_mu[i]*leptTrgEff_mu[i]);
    error_MET[i]=sqrt(2*Etsys[i]*Etsys[i]);
    error_pileup[i]=sqrt(2*pileUpsys[i]*pileUpsys[i]);
    error_DD[i]=sqrt(2*dataDrivensys[i]*dataDrivensys[i]);
    error_stat_relative[i]=sqrt(pow(error_stat_relative_minus[i],2)+ pow(error_stat_relative_plus[i], 2));
    
    error_syst[i]=sqrt(error_qcd[i]*error_qcd[i]+error_pdf[i]*error_pdf[i]+error_elScale[i]*error_elScale[i]+
		       error_muScale[i]*error_muScale[i]+error_mc[i]*error_mc[i]+ error_leptTrgEff_el[i]*error_leptTrgEff_el[i]+
		       error_leptTrgEff_mu[i]*error_leptTrgEff_mu[i]+error_MET[i]*error_MET[i]+ error_pileup[i]*error_pileup[i]+ error_DD[i]*error_DD[i]);
    error_syst_minus_plus[i]=error_syst[i]*ratio_minus_plus[i];
    error_syst_plus_minus[i]=error_syst[i]*ratio_plus_minus[i];
    
    error_stat_minus_plus[i]=error_stat_relative[i]*ratio_minus_plus[i];
    error_stat_plus_minus[i]=error_stat_relative[i]*ratio_plus_minus[i];

  }


  ///BLUE METHOD FOR RATIO/////

  Double_t elmZagreb[16], elmZagreb_2[16];
  for (size_t elm=0; elm<16; elm++){
    elmZagreb[elm]=0;
    elmZagreb_2[elm]=0;
  }
  //common elements
  double commonSys[4][4];
  double lumi_matrix[4][4];
  for (int cha=0; cha<4; cha++){
    for (int chb=0; chb<4; chb++){
      commonSys[cha][chb]=error_MET[cha]*error_MET[chb]+ error_pileup[cha]*error_pileup[chb]+
	error_pdf[cha]*error_pdf[chb]+ error_qcd[cha]*error_qcd[chb]+error_mc[cha]*error_mc[chb]+ error_lumi[cha]*error_lumi[chb];
      lumi_matrix[cha][chb]=error_lumi[cha]*error_lumi[chb];
    }
  }

  double corr_matrix_dd[4][4]={
    {1, 0, 1, 0},
    {0, 1, 0, 1},
    {1, 0, 1, 0},
    {0, 1, 0, 1}
  };
  //diagonal elements
  elmZagreb[0]=pow(error_syst_minus_plus[0], 2)+ pow(error_stat_minus_plus[0],2)+pow(error_lumi[0]*ratio_minus_plus[0],2);
  elmZagreb[5]=pow(error_syst_minus_plus[1], 2)+ pow(error_stat_minus_plus[1],2)+pow(error_lumi[1]*ratio_minus_plus[1],2);
  elmZagreb[10]=pow(error_syst_minus_plus[2], 2)+ pow(error_stat_minus_plus[2],2)+pow(error_lumi[2]*ratio_minus_plus[2],2);
  elmZagreb[15]=pow(error_syst_minus_plus[3], 2)+ pow(error_stat_minus_plus[3],2)+pow(error_lumi[3]*ratio_minus_plus[3],2);

  elmZagreb_2[0]=pow(error_syst_plus_minus[0], 2)+ pow(error_stat_plus_minus[0],2)+pow(error_lumi[0]*ratio_plus_minus[0],2);
  elmZagreb_2[5]=pow(error_syst_plus_minus[1], 2)+ pow(error_stat_plus_minus[1],2)+pow(error_lumi[1]*ratio_plus_minus[1],2);
  elmZagreb_2[10]=pow(error_syst_plus_minus[2], 2)+ pow(error_stat_plus_minus[2],2)+pow(error_lumi[2]*ratio_plus_minus[2],2);
  elmZagreb_2[15]=pow(error_syst_plus_minus[3], 2)+ pow(error_stat_plus_minus[3],2)+pow(error_lumi[3]*ratio_plus_minus[3],2);
  

  //matrix is symmetric
  //channels 0 and 1:
  elmZagreb[4]=elmZagreb[1]= ratio_minus_plus[0]*ratio_minus_plus[1]*(commonSys[0][1]+error_elScale[0]*error_elScale[1]+
								      error_muScale[0]*error_muScale[1]+error_leptTrgEff_el[0]*error_leptTrgEff_el[1]+
								      error_leptTrgEff_mu[0]*error_leptTrgEff_mu[1]+
								      error_DD[0]*error_DD[1]*corr_matrix_dd[0][1]);

  elmZagreb_2[4]=elmZagreb_2[1]= ratio_plus_minus[0]*ratio_plus_minus[1]*(commonSys[0][1]+error_elScale[0]*error_elScale[1]+
									  error_muScale[0]*error_muScale[1]+error_leptTrgEff_el[0]*error_leptTrgEff_el[1]+
									  error_leptTrgEff_mu[0]*error_leptTrgEff_mu[1]+
									  error_DD[0]*error_DD[1]*corr_matrix_dd[0][1]);

  //channels 0 and 2:
  elmZagreb[8]=elmZagreb[2]= ratio_minus_plus[0]*ratio_minus_plus[2]*(commonSys[0][2]+ error_elScale[0]*error_elScale[2]+error_muScale[0]*error_muScale[2]+ 
								      error_leptTrgEff_el[0]*error_leptTrgEff_el[2]+ 
								      error_leptTrgEff_mu[0]*error_leptTrgEff_mu[2]+
								      error_DD[0]*error_DD[2]*corr_matrix_dd[0][2]);

  elmZagreb_2[8]=elmZagreb_2[2]=ratio_plus_minus[0]*ratio_plus_minus[2]*(commonSys[0][2]+error_elScale[0]*error_elScale[2]+error_muScale[0]*error_muScale[2]+
									 error_leptTrgEff_el[0]*error_leptTrgEff_el[2]+
									 error_leptTrgEff_mu[0]*error_leptTrgEff_mu[2]+
									 error_DD[0]*error_DD[2]*corr_matrix_dd[0][2]);

  //channels 0 and 3
  elmZagreb[12]=elmZagreb[3]= ratio_minus_plus[0]*ratio_minus_plus[3]*(commonSys[0][3]+ error_elScale[0]*error_elScale[3]+ error_muScale[0]*error_muScale[3]+
								       error_leptTrgEff_el[0]*error_leptTrgEff_el[3]+
								       error_leptTrgEff_mu[0]*error_leptTrgEff_mu[3]+
								       error_DD[0]*error_DD[3]*corr_matrix_dd[0][3]);

  elmZagreb_2[12]=elmZagreb_2[3]= ratio_plus_minus[0]*ratio_plus_minus[3]*(commonSys[0][3]+ error_elScale[0]*error_elScale[3]+ error_muScale[0]*error_muScale[3]+
									   error_leptTrgEff_el[0]*error_leptTrgEff_el[3]+
									   error_leptTrgEff_mu[0]*error_leptTrgEff_mu[3]+
									   error_DD[0]*error_DD[3]*corr_matrix_dd[0][3]);

  //channels 1 and 2 
  elmZagreb[9]=elmZagreb[6]= ratio_minus_plus[1]*ratio_minus_plus[2]*(commonSys[1][2]+ error_elScale[1]*error_elScale[2]+ error_muScale[1]*error_muScale[2]+
								      error_leptTrgEff_el[1]*error_leptTrgEff_el[2]+
								      error_leptTrgEff_mu[1]*error_leptTrgEff_mu[2]+
								      error_DD[1]*error_DD[2]*corr_matrix_dd[1][2]);

  elmZagreb_2[9]=elmZagreb_2[6]= ratio_plus_minus[1]*ratio_plus_minus[2]*(commonSys[1][2]+ error_elScale[1]*error_elScale[2]+ error_muScale[1]*error_muScale[2]+
									  error_leptTrgEff_el[1]*error_leptTrgEff_el[2]+
									  error_leptTrgEff_mu[1]*error_leptTrgEff_mu[2]+
									  error_DD[1]*error_DD[2]*corr_matrix_dd[1][2]);
  //channels 1 and 3
  elmZagreb[13]=elmZagreb[7]=ratio_minus_plus[1]*ratio_minus_plus[3]*(commonSys[1][3]+error_elScale[1]*error_elScale[3]+error_muScale[1]*error_muScale[3]+
								      error_leptTrgEff_el[1]*error_leptTrgEff_el[3]+
								      error_leptTrgEff_mu[1]*error_leptTrgEff_mu[3]+
								      error_DD[1]*error_DD[3]*corr_matrix_dd[1][3]);

elmZagreb_2[13]=elmZagreb_2[7]=ratio_plus_minus[1]*ratio_plus_minus[3]*(commonSys[1][3]+error_elScale[1]*error_elScale[3]+error_muScale[1]*error_muScale[3]+
								        error_leptTrgEff_el[1]*error_leptTrgEff_el[2]+
									error_leptTrgEff_mu[1]*error_leptTrgEff_mu[2]+
									error_DD[1]*error_DD[2]*corr_matrix_dd[1][2]);

  //channels 2 and 3
  elmZagreb[11]=elmZagreb[14]=ratio_minus_plus[2]*ratio_minus_plus[3]*(commonSys[2][3]+error_elScale[2]*error_elScale[3]+error_muScale[2]*error_muScale[3]+
								       error_leptTrgEff_el[2]*error_leptTrgEff_el[3]+
								       error_leptTrgEff_mu[2]*error_leptTrgEff_mu[3]+
								       error_DD[2]*error_DD[3]*corr_matrix_dd[2][3]);

  elmZagreb_2[11]=elmZagreb_2[14]=ratio_plus_minus[2]*ratio_plus_minus[3]*(commonSys[2][3]+error_elScale[2]*error_elScale[3]+error_muScale[2]*error_muScale[3]+
									   error_leptTrgEff_el[2]*error_leptTrgEff_el[3]+
									   error_leptTrgEff_mu[2]*error_leptTrgEff_mu[3]+
									   error_DD[2]*error_DD[3]*corr_matrix_dd[2][3]);;


  std::cout<<"Ratio minus/plus error matrix finished!!"<<std::endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
       std:: cout << "       " << elmZagreb[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl; 

  std::cout<<"Ratio plus/minus error matrix finished!!"<<std::endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
       std:: cout << "       " << elmZagreb_2[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl; 

  //inititalize matrix for inversion
  TMatrixD errMat(4,4,elmZagreb);
  TMatrixD errMatInv(4,4);

  TMatrixD errMatCopy(errMat);

  TMatrixD errMat_2(4,4,elmZagreb_2);
  TMatrixD errMatInv_2(4,4);

  TMatrixD errMatCopy_2(errMat_2);

  //invert 
  errMatInv=errMat.Invert();
  errMatInv_2=errMat_2.Invert();
  
  Double_t *mRefTest_2=errMat_2.GetMatrixArray();
  Double_t *mRefTest=errMat.GetMatrixArray();
  
  std::cout << endl << "Ratio minus/plus matrix inverse:" << endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      std::cout << "       " << mRefTest[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl;


  std::cout << endl << "Ratio plus/minus matrix inverse:" << endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      std::cout << "       " << mRefTest_2[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl;

  //get norm and alpha factors for each channel
  Double_t *mRef= errMat.GetMatrixArray();
  Double_t norm=0.;
  Double_t alphaCH[4]={0.,0.,0.,0.};

  Double_t *mRef_2= errMat_2.GetMatrixArray();
  Double_t norm_2=0.;
  Double_t alphaCH_2[4]={0.,0.,0.,0.};

  for (int imatrix=0;imatrix<16;imatrix++) norm+=mRef[imatrix];
  for (int imatrix=0;imatrix<16;imatrix++) norm_2+=mRef_2[imatrix];

  for (size_t im=0;im<4;im++) {
    for (size_t jm=0;jm<4;jm++) {
      alphaCH[im]+=mRef[im*4+jm];
    }
    alphaCH[im]/=norm;
  }
  
  for (size_t im=0;im<4;im++) {
    for (size_t jm=0;jm<4;jm++) {
      alphaCH_2[im]+=mRef_2[im*4+jm];
    }
    alphaCH_2[im]/=norm_2;
  }

  std::cout<<"RATIO MINUS/PLUS"<<std::endl;
  std::cout << "al0 " << alphaCH[0] << " al1 " << alphaCH[1] << " al2 " << alphaCH[2] << " al3 " << alphaCH[3] << endl;
  std::cout << "consistency check:" << alphaCH[0]+alphaCH[1]+alphaCH[2]+alphaCH[3] <<endl;
  std::cout << endl;

  std::cout<<"RATIO PLUS/MINUS"<<std::endl;
  std::cout << "al0 " << alphaCH_2[0] << " al1 " << alphaCH_2[1] << " al2 " << alphaCH_2[2] << " al3 " << alphaCH_2[3] << endl;
  std::cout << "consistency check:" << alphaCH_2[0]+alphaCH_2[1]+alphaCH_2[2]+alphaCH_2[3] <<endl;
  std::cout << endl;
  //  std::cout << "al0_plus " << alphaCH_plus[0] << " al1_plus " << alphaCH_plus[1] << " al2_plus " << alphaCH_plus[2] << " al3_plus " << alphaCH_plus[3] << endl;
  //std::cout << "consistency check:" << alphaCH_plus[0]+alphaCH_plus[1]+alphaCH_plus[2]+alphaCH_plus[3] <<endl;

  double final_ratio_minus_plus = alphaCH[0]*ratio_minus_plus[0] + alphaCH[1]*ratio_minus_plus[1] + alphaCH[2]*ratio_minus_plus[2] + alphaCH[3]*ratio_minus_plus[3];
  double final_ratio_plus_minus = alphaCH_2[0]*ratio_plus_minus[0] + alphaCH_2[1]*ratio_plus_minus[1] + alphaCH_2[2]*ratio_plus_minus[2] + alphaCH_2[3]*ratio_plus_minus[3];

  Double_t combined_error_minus_plus=0;
  Double_t combined_error_plus_minus=0;
  Double_t *copyRef_minus_plus = errMatCopy.GetMatrixArray();
  Double_t *copyRef_plus_minus = errMatCopy_2.GetMatrixArray();
  for (int ier=0;ier<4;ier++)
    for (int jer=0;jer<4;jer++) {
      combined_error_minus_plus+=alphaCH[ier]*alphaCH[jer]*copyRef_minus_plus[ier*4+jer];
      combined_error_plus_minus+=alphaCH_2[ier]*alphaCH_2[jer]*copyRef_plus_minus[ier*4+jer];
    }  
  Double_t stat_err_tot_minus_plus2 =  pow(alphaCH[0]*error_stat_minus_plus[0],2) + pow(alphaCH[1]*error_stat_minus_plus[1],2)
    +pow(alphaCH[2]*error_stat_minus_plus[2],2) + pow(alphaCH[3]*error_stat_minus_plus[3],2);

  Double_t stat_err_tot_minus_plus = sqrt(stat_err_tot_minus_plus2);
  Double_t syst_err_tot_minus_plus = sqrt (combined_error_minus_plus - stat_err_tot_minus_plus2);
  std::cout<<"Stat error squared: "<<stat_err_tot_minus_plus2<<std::endl;
  std::cout<<"Sta error: "<<stat_err_tot_minus_plus<<std::endl;

  combined_error_plus_minus=sqrt(combined_error_minus_plus);
  std::cout<<"Stat Error:   "<<stat_err_tot_minus_plus<<std::endl;
  std::cout<<"Syst Error:   "<<syst_err_tot_minus_plus<<std::endl;
  std::cout<<"Total Error: "<<combined_error_minus_plus<<std::endl;

  Double_t stat_err_tot_plus_minus2 =  pow(alphaCH_2[0]*error_stat_plus_minus[0],2) + pow(alphaCH_2[1]*error_stat_plus_minus[1],2)
    +pow(alphaCH_2[2]*error_stat_plus_minus[2],2) + pow(alphaCH_2[3]*error_stat_plus_minus[3],2);
  Double_t stat_err_tot_plus_minus = sqrt(stat_err_tot_plus_minus2);
  Double_t syst_err_tot_plus_minus = sqrt (combined_error_plus_minus - stat_err_tot_plus_minus2);
  
  
  combined_error_plus_minus=sqrt(combined_error_plus_minus);

  std::cout << "combined ratio(W-Z)/(W+Z) = "<<fixed<<setprecision(3) << final_ratio_minus_plus << " $\\pm$ " << stat_err_tot_minus_plus << "(stat.)$\\pm$ " << syst_err_tot_minus_plus << "(syst) " << endl;
  std::cout << "combined ratio(W+Z)/(W-Z) = "<<fixed<<setprecision(3) << final_ratio_plus_minus << " $\\pm$ " << stat_err_tot_plus_minus << "(stat.)$\\pm$ " << syst_err_tot_plus_minus << "(syst) " << endl;
  
  std::cout<<"CONSISTENCY CHECK:"<<std::endl;
  for (int check=0; check<4; check++){
    std::cout<<crossSectionZagreb_minus[check]+crossSectionZagreb_plus[check]<<std::endl;
  }

  if (latexOutput){
    std::cout<<"latex output"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"crossSection"<<std::endl;
    std::cout<<"\\begin{center}"<<std::endl;
    std::cout<<"\\begin{tabular}{c|c|c}"<<std::endl;
    std::cout<<"channel &  crossSection W-Z & crossSection W+Z\\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"  & $"<<crossSectionZagreb_minus[a]<<"\\pm "<<errorCsZagreb_minus[a]<<"\\pm "<<systematicErrorZagreb_minus[a]<<"$  & $"<<crossSectionZagreb_plus[a]<<"\\pm "<<errorCsZagreb_plus[a]<<"\\pm"<<systematicErrorZagreb_plus[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
    
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{center}"<<std::endl;
    std::cout<<"*****************"<<std::endl;

    std::cout<<"*****************"<<std::endl;
    std::cout<<"crossSection"<<std::endl;
    std::cout<<"\\begin{center}"<<std::endl;
    std::cout<<"\\begin{tabular}{c|c|c}"<<std::endl;
    std::cout<<"channel &  ratio W-Z/W+Z & ratio W+Z/W-Z\\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"  & $"<<fixed<<setprecision(2)<<ratio_minus_plus[a]<<"\\pm "<<error_syst_minus_plus[a]<<"\\pm "<<error_stat_minus_plus[a]<<"$  & $"<<ratio_plus_minus[a]<<"\\pm "<<error_syst_plus_minus[a]<<"\\pm"<<error_stat_plus_minus[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
    
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{center}"<<std::endl;
    std::cout<<"*****************"<<std::endl;

  }

  /*
/////BLUE METHOD////////
  //TO MAKNUTI!!!
  
  //compose Error Matrix [i = column*4+row] [0-3][0-3] indices
//matrix is symmetric
//channel legend: 0 = 3e, 1 = 2e1mu, 2 = 2mu1e,  3 = 3mu
  

  Double_t elmZagreb_minus[16], elmZagreb_plus[16];
  for (size_t elm=0;elm<16;elm++){
    elmZagreb_minus[elm]=0;
    elmZagreb_plus[elm]=0;
  }
  
  //common elements
  // [Et_miss, pileUp, PDF, QCDscales, backg]
  double commonSys[4][4];
  for (int cha=0; cha<4; cha++){
    for (int chb=0; chb<4; chb++){
      commonSys[cha][chb]= Etsys[cha]*Etsys[chb]+ pileUpsys[cha]*pileUpsys[chb]
	+ PDFsys[cha]*PDFsys[chb] + qcdScale[cha]*qcdScale[chb] +
	bckgSys[cha]*bckgSys[chb];
    }
  }
  //diagonals elements:
  elmZagreb_minus[0]=pow(systematicErrorZagreb_minus[0],2) + pow(errorCsZagreb_minus[0],2);
  elmZagreb_minus[5]=pow(systematicErrorZagreb_minus[1],2) + pow(errorCsZagreb_minus[1],2);
  elmZagreb_minus[10]=pow(systematicErrorZagreb_minus[2],2) + pow(errorCsZagreb_minus[2],2);
  elmZagreb_minus[15]=pow(systematicErrorZagreb_minus[3],2) + pow(errorCsZagreb_minus[3],2);

  elmZagreb_plus[0]=pow(systematicErrorZagreb_plus[0],2) + pow(errorCsZagreb_plus[0],2);
  elmZagreb_plus[5]=pow(systematicErrorZagreb_plus[1],2) + pow(errorCsZagreb_plus[1],2);
  elmZagreb_plus[10]=pow(systematicErrorZagreb_plus[2],2) + pow(errorCsZagreb_plus[2],2);
  elmZagreb_plus[15]=pow(systematicErrorZagreb_plus[3],2) + pow(errorCsZagreb_plus[3],2);
  
  //matrix is symmetric
  //channels 0 and 1 : [ trigger&scaleFactors, electronEnergyScale, muonEnergyScale]
  elmZagreb_minus[4]=elmZagreb_minus[1]= crossSectionZagreb_minus[0]*crossSectionZagreb_minus[1]* (commonSys[0][1] + elEnScale[0]*elEnScale[1] + 
							 muMomScale[0]*muMomScale[1]+leptTrgEff[0]*(sqrt(2/3)*leptTrgEff[1]));
  elmZagreb_plus[4]=elmZagreb_plus[1]= crossSectionZagreb_plus[0]*crossSectionZagreb_plus[1]* (commonSys[0][1] + elEnScale[0]*elEnScale[1] + 
                                                         muMomScale[0]*muMomScale[1]+leptTrgEff[0]*(sqrt(2/3)*leptTrgEff[1]));
  //channels 0 and 2: 
  elmZagreb_minus[8]=elmZagreb_minus[2]=  crossSectionZagreb_minus[0]*crossSectionZagreb_minus[2]* (commonSys[0][2] + elEnScale[0]*elEnScale[2] + 
                                                        muMomScale[0]*muMomScale[2]+leptTrgEff[0]*(sqrt(1/3)*leptTrgEff[2]));
  elmZagreb_plus[8]=elmZagreb_plus[2]=  crossSectionZagreb_plus[0]*crossSectionZagreb_plus[2]* (commonSys[0][2] + elEnScale[0]*elEnScale[2] + 
                                                       muMomScale[0]*muMomScale[2]+leptTrgEff[0]*(sqrt(1/3)*leptTrgEff[2]));

  //channels 0 and 3
  elmZagreb_minus[12]=elmZagreb_minus[3]= crossSectionZagreb_minus[0]*crossSectionZagreb_minus[3]* (commonSys[0][3] + elEnScale[0]*elEnScale[3] +
                                             muMomScale[0]*muMomScale[3]);
  elmZagreb_plus[12]=elmZagreb_plus[3]= crossSectionZagreb_plus[0]*crossSectionZagreb_plus[3]* (commonSys[0][3] + elEnScale[0]*elEnScale[3] + 
                                             muMomScale[0]*muMomScale[3]);

  //channels 1 and 2
  elmZagreb_minus[9]=elmZagreb_minus[6] =crossSectionZagreb_minus[1]*crossSectionZagreb_minus[2]* (commonSys[1][2] +elEnScale[1]*elEnScale[2] + 
                                                muMomScale[1]*muMomScale[2]+2*sqrt(1/3)*leptTrgEff[1]*sqrt(2/3)*leptTrgEff[2]);
  elmZagreb_plus[9]=elmZagreb_plus[6] =crossSectionZagreb_plus[1]*crossSectionZagreb_plus[2]* (commonSys[1][2] +elEnScale[1]*elEnScale[2] +
                                                muMomScale[1]*muMomScale[2]+ 2*sqrt(1/3)*leptTrgEff[1]*sqrt(2/3)*leptTrgEff[2]);
  //channels 1 and 3
  elmZagreb_minus[13]=elmZagreb_minus[7] =crossSectionZagreb_minus[1]*crossSectionZagreb_minus[3] *(commonSys[1][3] + elEnScale[1]*elEnScale[3] + 
                                                muMomScale[1]*muMomScale[3]+leptTrgEff[1]*(sqrt(1/3)*leptTrgEff[3]));
  elmZagreb_plus[13]=elmZagreb_plus[7] =crossSectionZagreb_plus[1]*crossSectionZagreb_plus[3] *(commonSys[1][3] + elEnScale[1]*elEnScale[3] + 
                                                muMomScale[1]*muMomScale[3]+leptTrgEff[1]*(sqrt(1/3)*leptTrgEff[3]));
  //channels 2 and 3
  elmZagreb_minus[11]=elmZagreb_minus[14] =crossSectionZagreb_minus[2]*crossSectionZagreb_minus[3] *(commonSys[2][3] + elEnScale[2]*elEnScale[3] + 
                                                muMomScale[2]*muMomScale[3]+leptTrgEff[2]*(sqrt(2/3)*leptTrgEff[3]));
  elmZagreb_plus[11]=elmZagreb_plus[14] =crossSectionZagreb_plus[2]*crossSectionZagreb_plus[3] *(commonSys[2][3] + elEnScale[2]*elEnScale[3] + 
                                                muMomScale[2]*muMomScale[3]+leptTrgEff[2]*(sqrt(2/3)*leptTrgEff[3]));
 
  std::cout<<"W- error matrix finished!!"<<std::endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      //      std::cout<<"!!"<<std::endl;
      std:: cout << "       " << elmZagreb_minus[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl; 


  std::cout<<"W+ error matrix finished!!"<<std::endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      //      std::cout<<"!!"<<std::endl;
      std:: cout << "       " << elmZagreb_plus[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl; 
  
//initialize matrix for inversion
  TMatrixD errMat_minus(4,4,elmZagreb_minus);
  TMatrixD errMatInv_minus(4,4);
  
  TMatrixD errMatCopy_minus(errMat_minus);

  TMatrixD errMat_plus(4,4,elmZagreb_plus);
  TMatrixD errMatInv_plus(4,4);
  
  TMatrixD errMatCopy_plus(errMat_plus);
  
  //invert !
  errMatInv_minus = errMat_minus.Invert();
  errMatInv_plus = errMat_plus.Invert();
  
  
  Double_t *mRefTest_minus = errMat_minus.GetMatrixArray();
  Double_t *mRefTest_plus = errMat_plus.GetMatrixArray();
  
  std::cout << endl << "W- matrix Inverse:" << endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      std::cout << "       " << mRefTest_minus[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl;

  std::cout << endl << "W+ matrix Inverse:" << endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      std::cout << "       " << mRefTest_plus[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl;


  //get norm, and alpha factors for each channel
  Double_t *mRef_minus= errMat_minus.GetMatrixArray();
  Double_t norm_minus=0.;
  Double_t alphaCH_minus[4]={0.,0.,0.,0.};
  Double_t *mRef_plus= errMat_plus.GetMatrixArray();
  Double_t norm_plus=0.;
  Double_t alphaCH_plus[4]={0.,0.,0.,0.};

  for (int imatrix=0;imatrix<16;imatrix++) norm_minus+=mRef_minus[imatrix];
  for (int imatrix=0;imatrix<16;imatrix++) norm_plus+=mRef_plus[imatrix];
  
  for (size_t im=0;im<4;im++) {
    for (size_t jm=0;jm<4;jm++) {
      alphaCH_minus[im]+=mRef_minus[im*4+jm];
    }
    alphaCH_minus[im]/=norm_minus;
  }

  for (size_t im=0;im<4;im++) {
    for (size_t jm=0;jm<4;jm++) {
      alphaCH_plus[im]+=mRef_plus[im*4+jm];
    }
    alphaCH_plus[im]/=norm_plus;
  }

  std::cout << "al0_minus " << alphaCH_minus[0] << " al1_minuns " << alphaCH_minus[1] << " al2_minus " << alphaCH_minus[2] << " al3_minus " << alphaCH_minus[3] << endl;
  std::cout << "consistency check:" << alphaCH_minus[0]+alphaCH_minus[1]+alphaCH_minus[2]+alphaCH_minus[3] <<endl;
  std::cout << endl;
  std::cout << "al0_plus " << alphaCH_plus[0] << " al1_plus " << alphaCH_plus[1] << " al2_plus " << alphaCH_plus[2] << " al3_plus " << alphaCH_plus[3] << endl;
  std::cout << "consistency check:" << alphaCH_plus[0]+alphaCH_plus[1]+alphaCH_plus[2]+alphaCH_plus[3] <<endl;
  std::cout << endl;
  double final_Xsec_minus = alphaCH_minus[0]*crossSectionZagreb_minus[0] + alphaCH_minus[1]*crossSectionZagreb_minus[1] + alphaCH_minus[2]*crossSectionZagreb_minus[2] + alphaCH_minus[3]*crossSectionZagreb_minus[3];
  double final_Xsec_plus = alphaCH_plus[0]*crossSectionZagreb_plus[0] + alphaCH_plus[1]*crossSectionZagreb_plus[1] + alphaCH_plus[2]*crossSectionZagreb_plus[2] + alphaCH_plus[3]*crossSectionZagreb_plus[3];

  Double_t combined_error_minus=0;
  Double_t combined_error_plus=0;
  Double_t *copyRef_minus = errMatCopy_minus.GetMatrixArray();
  Double_t *copyRef_plus = errMatCopy_plus.GetMatrixArray();
  for (int ier=0;ier<4;ier++)
    for (int jer=0;jer<4;jer++) {
      combined_error_minus+=alphaCH_minus[ier]*alphaCH_minus[jer]*copyRef_minus[ier*4+jer];
      combined_error_plus+=alphaCH_plus[ier]*alphaCH_plus[jer]*copyRef_plus[ier*4+jer];
    }
  Double_t stat_err_tot2_minus =  pow(alphaCH_minus[0]*errorCsZagreb_minus[0],2) + pow(alphaCH_minus[1]*errorCsZagreb_minus[1],2)
    +pow(alphaCH_minus[2]*errorCsZagreb_minus[2],2) + pow(alphaCH_minus[3]*errorCsZagreb_minus[3],2);
  Double_t stat_err_tot2_plus =  pow(alphaCH_plus[0]*errorCsZagreb_plus[0],2) + pow(alphaCH_plus[1]*errorCsZagreb_plus[1],2)
    +pow(alphaCH_plus[2]*errorCsZagreb_plus[2],2) + pow(alphaCH_plus[3]*errorCsZagreb_plus[3],2);
  
  Double_t stat_err_tot_minus = sqrt(stat_err_tot2_minus);
  Double_t syst_err_tot_minus = sqrt (combined_error_minus - stat_err_tot2_minus);
  Double_t stat_err_tot_plus = sqrt(stat_err_tot2_plus);
  Double_t syst_err_tot_plus = sqrt (combined_error_plus - stat_err_tot2_plus);
  
  combined_error_minus=sqrt(combined_error_minus);
  combined_error_plus=sqrt(combined_error_plus);


  std::cout << "combined sigma(W-Z) = " << final_Xsec_minus << " +- " << stat_err_tot_minus << "(stat.) +- " << syst_err_tot_minus << "(syst) +- "
	    << 0.026*final_Xsec_minus << " (lumi) pb " << endl;
  std::cout << "combined sigma(W+Z) = " << final_Xsec_plus << " +- " << stat_err_tot_plus << "(stat.) +- " << syst_err_tot_plus << "(syst) +- "
	    << 0.026*final_Xsec_plus << " (lumi) pb " << endl;
  
  std::cout << endl;
 //systematics numbers
 

  ////**********outputs**********************
  std::cout<<"*****************************************"<<std::endl;
  std::cout<<"*********W+Z cross section(inclusive)************"<<std::endl;
  for (int i2=0; i2<4; i2++){
    std::cout<<i2<<" : "<<crossSectionZagreb_plus[i2] <<"+/-"<<errorCsZagreb_plus[i2] <<"+/-"<<systematicErrorZagreb_plus[i2]<<std::endl;
  }
  std::cout<<"*********W-Z cross section(inclusive)************"<<std::endl;
  for (int i3=0; i3<4; i3++){
    std::cout<<i3<<" : "<<crossSectionZagreb_minus[i3] <<"+/-"<<errorCsZagreb_minus[i3] <<"+/-"<<systematicErrorZagreb_minus[i3]<<std::endl;
  }
  */ 
   ////**************latexOutput
  

  /*
  std::cout<<"Matrix method results"<<std::endl;
  std::cout<<"\\begin{center}"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c|c|c}"<<std::endl;
  std::cout<<"channel & Ndata & Ngood & Nfake\\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  for (int mm=0; mm<4; mm++){
    std::cout<<names[mm]<<" & " <<N_data[mm]<<" &  $"<<N_good[mm]<<"\\pm "<<sN_good[mm] <<"$ & $"<<N_fake[mm]<<"\\pm"<<sN_fake[mm]<<"$ \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
  }
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{center}"<<std::endl;

  std::cout<<"*****************"<<std::endl;
  std::cout<<"Matrix method results WITH SYSTEMATICS"<<std::endl;
  std::cout<<"\\begin{center}"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c|c|c}"<<std::endl;
  std::cout<<"channel & Ndata  & Nfake\\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  for (int mm=0; mm<4; mm++){
    std::cout<<names[mm]<<" & " <<N_data[mm]<<"$ & $"<<N_fake[mm]<<"\\pm"<<sN_fake_all[mm]<<"$ \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
  }
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{center}"<<std::endl;

  std::cout<<"*****************"<<std::endl;
  std::cout<<"MC  yields"<<std::endl;
  std::cout<<"\\begin{center}"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c|c|c|c}"<<std::endl;
  std::cout<<" & 3e & 2e1mu & 1e2mu & 3mu \\\\"<<std::endl;
  std::cout<<"ZZ & $"<<NZZ[0]<<"\\pm "<<sNZZ[0]<<"$ & $"<<NZZ[1]<<"\\pm "<<sNZZ[1]<<"$ & $"<<NZZ[2]<<"\\pm "<<sNZZ[2]<<"$ & $"<<NZZ[3]<<"\\pm "<<sNZZ[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"Zgamma & $"<<NZgamma[0]<<"\\pm "<<sNZgamma[0]<<"$ & $"<<NZgamma[1]<<"\\pm "<<sNZgamma[1]<<"$ & $"<<NZgamma[2]<<"\\pm "<<sNZgamma[2]<<"$ & $"<<NZgamma[3]<<"\\pm "<<sNZgamma[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"WV & $"<<NWV[0]<<"\\pm "<<sNWV[0]<<"$ & $"<<NWV[1]<<"\\pm "<<sNWV[1]<<"$ & $"<<NWV[2]<<"\\pm "<<sNWV[2]<<"$ & $"<<NWV[3]<<"\\pm "<<sNWV[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"VVV & $"<<NVVV[0]<<"\\pm "<<sNVVV[0]<<"$ & $"<<NVVV[1]<<"\\pm "<<sNVVV[1]<<"$ & $"<<NVVV[2]<<"\\pm "<<sNVVV[2]<<"$ & $"<<NVVV[3]<<"\\pm "<<sNVVV[3]<<"$ \\\\"<<std::endl;
  
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{center}"<<std::endl;



  std::cout<<"*****************"<<std::endl;
  std::cout<<"MC  yields WITH SYSTEMATICS"<<std::endl;
  std::cout<<"\\begin{center}"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c|c|c|c}"<<std::endl;
  std::cout<<" & 3e & 2e1mu & 1e2mu & 3mu \\\\"<<std::endl;
  std::cout<<"ZZ & $"<<NZZ[0]<<"\\pm "<<sNZZ_all[0]<<"$ & $"<<NZZ[1]<<"\\pm "<<sNZZ_all[1]<<"$ & $"<<NZZ[2]<<"\\pm "<<sNZZ_all[2]<<"$ & $"<<NZZ[3]<<"\\pm "<<sNZZ_all[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"Zgamma & $"<<NZgamma[0]<<"\\pm "<<sNZgamma_all[0]<<"$ & $"<<NZgamma[1]<<"\\pm "<<sNZgamma_all[1]<<"$ & $"<<NZgamma[2]<<"\\pm "<<sNZgamma_all[2]<<"$ & $"<<NZgamma[3]<<"\\pm "<<sNZgamma_all[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"WV & $"<<NWV[0]<<"\\pm "<<sNWV_all[0]<<"$ & $"<<NWV[1]<<"\\pm "<<sNWV_all[1]<<"$ & $"<<NWV[2]<<"\\pm "<<sNWV_all[2]<<"$ & $"<<NWV[3]<<"\\pm "<<sNWV_all[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"VVV & $"<<NVVV[0]<<"\\pm "<<sNVVV_all[0]<<"$ & $"<<NVVV[1]<<"\\pm "<<sNVVV_all[1]<<"$ & $"<<NVVV[2]<<"\\pm "<<sNVVV_all[2]<<"$ & $"<<NVVV[3]<<"\\pm "<<sNVVV_all[3]<<"$ \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{center}"<<std::endl;
  //    std::cout<<"WZ & $"<<NWZ[0]<<"\\pm "<<sNWZ[0]<<"$ &"<<NWZ[1]<<"\\pm "<<sNWZ[1]<<"$ &"<<NWZ[2]<<"\\pm "<<sNWZ[2]<<"$ &"<<NWZ[3]<<"\\pm "<<sNWZ[3]<<"$ \\\\"<<std::endl;
  
  std::cout<<"Acceptance times efficiency"<<std::endl;
  std::cout<<"\\begin{center}"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c|c|c}"<<std::endl;    
  std::cout<<"channel & AxeSpanish  & AxeZagreb \\\\"<<std::endl;
  for (int a=0; a<4; a++){
    std::cout<<names[a]<<fixed<<setprecision(4)<<"  & $"<<AxeS[a]<<"\\pm "<<sAxeS[a]<<"$ & $"<<AxeZ[a]<<"\\pm "<<sAxeZ[a]<<"$\\\\"<<std::endl;  
    std::cout<<"\\hline"<<std::endl;
  }
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{center}"<<std::endl;
  std::cout<<"*****************"<<std::endl;
    */
  /*
    
    std::cout<<"Matrix method"<<std::endl;
    std::cout<<"channel & Ngood  & AxeNfake"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<a<<"  & $"<<AxeS[a]<<"\pm "<<sAxeS[a]<<"$ & $"<<AxeZ[a]<<"/pm "<<sAxeZ[a]<<"\\\\"<<std::endl;  
      std::cout<<"\\hline"<<std::endl;
    }
  */
  /*
  std::cout<<"*****************"<<std::endl;
  std::cout<<"*****************"<<std::endl;
  std::cout<<"tau factor"<<std::endl;
  std::cout<<"\\begin{center}"<<std::endl;
  std::cout<<"\\begin{tabular}{c|c}"<<std::endl;
  std::cout<<"channel &  tauFactor\\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"  & $"<<tauFactor[a]<<"\\pm "<<stauFactor[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{center}"<<std::endl;
    std::cout<<"*****************"<<std::endl;

  if (outputNumbers){
    for (int cross=0; cross<4;cross++){
      fileNum<<crossSectionZagreb[cross]<<std::endl;
    }
    for (int tau=0; tau<4; tau++){
      fileNum<<tauFactor[tau]<<std::endl;
    }
    for (int axe=0; axe<4; axe++){
      fileNum<<AxeZ[axe]<<std::endl;
    }
  }
*/

  }
