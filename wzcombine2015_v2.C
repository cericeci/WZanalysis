{
#include <algorithm>
//#include <iostream>
cout << std::endl;

TMatrixD errMat(4,4);
//#include <iomanip>

//corr-factors {1.00352, 0.992499, 1.00634, 0.987805};

Double_t rhoEle_ScaleErrUp0[]={0.0134,0.0070,0.0086,0.0000};
Double_t rhoEle_ScaleErrDown0[]={0.0190,0.0065,0.0114,0.0000};

Double_t rhoMu_ScaleErrUp0[] ={0.0000,0.0057,0.0041,0.0102};
Double_t rhoMu_ScaleErrDown0[] ={0.0000,0.0047,0.0037,0.0085};

Double_t rhoMET_Err0[4] = {0.0313,0.0286,0.0341,0.0317};
Double_t rhoPDF_Err0[4] = {0.0140,0.0140,0.0140,0.0140};
Double_t rhoQFScale_Err0[4] = {0.0130,0.0130,0.0130,0.0130};

//PU (up,down)
Double_t PURW_ErrUp0[4] = {0.00016,0.0024,0.0033,0.0039};
Double_t PURW_ErrDown0[4] = {0.00126,0.0040,0.0051,0.0053};


/********************************************************************************************/
Double_t rhoEle_ScaleErr[]={0.,0.,0.,0.0}; //electron energy scale
Double_t rhoEle_ScaleErrUp[]={0.0,0.0,0.0,0.0};
Double_t rhoEle_ScaleErrDown[]={0.0,0.0,0.0,0.0};

Double_t rhoMu_ScaleErr[] ={0.0,0.0,0.0,0.0};  //muon momentum scale
Double_t rhoMu_ScaleErrUp[] ={0.0,0.0,0.0,0.0};
Double_t rhoMu_ScaleErrDown[] ={0.0,0.0,0.0,0.0};

Double_t rhoMET_Err[4] = {0.0,0.0,0.0,0.0}; //MET
Double_t rhoPDF_Err[4] = {0.0,0.0,0.0,0.0}; //PDF syst

//PU (up,down)
Double_t PURW_Err[4] = {0.00,0.0,0.0,0.0};
Double_t PURW_ErrUp[4] = {0.0,0.0,0.0,0.0};
Double_t PURW_ErrDown[4] = {0.0,0.0,0.0,0.0};

/********************************************************************************************/
/* ACCEPTANCE x EFFICIENCY */
//Rescaled acceptance
Double_t acceptance[4]=         { 0.136247, 0.152918,0.170436,0.226478 };
Double_t acceptance_statErr[4]= { 0.00117, 0.00124,0.00128,0.00144 };
Double_t acceptance_systErr0[4]= { 0.00378, 0.00373,0.00314,0.00294 };//SF and trigger
Double_t acceptance_systErr[4];
Double_t acceptance_Err[4];

Double_t acceptance_method2[4]=         {0.0160682, 0.0178895, 0.0201011, 0.0265274 };

Double_t intLumi=4922.;
Double_t intLumiErr = 0.022*intLumi; //2.2%

//YIELD from Matrix Method
Double_t dataYield[4]={61.83,60.48,67.59,95.23};
Double_t dataYield_Err[4]={8.08,7.93,8.42,9.89};

Double_t dataYield_forSyst[4]={60.32,57.06,63.78,90.23};
//MM syst error
//old
//Double_t dataYield_systErr[4]={0.024*dataYield[0],0.009*dataYield[1],0.015*dataYield[2],0.009*dataYield[3]};
//new
Double_t dataYield_systErr[4];
for (int kp=0;kp<4;kp++) {dataYield_systErr[kp] = fabs(dataYield[kp]-dataYield_forSyst[kp]); std::cout<<"dy " << fabs(dataYield[kp]-dataYield_forSyst[kp])/dataYield[kp] << endl;}
//Double_t dataYield_systErr[4]={0.024*dataYield[0],0.009*dataYield[1],0.015*dataYield[2],0.009*dataYield[3]};

Double_t dataYield_systErrREL[4];

//-----------------------------ZZ,Zg-------------------------------------------------------------------
//MC yield for ZZ (7.5% syst)
Double_t NbkgEstZZ[4]={2.01,3.47,2.66,5.09};
Double_t ZZ_statErr[4] = {0.04, 0.04, 0.03, 0.06};
//Double_t NbkgEstZZErr[4]={0.075*NbkgEstZZ[0],0.075*NbkgEstZZ[1],0.075*NbkgEstZZ[2],0.075*NbkgEstZZ[3]};
//Double_t NbkgEstZZErr[4]={0.075*NbkgEstZZ[0],0.075*NbkgEstZZ[1],0.075*NbkgEstZZ[2],0.075*NbkgEstZZ[3]};
Double_t NbkgEstZZErr[4]={0.14*NbkgEstZZ[0],0.14*NbkgEstZZ[1],0.14*NbkgEstZZ[2],0.14*NbkgEstZZ[3]};

Double_t zzPart[4];
Double_t zzPartSF[4];

//ZZ acceptance systematics
Double_t ZZ_muscaleErrUp[4]={0.,0.00536498,0.00410551,0.00902297};
Double_t ZZ_muscaleErrDown[4]={0.,0.00547153,0.00325055,0.00833};

Double_t ZZ_escaleErrUp[4]={0.0120834,0.00934668,0.0081268,0.};
Double_t ZZ_escaleErrDown[4]={0.0187905,0.0101121,0.0113352,0.};

Double_t ZZ_puErrUp[4]={0.04155,0.01425,0.0456,0.0181};
Double_t ZZ_puErrDown[4]={0.04187,0.01537,0.0460,0.0194};

Double_t ZZ_sfErr[4]={0.0268,0.0243,0.0186,0.0129};

Double_t ZZ_METErr[4]={0.13986,0.0774586,0.18889,0.070955};

//Zgamma MC: 13%
Double_t NbkgEstZG[4]={0.,0.,0.538175,0.};
Double_t ZG_statErr[4] = {0., 0., 0.538175, 0.};
//Double_t NbkgEstZGErr[4]={0.13*NbkgEstZG[0],0.13*NbkgEstZG[1],0.13*NbkgEstZG[2],0.13*NbkgEstZG[3]};
Double_t NbkgEstZGErr[4]={0.07*NbkgEstZG[0],0.07*NbkgEstZG[1],0.07*NbkgEstZG[2],0.07*NbkgEstZG[3]};

Double_t zgPartSF[4];

Double_t ZG_muscaleErr[4]={0.,0,0.00005,0};
Double_t ZG_escaleErr[4]={0.,0,0,0};

Double_t ZG_puErr[4]={0,0,0.35,0};

Double_t ZG_sfErr[4]={0,0,0.0169,0};

Double_t ZG_METErr[4]={0.00,0.00,0.0053,0.0};


//ZZ and ZG xsec.: correlated between channel combinations
Double_t NbkgEstZZErrREL[4];
Double_t NbkgEstZGErrREL[4];

Double_t NbkgEstWWZErrREL[4];
Double_t NbkgEstWWWErrREL[4];
Double_t NbkgEstZZZErrREL[4];
Double_t NbkgEstWZZErrREL[4];
Double_t NbkgEstttWErrREL[4];
Double_t NbkgEstttZErrREL[4];

Double_t NbkgEstVVVErrREL[4];
Double_t NbkgEstttVErrREL[4];

Double_t YieldVVVErrREL[4];
Double_t YieldttVErrREL[4];
//new backgrounds!

/*
Double_t NbkgEst_WWZ[4]={0,0,0,0};
Double_t NbkgEst_WWW[4]={0,0,0,0};
Double_t NbkgEst_ZZZ[4]={0,0,0,0};
Double_t NbkgEst_WZZ[4]={0,0,0,0};
Double_t NbkgEst_ttW[4]={0,0,0,0};
Double_t NbkgEst_ttZ[4]={0,0,0,0};
*/

Double_t NbkgEst_WWZ[4]={0.295,0.392,0.503,0.623};
Double_t NbkgEst_WWW[4]={0.044,0.074,0.076,0.088};
Double_t NbkgEst_ZZZ[4]={0.010,0.009,0.019,0.016};
Double_t NbkgEst_WZZ[4]={0.106,0.134,0.162,0.224};
Double_t NbkgEst_ttW[4]={0.092,0.176,0.228,0.256};
Double_t NbkgEst_ttZ[4]={1.055,1.253,1.373,1.756};
//Double_t NbkgEst_ttZ[4]={0.836,1.040,1.404,1.876};




Double_t NbkgEst_WWZErr[4]={NbkgEst_WWZ[0]*0.5,NbkgEst_WWZ[1]*0.5,NbkgEst_WWZ[2]*0.5,NbkgEst_WWZ[3]*0.5};
Double_t NbkgEst_WWWErr[4]={NbkgEst_WWW[0]*0.5,NbkgEst_WWW[1]*0.5,NbkgEst_WWW[2]*0.5,NbkgEst_WWW[3]*0.5};
Double_t NbkgEst_ZZZErr[4]={NbkgEst_ZZZ[0]*0.5,NbkgEst_ZZZ[1]*0.5,NbkgEst_ZZZ[2]*0.5,NbkgEst_ZZZ[3]*0.5};
Double_t NbkgEst_WZZErr[4]={NbkgEst_WZZ[0]*0.5,NbkgEst_WZZ[1]*0.5,NbkgEst_WZZ[2]*0.5,NbkgEst_WZZ[3]*0.5};
Double_t NbkgEst_ttWErr[4]={NbkgEst_ttW[0]*0.5,NbkgEst_ttW[1]*0.5,NbkgEst_ttW[2]*0.5,NbkgEst_ttW[3]*0.5};
Double_t NbkgEst_ttZErr[4]={NbkgEst_ttZ[0]*0.5,NbkgEst_ttZ[1]*0.5,NbkgEst_ttZ[2]*0.5,NbkgEst_ttZ[3]*0.5};

/*
Double_t NbkgEst_WWZErr[4]={0,0,0,0};
Double_t NbkgEst_WWWErr[4]={0,0,0,0};
Double_t NbkgEst_ZZZErr[4]={0,0,0,0};
Double_t NbkgEst_WZZErr[4]={0,0,0,0};
Double_t NbkgEst_ttWErr[4]={0,0,0,0};
Double_t NbkgEst_ttZErr[4]={0,0,0,0};
*/

Double_t NbkgEst_VVVErr[4]= {0,0,0,0};
Double_t NbkgEst_ttVErr[4]= {0,0,0,0};


Double_t Yield_VVV[4]={0,0,0,0};
Double_t Yield_ttV[4]={0,0,0,0};
Double_t Yield_VVVErr[4]={0,0,0,0};
Double_t Yield_ttVErr[4]={0,0,0,0};

for (size_t nw=0;nw<4;nw++) {
  NbkgEst_VVVErr[nw]=NbkgEst_WWZErr[nw]+NbkgEst_WWWErr[nw]+NbkgEst_ZZZErr[nw]+NbkgEst_WZZErr[nw];
  NbkgEst_ttVErr[nw]=NbkgEst_ttWErr[nw]+NbkgEst_ttZErr[nw];

  Yield_VVV[nw]=NbkgEst_WWZ[nw]+NbkgEst_WWW[nw]+NbkgEst_ZZZ[nw]+NbkgEst_WZZ[nw];
  Yield_ttV[nw]=NbkgEst_ttW[nw]+NbkgEst_ttZ[nw];

  Yield_VVVErr[nw]=0.25*Yield_VVV[nw];
  Yield_ttVErr[nw]=0.2*Yield_ttV[nw];
}
//------------------------------Tau fraction ---------------------------------------------------------
//WZ madgraph Fall11 numbers (signal and tau)
Double_t NCountWZMC[4]={12563.8,13987.9,15717.2,20741.9};
Double_t NCountWZMC_forErr[4]={12448,14706,17940,22705};

Double_t WZTauFrac[4]={813.58,1023.5,1004.52,1526.76};
Double_t WZTauHadronic[4]={0.02,0.004,0.02,0.002};
Double_t WZTauFrac_forErr[4]={796,1059,1144,1603};
Double_t WZTauHadronic_forErr[4]={4,1,6,2};

Double_t WZTauFrac_Err2[4]={0,0,0,0};

for (size_t i=0;i<4;i++) {
  WZTauFrac[i]=(WZTauFrac[i]-WZTauHadronic[i])/NCountWZMC[i];
  //WZTauFrac_forErr[i]=(WZTauFrac_forErr[i]-WZTauHadronic_forErr[i])/NCountWZMC_forErr[i];
  WZTauFrac_forErr[i]=(WZTauFrac_forErr[i]/*-WZTauHadronic_forErr[i]*/)/NCountWZMC_forErr[i];
  WZTauFrac_Err2[i]=(WZTauFrac_forErr[i]*(1.-WZTauFrac_forErr[i]))/NCountWZMC_forErr[i];//binomial error
  cout << " tau fraction " << i << " " << WZTauFrac[i] << " Err:" << sqrt(WZTauFrac_Err2[i]) << endl;
}

Double_t oneMinusTauFrac[4];
Double_t oneMinusTauFrac_Err2[4];

for (i=0;i<4;i++) {
 oneMinusTauFrac[i] = 1. - WZTauFrac[i];
 oneMinusTauFrac_Err2[i] = WZTauFrac_Err2[i];
}
//------------------------------Tau fraction end------------------------------------------------------

Double_t PURW_Err2[4];
Double_t rhoPDF_Err2[4];
Double_t rhoMET_Err2[4];
Double_t rhoMu_ScaleErr2[4];
Double_t rhoEle_ScaleErr2[4];
Double_t rho_totalErr[4];


////////////////////////////
/* cross section */

//numerator of cross section
Double_t xsection_Num[4];

Double_t xsection_Num_Err[4];
Double_t xsection_Num_Err2[4];

Double_t xsection_Num_ErrSyst[4];
Double_t xsection_Num_ErrSyst2[4];

//denominator of cross section
Double_t xsection_Den[4];
Double_t xsection_Den_Err[4];
Double_t xsection_Den_Err2[4];

Double_t xsection_Den_ErrSyst[4];
Double_t xsection_Den_ErrSyst2[4];

//xsection
Double_t xsection[4];
Double_t xsection_ErrStat[4];
Double_t xsection_ErrSyst[4];
Double_t xsection_ErrLumi[4];

Double_t xsection_method2[4];

//BRs
double BR_ratios[4]={0.03363*0.1075,0.03363*0.1057,0.03366*0.1075,0.03366*0.1057};
double BR_method2 = 3.*0.033658*(0.1125+0.1075+0.1057);

double bkgHelper[4]={0,0,0,0};

//xsection
for (size_t i=0;i<4;i++) {

  //numerator
  //
  bkgHelper[i]=NbkgEst_WWZ[i] + NbkgEst_WWW[i]  + NbkgEst_ZZZ[i] + NbkgEst_WZZ[i] + NbkgEst_ttW[i] + NbkgEst_ttZ[i];
  std::cout <<"added backgrounds " << bkgHelper[i] << std::endl;
  //
  xsection_Num[i] = (dataYield[i] - NbkgEstZZ[i] - NbkgEstZG[i] - bkgHelper[i])*oneMinusTauFrac[i];
  //stat.error
  xsection_Num_Err2[i]=pow(dataYield_Err[i]*oneMinusTauFrac[i],2); 
  xsection_Num_Err[i]=sqrt(xsection_Num_Err2[i]);

  //syst.err
  //xsection_Num_ErrSyst2[i]=pow(oneMinusTauFrac[i],2) * (
  //                              pow(NbkgEstZZErr[i],2) // MC
  //                              + pow(NbkgEstZGErr[i],2) // MC
  //                              + pow(dataYield_systErr[i],2) //MM syst
  //                         )
  //                         +pow(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i],2)*oneMinusTauFrac_Err2[i];
  //
  //xsection_Num_ErrSyst[i]=sqrt(xsection_Num_ErrSyst2[i]);


  //cout << i << " systeff. ZZ: " << 100*NbkgEstZZErr[i]/xsection_Num[i] <<"% \n";
  //cout << i << " systeff. ZG: " << 100*NbkgEstZGErr[i]/xsection_Num[i] <<"% \n";


  cout<<  i << "  MM: "         << 100*dataYield_systErr[i]/xsection_Num[i] << "%"<< endl;


  //denominator
  xsection_Den[i]        = acceptance[i]*intLumi;
  xsection_Den_Err[i]    = 0;//was stat.err
  xsection_Den_Err2[i]   = 0;

  xsection[i]=xsection_Num[i]/xsection_Den[i];


  xsection_ErrStat[i] = sqrt ( pow( xsection_Num_Err[i]*xsection[i]/xsection_Num[i],2) 
                             + pow( xsection_Den_Err[i]*xsection[i]/xsection_Den[i],2));


  xsection_ErrLumi[i]=xsection[i]*intLumiErr/intLumi;
  //systematics-------------------------------------------------------------------------

  PURW_ErrUp[i] = fabs(NbkgEstZZ[i]*(PURW_ErrUp0[i]-ZZ_puErrUp[i])+ NbkgEstZG[i]*(PURW_ErrUp0[i]-ZG_puErr[i])
                       +(bkgHelper[i]-dataYield[i])*PURW_ErrUp0[i])/((1+PURW_ErrUp0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]- bkgHelper[i]));
  PURW_ErrDown[i] = fabs(NbkgEstZZ[i]*(PURW_ErrDown0[i]-ZZ_puErrDown[i])+ NbkgEstZG[i]*(PURW_ErrDown0[i]-ZG_puErr[i])
                       +(bkgHelper[i]-dataYield[i])*PURW_ErrDown0[i])/((1-PURW_ErrDown0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i] - bkgHelper[i]));
  PURW_Err[i] = std::max(PURW_ErrUp[i],PURW_ErrDown[i]);

  //std::cout << " pu " << PURW_Err[i] << endl; 
 
  rhoMu_ScaleErrUp[i] = fabs(NbkgEstZZ[i]*(rhoMu_ScaleErrUp0[i]-ZZ_muscaleErrUp[i])+ NbkgEstZG[i]*(rhoMu_ScaleErrUp0[i]-ZG_muscaleErr[i])
                       +(bkgHelper[i]-dataYield[i])*rhoMu_ScaleErrUp0[i])/((1+rhoMu_ScaleErrUp0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]- bkgHelper[i]));
  rhoMu_ScaleErrDown[i] = fabs(NbkgEstZZ[i]*(rhoMu_ScaleErrDown0[i]-ZZ_muscaleErrDown[i])+ NbkgEstZG[i]*(rhoMu_ScaleErrDown0[i]-ZG_muscaleErr[i])
                       +(bkgHelper[i]-dataYield[i])*rhoMu_ScaleErrDown0[i])/((1+rhoMu_ScaleErrDown0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]- bkgHelper[i]));
  rhoMu_ScaleErr[i] = std::max(rhoMu_ScaleErrUp[i],rhoMu_ScaleErrDown[i]);

  //std::cout << " mu " << rhoMu_ScaleErr[i] << endl; 

  rhoEle_ScaleErrUp[i] = fabs(NbkgEstZZ[i]*(rhoEle_ScaleErrUp0[i]-ZZ_escaleErrUp[i])+ NbkgEstZG[i]*(rhoEle_ScaleErrUp0[i]-ZG_escaleErr[i])
                       +(bkgHelper[i]-dataYield[i])*rhoEle_ScaleErrUp0[i])/((1+rhoEle_ScaleErrUp0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]- bkgHelper[i]));
  rhoEle_ScaleErrDown[i] = fabs(NbkgEstZZ[i]*(rhoEle_ScaleErrDown0[i]-ZZ_escaleErrDown[i])+ NbkgEstZG[i]*(rhoEle_ScaleErrDown0[i]-ZG_escaleErr[i])
                       +(bkgHelper[i]-dataYield[i])*rhoEle_ScaleErrDown0[i])/((1+rhoEle_ScaleErrDown0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]- bkgHelper[i]));
  rhoEle_ScaleErr[i] = std::max(rhoEle_ScaleErrUp[i],rhoEle_ScaleErrDown[i]);

  //std::cout << " ele " << rhoEle_ScaleErr[i] << endl; 

  //down variation larger for symmetrical errors (ZZ)
  //double combpdf = sqrt(pow(rhoPDF_Err0[i],2) + pow(rhoQFScale_Err0[i],2));
  //rhoPDF_Err[i] = fabs(dataYield[i]*combpdf)/((1-combpdf)*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]));
  rhoPDF_Err[i] = sqrt(pow(rhoPDF_Err0[i],2) + pow(rhoQFScale_Err0[i],2));//fabs(dataYield[i]*rhoPDF_Err0[i])/((1-rhoPDF_Err0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]));//hack
  rhoMET_Err[i] = fabs(NbkgEstZZ[i]*(rhoMET_Err0[i]-ZZ_METErr[i])+ NbkgEstZG[i]*(rhoMET_Err0[i]-ZG_METErr[i])
                       +(bkgHelper[i]-dataYield[i])*rhoMET_Err0[i])/((1-rhoMET_Err0[i])*(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]- bkgHelper[i]));

  //std::cout <<"pu:"<< PURW_Err[i] <<" musc:" << rhoMu_ScaleErr[i] <<" esc:"<< rhoEle_ScaleErr[i]
  //          <<" pdf:"<<  rhoPDF_Err[i] << " met:"<< rhoMET_Err[i] << endl;
  //std::cout << " pdf " << rhoPDF_Err[i] << endl; 
  //std::cout << " met " << rhoMET_Err[i] << endl; 
  PURW_Err2[i]=pow(PURW_Err[i],2);
  rhoPDF_Err2[i]=pow(rhoPDF_Err[i],2);
  rhoMET_Err2[i]=pow(rhoMET_Err[i],2);
  rhoMu_ScaleErr2[i]=pow(rhoMu_ScaleErr[i],2);
  rhoEle_ScaleErr2[i]=pow(rhoEle_ScaleErr[i],2);

  Double_t rho_totalErr2 = rhoPDF_Err2[i]+rhoMET_Err2[i]+rhoMu_ScaleErr2[i]+rhoEle_ScaleErr2[i]+PURW_Err2[i];
  rho_totalErr[i]=sqrt(rho_totalErr2);

  Double_t xsec_ZZ_sfErr=oneMinusTauFrac[i]*NbkgEstZZ[i]*ZZ_sfErr[i]/xsection_Den[i];
  Double_t xsec_ZG_sfErr=oneMinusTauFrac[i]*NbkgEstZG[i]*ZG_sfErr[i]/xsection_Den[i];
  zzPartSF[i]=(xsec_ZZ_sfErr)/(xsection[i]*acceptance_systErr0[i]/acceptance[i]);
  zgPartSF[i]=(xsec_ZG_sfErr)/(xsection[i]*acceptance_systErr0[i]/acceptance[i]);

  //ZZ,Zgamma stat err + tau(not included?) //for table...
  Double_t xsec_ZZ_statErr=oneMinusTauFrac[i]*ZZ_statErr[i]/xsection_Den[i];
  Double_t xsec_ZG_statErr=oneMinusTauFrac[i]*ZG_statErr[i]/xsection_Den[i];

  //Double_t xsec_VVV_accErr=oneMinusTauFrac[i]*Yield_VVVErr/xsection_Den[i];
  //Double_t xsec_ttV_accErr=oneMinusTauFrac[i]*Yield_ttVErr/xsection_Den[i];

  //Double_t otau = (sqrt(oneMinusTauFrac_Err2[i])/oneMinusTauFrac[i])*xsection[i];
  //std::cout << " tau+ZZ+Zg stat.: " << sqrt(xsec_ZZ_statErr*xsec_ZZ_statErr + xsec_ZG_statErr*xsec_ZG_statErr + otau*otau)/xsection[i] << endl;

  double accsystErrRelXS = (dataYield[i] - NbkgEstZZ[i]*(1-ZZ_sfErr[i]) - NbkgEstZG[i]*(1-ZG_sfErr[i])-bkgHelper[i])*oneMinusTauFrac[i];
  accsystErrRelXS/=intLumi*acceptance[i]*(1.-acceptance_systErr0[i]/acceptance[i]);
  acceptance_systErr[i]=acceptance[i]*fabs(accsystErrRelXS-xsection[i])/xsection[i];
  //std::cout << " SF " << fabs(accsystErrRelXS-xsection[i])/xsection[i]  << " input:" << acceptance_systErr0[i]/acceptance[i] << endl; 

  acceptance_Err[i]= sqrt(pow(acceptance_statErr[i],2)+pow(acceptance_systErr[i],2));

  xsection_Num_ErrSyst2[i]=pow(oneMinusTauFrac[i],2) * (
                                pow(NbkgEstZZErr[i],2) // MC xsection
                                + pow(NbkgEstZGErr[i],2) // MC
                                +pow(NbkgEst_WWZErr[i],2) //new bkgs 
                                +pow(NbkgEst_WWWErr[i],2)
                                +pow(NbkgEst_ZZZErr[i],2)
                                +pow(NbkgEst_WZZErr[i],2)
                                +pow(NbkgEst_ttWErr[i],2)
                                +pow(NbkgEst_ttZErr[i],2)

                                +pow(ZZ_statErr[i],2) // MC stat
                                +pow(ZG_statErr[i],2) // MC stat
 
                                +pow(Yield_VVVErr[i],2)
                                +pow(Yield_ttVErr[i],2)

                                + pow(dataYield_systErr[i],2) //MM syst
                           )
                           +pow(dataYield[i]-NbkgEstZZ[i]-NbkgEstZG[i]-bkgHelper[i],2)*oneMinusTauFrac_Err2[i]
                           //+xsec_ZZ_statErr*xsec_ZZ_statErr + xsec_ZG_statErr*xsec_ZG_statErr
                           //+xsec_VVV_accErr*xsec_VVV_accErr+xsec_ttV_accErr*xsec_ttV_accErr;
                           ;

  xsection_Num_ErrSyst[i]=sqrt(xsection_Num_ErrSyst2[i]);


  xsection_Den_ErrSyst2[i]= pow(rho_totalErr[i]*xsection_Den[i],2) + pow(acceptance_Err[i]*xsection_Den[i]/acceptance[i],2);
  xsection_Den_ErrSyst[i] = sqrt( xsection_Den_ErrSyst2[i]);

  xsection_ErrSyst[i] = sqrt ( pow( xsection_Num_ErrSyst[i]*xsection[i]/xsection_Num[i],2)
                             + pow( xsection_Den_ErrSyst[i]*xsection[i]/xsection_Den[i],2));


  cout << "|ch-only       | " << xsection[i] << " +- " << xsection_ErrStat[i]<< "(stat) +- " 
       << xsection_ErrSyst[i] << "(syst) +- "<< xsection_ErrLumi[i] << "(lumi)" << endl;

  cout << "|inclusive st.1| " << xsection[i]/BR_ratios[i] << " +- "<<xsection_ErrStat[i]/BR_ratios[i]<<" +- " << xsection_ErrSyst[i]/BR_ratios[i]
       << " +- "<< xsection_ErrLumi[i]/BR_ratios[i] << std::endl;

  cout << "|inclusive st.2| " << (dataYield[i] - NbkgEstZZ[i] - NbkgEstZG[i]-bkgHelper[i])/(intLumi*acceptance_method2[i]*BR_method2) << std::endl;
  cout << "ratio: " << xsection_ErrStat[i]/xsection[i] << "(stat.) " << xsection_ErrSyst[i]/xsection[i] << "(syst.)" << endl;
  //inclusive
  xsection_method2[i]=(dataYield[i] - NbkgEstZZ[i] - NbkgEstZG[i]-bkgHelper[i])/(intLumi*acceptance_method2[i]*BR_method2);

  //calculate propagated parts for syst error of ZZ and Zg:
  NbkgEstZZErrREL[i]=(oneMinusTauFrac[i]*NbkgEstZZErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstZGErrREL[i]=(oneMinusTauFrac[i]*NbkgEstZGErr[i])/(xsection_Den[i]*xsection[i]);
/*
  NbkgEstWWZErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_WWZErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstWWWErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_WWWErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstZZZErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_ZZZErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstWZZErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_WZZErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstttWErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_ttWErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstttZErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_ttZErr[i])/(xsection_Den[i]*xsection[i]);
*/
  NbkgEstVVVErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_VVVErr[i])/(xsection_Den[i]*xsection[i]);
  NbkgEstttVErrREL[i]=(oneMinusTauFrac[i]*NbkgEst_ttVErr[i])/(xsection_Den[i]*xsection[i]);


  YieldVVVErrREL[i]=(oneMinusTauFrac[i]*Yield_VVVErr[i])/(xsection_Den[i]*xsection[i]);
  YieldttVErrREL[i]=(oneMinusTauFrac[i]*Yield_ttVErr[i])/(xsection_Den[i]*xsection[i]);


//std::cout << "zz:"<<NbkgEstZZErrREL[i]<< " zg:"<<NbkgEstZGErrREL[i]<<endl;
  dataYield_systErrREL[i]=(oneMinusTauFrac[i]*dataYield_systErr[i])/(xsection_Den[i]*xsection[i]);
  //std::cout << dataYield_systErrREL[i] << endl;
  xsection_ErrSyst[i]/=BR_ratios[i];
  xsection_ErrStat[i]/=BR_ratios[i];
  xsection[i]/=BR_ratios[i];
  cout << std::endl;
}


/*
  BLUE Method
*/

//scale factor components
double sb_ereco[4]={0.995,0.914,0.752,0};
double sb_mubase[4]={0,0.265,0.425,0.680};
double sb_mutrg[4]={0,0,0.0704,0.00647};
double sb_etrg[4]={0.000088,0.0108,0,0};
double accyield[4]={41.82,46.15,52.36,68.41};

double zzsb_ereco[4]={0.0478,0.06945,0.0382,0};
double zzsb_mubase[4]={0,0.02,0.029,0.059};
double zzsb_mutrg[4]={0,0,0.00367,0.000511};
double zzsb_etrg[4]={0.0000042,0.00081,0,0};

//lumi systematics
 double lumiSyst[4]={0.022, 0.022, 0.022, 0.022};
if (0) //tiny effect
for (size_t k=0;k<4;k++) {
  sb_ereco[k]+=zzPart[k]*zzsb_ereco[k]/NbkgEstZZ[k];
  sb_mubase[k]+=zzPart[k]*zzsb_mubase[k]/NbkgEstZZ[k];
  sb_mutrg[k]+=zzPart[k]*zzsb_mutrg[k]/NbkgEstZZ[k];
  sb_etrg[k]+=zzPart[k]*zzsb_etrg[k]/NbkgEstZZ[k];
}
for (size_t k=0;k<4;k++) {
  sb_ereco[k]*=1./accyield[k];
  sb_mubase[k]*=1./accyield[k];
  sb_mutrg[k]*=1./accyield[k];
  sb_etrg[k]*=1./accyield[k];

}

//compose Error Matrix [i = column*4+row] [0-3][0-3] indices
//matrix is symmetric
//channel legend: 0 = 3e, 1 = 2e1mu, 2 = 2mu1e,  3 = 3mu

Double_t elm[16];
for (size_t i=0;i<16;i++) elm[i]=0;
//diagonals:
elm[0]=pow(xsection_ErrSyst[0],2) + pow(xsection_ErrStat[0],2);
elm[5]=pow(xsection_ErrSyst[1],2) + pow(xsection_ErrStat[1],2);
elm[10]=pow(xsection_ErrSyst[2],2) + pow(xsection_ErrStat[2],2);
elm[15]=pow(xsection_ErrSyst[3],2) + pow(xsection_ErrStat[3],2);

//non-diagonal elements containing correlated systematics (scaled to cross-section)

double commonSystMat[4][4];
double PUCombMat[4][4];
double EleScaleCombMat[4][4];
double MuScaleCombMat[4][4];
 double lumi_matrix[4][4];
for (int cha=0;cha<4;cha++) {
  for (int chb=0;chb<4;chb++) {
    PUCombMat[cha][chb]=max(fabs(PURW_ErrUp[cha]*PURW_ErrUp[chb]),fabs(PURW_ErrDown[cha]*PURW_ErrDown[chb]));
    EleScaleCombMat[cha][chb]= max(fabs(rhoEle_ScaleErrUp[cha]*rhoEle_ScaleErrUp[chb]),
                                        fabs(rhoEle_ScaleErrDown[cha]*rhoEle_ScaleErrDown[chb]));
    MuScaleCombMat[cha][chb]= max(fabs(rhoMu_ScaleErrUp[cha]*rhoMu_ScaleErrUp[chb]),
                                        fabs(rhoMu_ScaleErrDown[cha]*rhoMu_ScaleErrDown[chb]));
    lumi_matrix[cha][chb]=lumiSyst[cha]*lumiSyst[chb]*xsection[cha]*xsection[chb];
   commonSystMat[cha][chb]=rhoMET_Err[cha]*rhoMET_Err[chb]
                          +rhoPDF_Err[chb]*rhoPDF_Err[cha]
                          +NbkgEstZZErrREL[cha]*NbkgEstZZErrREL[chb]
                          +NbkgEstVVVErrREL[cha]*NbkgEstVVVErrREL[chb]
                          +NbkgEstttVErrREL[cha]*NbkgEstttVErrREL[chb]
                          +YieldVVVErrREL[cha]*YieldVVVErrREL[chb]
                          +YieldttVErrREL[cha]*YieldttVErrREL[chb]
                          +PUCombMat[cha][chb]
                          +lumiSyst[cha]*lumiSyst[chb];  
                          /*+NbkgEstWWZErrREL[cha]*NbkgEstWWZErrREL[chb]
                          +NbkgEstWWWErrREL[cha]*NbkgEstWWWErrREL[chb]
                          +NbkgEstZZZErrREL[cha]*NbkgEstZZZErrREL[chb]
                          +NbkgEstWZZErrREL[cha]*NbkgEstWZZErrREL[chb]
                          +NbkgEstttWErrREL[cha]*NbkgEstttWErrREL[chb]
                          +NbkgEstttZErrREL[cha]*NbkgEstttZErrREL[chb]*/
  }
}

//channels 0 and 1
elm[4]=elm[1]=xsection[0]*xsection[1] * ( sb_ereco[0]*sb_ereco[1]
                                          +sb_etrg[0]*sb_etrg[1]
                                          +EleScaleCombMat[0][1]
                                          +commonSystMat[0][1]);

//channels 0 and 2
elm[8]=elm[2]=xsection[0]*xsection[2]* (  sb_ereco[0]*sb_ereco[2]
                                          +EleScaleCombMat[0][2]
                                          +commonSystMat[0][2]);

//channels 0 and 3
elm[12]=elm[3]=xsection[0]*xsection[3]*(
                                          commonSystMat[0][3]);

//element channels 1 and 2
elm[9]=elm[6]=xsection[1]*xsection[2]* (  sb_ereco[1]*sb_ereco[2]
                                          +sb_mubase[1]*sb_mubase[2]
                                          +EleScaleCombMat[1][2]
                                          +MuScaleCombMat[1][2]
                                          +commonSystMat[1][2]);

//channels 1 and 3
elm[13]=elm[7]=xsection[1]*xsection[3]*(  sb_mubase[1]*sb_mubase[3]
                                          +MuScaleCombMat[1][3]
                                          +commonSystMat[1][3]);

//channels 2 and 3
elm[14]=elm[11]=xsection[2]*xsection[3]*( sb_mubase[2]*sb_mubase[3]
                                          +sb_mutrg[2]*sb_mutrg[3]
                                          +MuScaleCombMat[2][3]
                                          +commonSystMat[2][3]);

//for (size_t i=0;i<16;i++) if (i!=0 && i!=5 && i!=10 && i!=15) elm[i]=0;
cout << endl;

for (size_t i=0;i<4;i++) {
  for (size_t j=0;j<4;j++) {
  cout << "       " << elm[j*4+i] << "       ";
  }
  cout << endl;
}
cout << endl;

//initialize matrix for inversion
TMatrixD errMat(4,4,elm);
TMatrixD errMatInv(4,4);

TMatrixD errMatCopy(errMat);

//invert !
errMatInv = errMat.Invert();

//print inverse
Double_t *mRefTest = errMat.GetMatrixArray();

cout << endl << " Matrix Inverse:" << endl;
for (size_t i=0;i<4;i++) {
  for (size_t j=0;j<4;j++) {
  cout << "       " << mRefTest[j*4+i] << "       ";
  }
  cout << endl;
}
cout << endl;


//get norm, and alpha factors for each channel
Double_t *mRef= errMat.GetMatrixArray();
Double_t norm=0.;
Double_t alphaCH[4]={0.,0.,0.,0.};

for (size_t i=0;i<16;i++) norm+=mRef[i];

for (size_t i=0;i<4;i++) {
  for (size_t j=0;j<4;j++) {
    alphaCH[i]+=mRef[i*4+j];
  }
  alphaCH[i]/=norm;
}

cout << "al0 " << alphaCH[0] << " al1 " << alphaCH[1] << " al2 " << alphaCH[2] << " al3 " << alphaCH[3] << endl;
//cout << "consistency check:" << alphaCH[0]+alphaCH[1]+alphaCH[2]+alphaCH[3] <<endl;
cout << endl;

Double_t final_Xsec = alphaCH[0]*xsection[0] + alphaCH[1]*xsection[1] + alphaCH[2]*xsection[2] + alphaCH[3]*xsection[3];

Double_t final_Xsec_method2 =   alphaCH[0]*xsection_method2[0] + alphaCH[1]*xsection_method2[1]
                              +alphaCH[2]*xsection_method2[2] + alphaCH[3]*xsection_method2[3];

Double_t combined_error=0;
Double_t *copyRef = errMatCopy.GetMatrixArray();
for (size_t i=0;i<4;i++)
  for (size_t j=0;j<4;j++) {
    combined_error+=alphaCH[i]*alphaCH[j]*copyRef[i*4+j];
  }

double lumi_error(0);
 for (int ier=0; ier<4; ier++){
   for (int jer=0; jer<4; jer++){
     lumi_error+=alphaCH[ier]*alphaCH[jer]*lumi_matrix[ier][jer];
   }
 }

Double_t stat_err_tot2 =  pow(alphaCH[0]*xsection_ErrStat[0],2) + pow(alphaCH[1]*xsection_ErrStat[1],2)
                         +pow(alphaCH[2]*xsection_ErrStat[2],2) + pow(alphaCH[3]*xsection_ErrStat[3],2);


Double_t stat_err_tot = sqrt(stat_err_tot2);
Double_t syst_err_tot = sqrt (combined_error - stat_err_tot2 - lumi_error);

combined_error=sqrt(combined_error);


//cout << "combined sigma(WZ) = " << final_Xsec << " +- " << stat_err_tot << "(stat.) +- " << syst_err_tot << "(syst) +- "
//   << 0.022*final_Xsec << " (lumi) pb " << endl;
cout << "combined sigma(WZ) = " << final_Xsec << " +- " << stat_err_tot << "(stat.) +- " << syst_err_tot << "(syst) +- "
     << sqrt(lumi_error) << " (lumi) pb " << endl;

cout << endl;
cout << "combined method2 = " << final_Xsec_method2 << endl; 


//old cout << "THEORETICAL NLO sigma~=" << (1/3.)*(1/3.)*0.594 << " pb " << endl;

cout << endl << endl;

return;
/*
Double_t *mRefTest2 = errMatCopy.GetMatrixArray();

cout << endl << " Matrix Print:" << endl;
for (size_t i=0;i<4;i++) {
  for (size_t j=0;j<4;j++) {
  cout << "       " << mRefTest2[j*4+i]*10000 << "       ";
  }
  cout << endl;
}
cout << " * 10^-4";
cout << endl;

return;
*/


}

