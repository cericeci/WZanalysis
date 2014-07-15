#include "numbers.h"
#include "numMC.h"
#include "numMM_met.h"
//#include "numMM.h"
//#include "numGEN.h"
//#include "numData.h"
#include "numGEN_met.h"
#include "numData_met.h"
#include "syst.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

void crossSection()
{
  //  bool latexOutput(false);
  bool latexOutput(false);
  bool outputNumbers(false);
  //0: eee, 1:eem, 2:emm, 3:mmm
  ofstream fileNum;

  if (outputNumbers){
    fileNum.open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/numFinalWoPu.txt");
  }
  double Axe3eS(dAxe3eS), Axe2e1muS(dAxe2e1muS), Axe1e2muS(dAxe1e2muS), Axe3muS(dAxe3muS);    //Spanish acceptance times eff
  double AxeSJ[4]={dAxe3eSJ, dAxe2e1muSJ, dAxe1e2muSJ, dAxe3muSJ};
  //double sAxeSJ[4]={dsAxe3eSJ, dsAxe2e1muSJ, dsAxe1e2muSJ, dsAxe3muSJ};
  double AxeS[4]={dAxe3eS, dAxe2e1muS, dAxe1e2muS, dAxe3muS};
  double sAxeS[4]={dsAxe3eS, dsAxe2e1muS, dsAxe1e2muS, dsAxe3muS};
  //  double AxeZ[4]={dAxe3eZ, dAxe2e1muZ, dAxe1e2muZ, dAxe3muZ};
  double AxeZ[4]={dAxe3eZ_2, dAxe2e1muZ_2, dAxe1e2muZ_2, dAxe3muZ_2};
  double sAxeZ[4]={dsAxe3eZ, dsAxe2e1muZ, dsAxe1e2muZ, dsAxe3muZ};
  double N_good[4]={dN_good3e, dN_good2e1mu, dN_good1e2mu, dN_good3mu};
  double sN_good[4]={dsN_good3e, dsN_good2e1mu, dsN_good1e2mu, dsN_good3mu};
  double N_fake[4]={dN_fake3e, dN_fake2e1mu, dN_fake1e2mu, dN_fake3mu};
  double sN_fake[4]={dsN_fake3e, dsN_fake2e1mu, dsN_fake1e2mu, dsN_fake3mu};
  double N_data[4]={dN_data3e, dN_data2e1mu, dN_data1e2mu, dN_data3mu};
  double NZZ[4]={dNZZ_3e, dNZZ_2e1mu, dNZZ_1e2mu, dNZZ_3mu};
  double sNZZ[4]={dsNZZ_3e, dsNZZ_2e1mu, dsNZZ_1e2mu, dsNZZ_3mu};
  double NZgamma[4]={dNZgamma_3e, dNZgamma_2e1mu, dNZgamma_1e2mu, dNZgamma_3mu};
  double sNZgamma[4]={dsNZgamma_3e, dsNZgamma_2e1mu, dsNZgamma_1e2mu, dsNZgamma_3mu};
  double NWV[4]={dNWV_3e, dNWV_2e1mu, dNWV_1e2mu, dNWV_3mu};
  double sNWV[4]={dsNWV_3e, dsNWV_2e1mu, dsNWV_1e2mu, dsNWV_3mu};
  double NVVV[4]={dNVVV_3e, dNVVV_2e1mu, dNVVV_1e2mu, dNVVV_3mu};
  double sNVVV[4]={dsNVVV_3e, dsNVVV_2e1mu, dsNVVV_1e2mu, dsNVVV_3mu};
  double tauFactor[4]={dtauFactor3e, dtauFactor2e1mu, dtauFactor1e2mu, dtauFactor3mu};
  double stauFactor[4]={dstauFactor3e, dstauFactor2e1mu, dstauFactor1e2mu, dstauFactor3mu};
  double crossSectionSpanish[4], crossSectionZagreb[4],crossSectionSpanishJonatan[4], errorCsSpanish[4], errorCsZagreb[4];
  double errorCsSpNum[4], errorCsSpDen[4], errorCsZgNum[4], errorCsZgDen[4], errorCsSpanish[4], errorCsZagreb[4], csSpNum[4], csSpDen[4], csSpJonDen[4], csZgNum[4], csZgDen[4];
  double Z2ll(dZ2ll), W2e(dW2e), W2mu(dW2mu), W2tau(dW2tau);
  // double WZbr[4]={dWZ23e, dWZ22e1mu, dWZ21e2mu, dWZ23mu};
  double WZbr[4]={0.03363*0.1075,0.03363*0.1057,0.03366*0.1075,0.03366*0.1057};
  double Nsig[4];
  double  systematicErrorZagreb[4]={0,0,0,0};
  double  systematicErrorSpanish[4]={0,0,0,0};
  double sys[4]={0.05, 0.047, 0.054, 0.052};  ///all together Spanish numbers
  double luminosity (LUMI);
  
  //  double WZ23lnu=3*Z2ll*(W2e+W2mu+W2tau);
  double WZ23lnu=3.*0.033658*(0.1125+0.1075+0.1057);
  TString names[4]={"3e", "2e1mu","1e2mu", "3mu"};

  for (int i=0; i<4; i++){
    Nsig[i]=(N_good[i]-NZgamma[i]-NWV[i]-NVVV[i]-NZZ[i]);
    csSpNum[i]=Nsig[i];

    csSpDen[i]=(AxeS[i]*luminosity*WZ23lnu);
    csSpJonDen[i]=(AxeSJ[i]*luminosity);
    crossSectionSpanish[i]=csSpNum[i]/csSpDen[i];
    crossSectionSpanishJonatan[i]=csSpNum[i]/csSpJonDen[i];
    

    csZgNum[i]=(1-tauFactor[i])*Nsig[i];
    csZgDen[i]=(AxeZ[i]*luminosity*WZbr[i]);
    crossSectionZagreb[i] =csZgNum[i]/csZgDen[i];


    errorCsSpanish[i]=sqrt(pow((sN_good[i]/csSpDen[i]),2)+pow(((Nsig[i]*luminosity*WZ23lnu*sAxeS[i])/(csSpDen[i]*csSpDen[i])), 2));
    errorCsZagreb[i]=sqrt(pow(((sN_good[i]*(1-tauFactor[i]))/csZgDen[i]),2)+pow((((1-tauFactor[i])*Nsig[i]*luminosity*WZbr[i]*sAxeZ[i])/(csZgDen[i]*csZgDen[i])), 2));

    //    errorCsZgNum[i]=sqrt(pow((N_good[i]-NZgamma[i]-NWV[i]-NVVV[i]-NZZ[i])*stauFactor[i],2)+pow((1-tauFactor[i])*sN_good[i],2)+ pow((1-tauFactor[i])*sNZZ[i],2)+pow((1-tauFactor[i])*sNZgamma[i],2) + pow((1-tauFactor[i])*sNWV[i],2) + pow((1-tauFactor[i])*sNVVV[i],2));
    errorCsZgDen[i]=luminosity*sAxeZ[i];


    //spanish cross section all derivation
    double denom=AxeS[i]*luminosity;

    double errorTest=sqrt(pow(sN_good[i]/(denom),2) + pow(sNZgamma[i]/denom,2)+ pow(NWV[i]/denom,2)+ pow(NVVV[i]/denom,2) + pow(NZZ[i]/denom,2)+ pow((N_good[i]-NZgamma[i]-NWV[i]-NVVV[i]-NZZ[i])*sAxeS[i]/(AxeS[i]*AxeS[i]*luminosity),2));
    std::cout<<"error test:"<<errorTest<<std::endl;
    
    //systematics
    systematicErrorZagreb[i]=sys[i]*crossSectionZagreb[i];
    systematicErrorSpanish[i]=sys[i]*crossSectionSpanish[i];
  }



  /////BLUE METHOD////////

  //spanish numbers- systematics

  double qcdScale[4]={qcdScale3e/100, qcdScale2e1mu/100, qcdScale1e2mu/100, qcdScale3mu/100};
  double PDFsys[4]={PDF3e/100, PDF2e1mu/100, PDF1e2mu/100, PDF3mu/100};
  double leptTrgEff[4] ={LepTrgEff3e/100, LepTrgEff2e1mu/100, LepTrgEff1e2mu/100, LepTrgEff3mu/100};
  double Etsys[4] = {Etmiss3e/100, Etmiss2e1mu/100, Etmiss1e2mu/100, Etmiss3mu/100};
  double muMomScale[4] = {muMomScale3e/100, muMomScale2e1mu/100, muMomScale1e2mu/100, muMomScale3mu/100};
  double elEnScale[4]= {elEnScale3e/100, elEnScale2e1mu/100, elEnScale1e2mu/100, elEnScale3mu/100};
  double pileUpsys[4]= {pileUp3e/100, pileUp2e1mu/100, pileUp1e2mu/100, pileUp3mu/100};
  double ZZcs[4]= {ZZcs3e/100, ZZcs2e1mu/100, ZZcs1e2mu/100, ZZcs3mu/100};
  double Zgammacs[4]= {Zgammacs3e/100, Zgammacs2e1mu/100, Zgammacs1e2mu/100, Zgammacs3mu/100};
  double dataDrivensys[4] = {dataDriven3e/100, dataDriven2e1mu/100, dataDriven1e2mu/100, dataDriven3mu/100};
  double bckgSys[4] = {back3e/100, back2e1mu/100, back1e2mu/100, back3mu/100};
  
  //compose Error Matrix [i = column*4+row] [0-3][0-3] indices
//matrix is symmetric
//channel legend: 0 = 3e, 1 = 2e1mu, 2 = 2mu1e,  3 = 3mu
  
  Double_t elmZagreb[16];
  for (size_t elm=0;elm<16;elm++) elmZagreb[elm]=0;
  
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
  elmZagreb[0]=pow(systematicErrorZagreb[0],2) + pow(errorCsZagreb[0],2);
  elmZagreb[5]=pow(systematicErrorZagreb[1],2) + pow(errorCsZagreb[1],2);
  elmZagreb[10]=pow(systematicErrorZagreb[2],2) + pow(errorCsZagreb[2],2);
  elmZagreb[15]=pow(systematicErrorZagreb[3],2) + pow(errorCsZagreb[3],2);
  
  //matrix is symmetric
  //channels 0 and 1 : [ trigger&scaleFactors, electronEnergyScale, muonEnergyScale]
  elmZagreb[4]=elmZagreb[1]= crossSectionZagreb[0]*crossSectionZagreb[1]* (commonSys[0][1] + elEnScale[0]*elEnScale[1] + muMomScale[0]*muMomScale[1]+
									   leptTrgEff[0]*(sqrt(2/3)*leptTrgEff[1]));
  /*
  std::cout<<"debugging"<<commonSys[0][1]<<std::endl;

  std::cout<<"debugging"<<elEnScale[0]*elEnScale[1]<<std::endl;
  std::cout<<"debugging"<<muMomScale[0]*muMomScale[1]<<std::endl;
  std::cout<<"debugging"<<leptTrgEff[0]*sqrt<<std::endl;
  */

  //channels 0 and 2: 
  elmZagreb[8]=elmZagreb[2]=  crossSectionZagreb[0]*crossSectionZagreb[2]* (commonSys[0][2] + elEnScale[0]*elEnScale[2] + muMomScale[0]*muMomScale[2]+
									   leptTrgEff[0]*(sqrt(1/3)*leptTrgEff[2]));

  //channels 0 and 3
  elmZagreb[12]=elmZagreb[3]= crossSectionZagreb[0]*crossSectionZagreb[3]* (commonSys[0][3] + elEnScale[0]*elEnScale[3] + muMomScale[0]*muMomScale[3]);

  //channels 1 and 2
  elmZagreb[9]=elmZagreb[6] =crossSectionZagreb[1]*crossSectionZagreb[2]* (commonSys[1][2] +elEnScale[1]*elEnScale[2] + muMomScale[1]*muMomScale[2]+
									   2*sqrt(1/3)*leptTrgEff[1]*sqrt(2/3)*leptTrgEff[2]);
  //channels 1 and 3
  elmZagreb[13]=elmZagreb[7] =crossSectionZagreb[1]*crossSectionZagreb[3] *(commonSys[1][3] + elEnScale[1]*elEnScale[3] + muMomScale[1]*muMomScale[3]+
									   leptTrgEff[1]*(sqrt(1/3)*leptTrgEff[3]));
  //channels 2 and 3
  elmZagreb[14]=elmZagreb[14] =crossSectionZagreb[2]*crossSectionZagreb[3] *(commonSys[2][3] + elEnScale[2]*elEnScale[3] + muMomScale[2]*muMomScale[3]+
									     leptTrgEff[2]*(sqrt(2/3)*leptTrgEff[3]));
 
  std::cout<<"Error matrix finished!!"<<std::endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      //      std::cout<<"!!"<<std::endl;
      std:: cout << "       " << elmZagreb[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl; 
  
//initialize matrix for inversion
  TMatrixD errMat(4,4,elmZagreb);
  TMatrixD errMatInv(4,4);
  
  TMatrixD errMatCopy(errMat);
  
  //invert !
  errMatInv = errMat.Invert();
  
  
  Double_t *mRefTest = errMat.GetMatrixArray();
  
  std::cout << endl << " Matrix Inverse:" << endl;
  for (int emi=0;emi<4;emi++) {
    for (int emj=0;emj<4;emj++) {
      std::cout << "       " << mRefTest[emj*4+emi] << "       ";
    }
    std::cout << endl;
  }
  std::cout << endl;

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

  std::cout << "al0 " << alphaCH[0] << " al1 " << alphaCH[1] << " al2 " << alphaCH[2] << " al3 " << alphaCH[3] << endl;
  std::cout << "consistency check:" << alphaCH[0]+alphaCH[1]+alphaCH[2]+alphaCH[3] <<endl;
  std::cout << endl;
  double final_Xsec = alphaCH[0]*crossSectionZagreb[0] + alphaCH[1]*crossSectionZagreb[1] + alphaCH[2]*crossSectionZagreb[2] + alphaCH[3]*crossSectionZagreb[3];

  Double_t combined_error=0;
  Double_t *copyRef = errMatCopy.GetMatrixArray();
  for (int ier=0;ier<4;ier++)
    for (int jer=0;jer<4;jer++) {
      combined_error+=alphaCH[ier]*alphaCH[jer]*copyRef[ier*4+jer];
    }
  Double_t stat_err_tot2 =  pow(alphaCH[0]*errorCsZagreb[0],2) + pow(alphaCH[1]*errorCsZagreb[1],2)
    +pow(alphaCH[2]*errorCsZagreb[2],2) + pow(alphaCH[3]*errorCsZagreb[3],2);
  
  Double_t stat_err_tot = sqrt(stat_err_tot2);
  Double_t syst_err_tot = sqrt (combined_error - stat_err_tot2);
  
  combined_error=sqrt(combined_error);


  std::cout << "combined sigma(WZ) = " << final_Xsec << " +- " << stat_err_tot << "(stat.) +- " << syst_err_tot << "(syst) +- "
	    << 0.044*final_Xsec << " (lumi) pb " << endl;
  
  std::cout << endl;
 //systematics numbers
  /********************************************************************************************/
//spanish numbers for now
  
  
  
  /*
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
  
 
*/


  ////**********outputs**********************
  std::cout<<"*********Spanish cross section (inclusive)***********"<<std::endl;
  for (int i1=0; i1<4; i1++){
    std::cout<<i1<<" : "<<crossSectionSpanish[i1] <<"+/-"<<errorCsSpanish[i1]<<"+/-"<<systematicErrorSpanish[i1]  <<std::endl;
  }
  std::cout<<"*****************************************"<<std::endl;
  std::cout<<"*********Zagreb cross section(inclusive)************"<<std::endl;
  for (int i2=0; i2<4; i2++){
    std::cout<<i2<<" : "<<crossSectionZagreb[i2] <<"+/-"<<errorCsZagreb[i2] <<"+/-"<<systematicErrorZagreb[i2]<<std::endl;
  }
  
  ////**************latexOutput
  if (latexOutput){

    std::cout<<"latex output"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"Acceptance times efficiency"<<std::endl;
    std::cout<<"channel & AxeSpanish  & AxeZagreb \\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"  & $"<<AxeS[a]<<"\\pm "<<sAxeS[a]<<"$ & $"<<AxeZ[a]<<"\\pm "<<sAxeZ[a]<<"$\\\\"<<std::endl;  
      std::cout<<"\\hline"<<std::endl;
    }
    std::cout<<"*****************"<<std::endl;
    /*

    std::cout<<"Matrix method"<<std::endl;
    std::cout<<"channel & Ngood  & AxeNfake"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<a<<"  & $"<<AxeS[a]<<"\pm "<<sAxeS[a]<<"$ & $"<<AxeZ[a]<<"/pm "<<sAxeZ[a]<<"\\\\"<<std::endl;  
      std::cout<<"\\hline"<<std::endl;
    }
    */
    std::cout<<"*****************"<<std::endl;
    std::cout<<"Matrix method results"<<std::endl;
    std::cout<<"channel & Ndata & Ngood & Nfake\\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
    for (int mm=0; mm<4; mm++){
      std::cout<<names[mm]<<" & " <<N_data[mm]<<" &  $"<<N_good[mm]<<"\\pm "<<sN_good[mm] <<"$ & $"<<N_fake[mm]<<"\\pm"<<sN_fake[mm]<<"$ \\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
    std::cout<<"*****************"<<std::endl;
    std::cout<<"MC  yields"<<std::endl;
    /*
    std::cout<<"channel &  ZZ & Zgamma & WV & VVV\\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"$  & $"<<NZZ[a]<<"\\pm "<<sNZZ[a]<<"$ & $"<<NZgamma[a]<<"\\pm "<<sNZgamma[a]<<"$ & $"<<NWV[a]<<"\\pm "<<sNWV[a]<<"$ & $"<<NVVV[a]<<"\\pm "<<sNVVV[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
    */
    
    std::cout<<" & 3e & 2e1mu & 1e2mu & 3mu \\\\"<<std::endl;
    std::cout<<"ZZ & $"<<NZZ[0]<<"\\pm "<<sNZZ[0]<<"$ & $"<<NZZ[1]<<"\\pm "<<sNZZ[1]<<"$ & $"<<NZZ[2]<<"\\pm "<<sNZZ[2]<<"$ & $"<<NZZ[3]<<"\\pm "<<sNZZ[3]<<"$ \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
    std::cout<<"Zgamma & $"<<NZgamma[0]<<"\\pm "<<sNZgamma[0]<<"$ & $"<<NZgamma[1]<<"\\pm "<<sNZgamma[1]<<"$ & $"<<NZgamma[2]<<"\\pm "<<sNZgamma[2]<<"$ & $"<<NZgamma[3]<<"\\pm "<<sNZgamma[3]<<"$ \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
    std::cout<<"WV & $"<<NWV[0]<<"\\pm "<<sNWV[0]<<"$ & $"<<NWV[1]<<"\\pm "<<sNWV[1]<<"$ & $"<<NWV[2]<<"\\pm "<<sNWV[2]<<"$ & $"<<NWV[3]<<"\\pm "<<sNWV[3]<<"$ \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
    std::cout<<"VVV & $"<<NVVV[0]<<"\\pm "<<sNVVV[0]<<"$ & $"<<NVVV[1]<<"\\pm "<<sNVVV[1]<<"$ & $"<<NVVV[2]<<"\\pm "<<sNVVV[2]<<"$ & $"<<NVVV[3]<<"\\pm "<<sNVVV[3]<<"$ \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;
    //    std::cout<<"WZ & $"<<NWZ[0]<<"\\pm "<<sNWZ[0]<<"$ &"<<NWZ[1]<<"\\pm "<<sNWZ[1]<<"$ &"<<NWZ[2]<<"\\pm "<<sNWZ[2]<<"$ &"<<NWZ[3]<<"\\pm "<<sNWZ[3]<<"$ \\\\"<<std::endl;

    std::cout<<"*****************"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"tau factor"<<std::endl;
    std::cout<<"channel &  tauFactor\\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"  & $"<<tauFactor[a]<<"\\pm "<<stauFactor[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
    std::cout<<"*****************"<<std::endl;
    std::cout<<"*****************"<<std::endl;
    std::cout<<"crossSection"<<std::endl;
    std::cout<<"channel &  crossSection Zg & crossSection Sp\\\\"<<std::endl;
    for (int a=0; a<4; a++){
      std::cout<<names[a]<<"  & $"<<crossSectionZagreb[a]<<"\\pm "<<errorCsZagreb[a]<<"\\pm "<<systematicErrorZagreb[a]<<"$  & $"<<crossSectionSpanish[a]<<"\\pm "<<errorCsSpanish[a]<<"\\pm"<<systematicErrorSpanish[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
  }

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
}
