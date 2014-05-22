#include "numbers.h"
#include "numMC.h"
#include "numMM_met.h"
//#include "numMM.h"
//#include "numGEN.h"
//#include "numData.h"
#include "numGEN_met.h"
#include "numData_met.h"

void crossSection()
{
  //  bool latexOutput(false);
  bool latexOutput(true);
  //0: eee, 1:eem, 2:emm, 3:mmm

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
  double luminosity(LUMI);
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
    

    //    std::cout<<"ZG "<<(1-tauFactor[i])/(AxeZ[i]*WZbr[i])<<std::endl;
    // std::cout<<"SP "<<1/(AxeSJ[i]*WZ23lnu)<<std::endl;
    }
  /*
  std::cout<<"*********Spanish cross section***********"<<std::endl;
  for (int i1=0; i1<4; i1++){
    std::cout<<i1<<" : "<<crossSectionSpanish[i1] <<"+/-"<<errorCsSpanish[i1]  <<std::endl;
  }
  std::cout<<"*****************************************"<<std::endl;
  std::cout<<"*********Zagreb cross section************"<<std::endl;
  for (int i2=0; i2<4; i2++){
    std::cout<<i2<<" : "<<crossSectionZagreb[i2] <<"+/-"<<errorCsZagreb[i2] <<std::endl;
  }
  */
  /*
  std::cout<<"*********Spanish cross section (inclusive)----Jonatan***********"<<std::endl;
  for (int i1=0; i1<4; i1++){
    std::cout<<i1<<" : "<<crossSectionSpanishJonatan[i1]/WZ23lnu <<"+/-"<<errorCsSpanish[i1]  <<std::endl;
  }
  std::cout<<"*****************************************"<<std::endl;
  */
  std::cout<<"*********Spanish cross section (inclusive)***********"<<std::endl;
  for (int i1=0; i1<4; i1++){
    std::cout<<i1<<" : "<<crossSectionSpanish[i1] <<"+/-"<<errorCsSpanish[i1]  <<std::endl;
  }
  std::cout<<"*****************************************"<<std::endl;
  std::cout<<"*********Zagreb cross section(inclusive)************"<<std::endl;
  for (int i2=0; i2<4; i2++){
    std::cout<<i2<<" : "<<crossSectionZagreb[i2] <<"+/-"<<errorCsZagreb[i2] <<std::endl;
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
      std::cout<<names[a]<<"  & $"<<crossSectionZagreb[a]<<"\\pm "<<errorCsZagreb[a]<<"$  & $"<<crossSectionSpanish[a]<<"\\pm "<<errorCsSpanish[a]<<"$\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
    }
  }
}
