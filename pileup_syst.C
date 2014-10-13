#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

void pileup_syst(){
  
  bool latexOutput(true);
  bool latexOutput2(false);

  double csNominal[4], csUp[4], csDown[4], tauNominal[4], tauUp[4], tauDown[4], axeNominal[4], axeUp[4], axeDown[4];
  fstream myfileUp("numWVup.txt", std::ios_base::in);
  fstream myfileNominal("numNominal.txt", std::ios_base::in);
  fstream myfileDown("numWVdown.txt", std::ios_base::in);

   myfileNominal>>csNominal[0]>>csNominal[1]>>csNominal[2]>>csNominal[3]>>tauNominal[0]>>tauNominal[1]>>tauNominal[2]>>tauNominal[3]>>axeNominal[0]>>axeNominal[1]>>axeNominal[2]>>axeNominal[3];

  myfileUp>>csUp[0]>>csUp[1]>>csUp[2]>>csUp[3]>>tauUp[0]>>tauUp[1]>>tauUp[2]>>tauUp[3]>>axeUp[0]>>axeUp[1]>>axeUp[2]>>axeUp[3];

  myfileDown>>csDown[0]>>csDown[1]>>csDown[2]>>csDown[3]>>tauDown[0]>>tauDown[1]>>tauDown[2]>>tauDown[3]>>axeDown[0]>>axeDown[1]>>axeDown[2]>>axeDown[3];

  //same thing for cross section
  
  double ldCs[4], ldTau[4],ldAxe[4];
  double csSyst[4], tauSyst[4], axeSyst[4];
  for (int scale=0; scale<4; scale++){
    //std::cout<<csNominal[scale]<<","<<csUp[scale]<<","<<csDown[scale]<<std::endl;
    ldCs[scale]=std::max(fabs(csNominal[scale]-csUp[scale]), fabs(csNominal[scale]-csDown[scale]));
    
    ldTau[scale]=std::max(fabs(tauNominal[scale]-tauUp[scale]), fabs(tauNominal[scale]-tauDown[scale]));
    ldAxe[scale]=std::max(fabs(axeNominal[scale]-axeUp[scale]), fabs(axeNominal[scale]-axeDown[scale]));
    csSyst[scale]=ldCs[scale]/csNominal[scale];
    tauSyst[scale]=ldTau[scale]/tauNominal[scale];
    axeSyst[scale]=ldAxe[scale]/axeNominal[scale];
  }
  std::cout<<"Cross section:"<<std::endl;
  for (int ispisEl=0; ispisEl<4; ispisEl++){
    std::cout<<ispisEl<<":"<<csSyst[ispisEl]<<std::endl;
  }
  std::cout<<"*******"<<std::endl;
  std::cout<<"tau factor:"<<std::endl;
  for (int ispisMu=0; ispisMu<4; ispisMu++){
    std::cout<<ispisMu<<":"<<tauSyst[ispisMu]<<std::endl;
  }
  std::cout<<"*******"<<std::endl;
  std::cout<<"Axe: "<<std::endl;
  for (int ispisMu2=0; ispisMu2<4; ispisMu2++){
    std::cout<<ispisMu2<<":"<<axeSyst[ispisMu2]<<std::endl;
  }
  TString labels[4]={"3e", "2e1mu", "1e2mu", "3mu"};
  if (latexOutput){
    std::cout<<"Cross section:"<<std::endl;
    std::cout<<"\\begin{center}"<<std::endl;
    std::cout<<"\\begin{tabular}{cc}"<<std::endl;
    std::cout<<"channel & pileupSyst\\\\"<<std::endl;
    for (int latex=0; latex<4; latex++){
      std::cout<<"\\hline"<<std::endl;
      std::cout<<labels[latex]<<" & "<<csSyst[latex]<<"\\\\"<<std::endl;
	} 
    std::cout<<"\\hline"<<std::endl;
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{center}"<<std::endl;

    std::cout<<"++++++ axe:"<<std::endl;
    std::cout<<"\\begin{center}"<<std::endl;
    std::cout<<"\\begin{tabular}{cc}"<<std::endl;
    std::cout<<"channel & pileupSyst\\\\"<<std::endl;
    for (int latex=0; latex<4; latex++){
      std::cout<<"\\hline"<<std::endl;
      std::cout<<labels[latex]<<" & "<<axeSyst[latex]<<"\\\\"<<std::endl;
	} 
    std::cout<<"\\hline"<<std::endl;
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{center}"<<std::endl;

    std::cout<<"++++++ tau:"<<std::endl;
    std::cout<<"\\begin{center}"<<std::endl;
    std::cout<<"\\begin{tabular}{cc}"<<std::endl;
    std::cout<<"channel & pileupSyst\\\\"<<std::endl;
    for (int latex=0; latex<4; latex++){
      std::cout<<"\\hline"<<std::endl;
      std::cout<<labels[latex]<<" & "<<tauSyst[latex]<<"\\\\"<<std::endl;
	} 
    std::cout<<"\\hline"<<std::endl;
    std::cout<<"\\end{tabular}"<<std::endl;
    std::cout<<"\\end{center}"<<std::endl;


  }
  if (latexOutput2){
    std::cout<<" & "<<csSyst[0]<<" & "<<csSyst[1]<<" & "<<csSyst[2]<<" & "<<csSyst[3]<<"\\\\"<<std::endl;
  }
}
