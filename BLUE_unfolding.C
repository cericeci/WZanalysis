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
#include <set>

#include "UnfoldingHistogramFactory.h"

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

  const int nChannels(4);
  TH1D * h_qcdScale[nChannels];
  TH1D * h_PDFsys[nChannels];
  TH1D * h_leptTrgEff_el[nChannels];
  TH1D * h_leptTrgEff_mu[nChannels];
  TH1D * h_Etsys[nChannels];
  TH1D * h_muMomScale[nChannels];
  TH1D * h_elEnScale[nChannels];
  TH1D * h_pileupSys[nChannels];
  TH1D * h_ZZxs[nChannels];
  TH1D * h_Zgammaxs[nChannels];
  TH1D * h_dataDrivensys[nChannels];
  TH1D * h_dataDrivensys_el[nChannels];
  TH1D * h_dataDrivensys_mu[nChannels];
  TH1D * h_bckgSys[nChannels];
  TH1D * h_JESsys[nChannels];
  TH1D * h_JERsys[nChannels];
  TH1D * h_crossSection[nChannels];
  TH1D * h_crossSection_final[nChannels];
  TH1D * h_crossSection_diff[nChannels];
  TH1D * h_totalSyst[nChannels];
  TH1D * h_totalSyst_diff[nChannels];
  TH1D * h_totalStat[nChannels];
  TH1D * h_totalStat_diff[nChannels];
  TH1D * h_unfSyst[nChannels];
  TH1D * h_lumi[nChannels];

  TH1D * h_crossSection_combination;
  TH1D * h_crossSection_comb_diff;
  TH1D * h_combStat;
  TH1D * h_combSyst;
  TH1D * h_combLumi;
  TH1D * h_combStat_diff;
  TH1D * h_combSyst_diff;
  double crossSection[4]={0,0,0,0};
  

  //ovo treba deklarirat!!!!!!!!!
  for (int fill=0; fill<nChannels; fill++){
    std::ostringstream totalSysName, totalSysDiffName, totalStatName, totalStatDiffName;
    totalSysName<<"h_totalSyst_"<<fill;
    totalSysDiffName<<"h_totalSyst_diff_"<<fill;
    totalStatName<<"h_totalStat_"<<fill;
    totalStatDiffName<<"h_totalStat_diff_"<<fill;
    if (variable=="Njets"){
      h_totalSyst[fill]=UnfoldingHistogramFactory::createNjetsHistogram(totalSysName.str().c_str(), totalSysName.str().c_str());
      h_totalSyst_diff[fill]=UnfoldingHistogramFactory::createNjetsHistogram(totalSysDiffName.str().c_str(), totalSysDiffName.str().c_str());
      h_totalStat[fill]=UnfoldingHistogramFactory::createNjetsHistogram(totalStatName.str().c_str(), totalStatName.str().c_str());
      h_totalStat_diff[fill]=UnfoldingHistogramFactory::createNjetsHistogram(totalStatDiffName.str().c_str(), totalStatDiffName.str().c_str());
    }
    if (variable=="LeadingJetPt"){
      h_totalSyst[fill]=UnfoldingHistogramFactory::createLeadingJetHistogram(totalSysName.str().c_str(), totalSysName.str().c_str());
      h_totalSyst_diff[fill]=UnfoldingHistogramFactory::createLeadingJetHistogram(totalSysDiffName.str().c_str(), totalSysDiffName.str().c_str());
      h_totalStat[fill]=UnfoldingHistogramFactory::createLeadingJetHistogram(totalStatName.str().c_str(), totalStatName.str().c_str());
      h_totalStat_diff[fill]=UnfoldingHistogramFactory::createLeadingJetHistogram(totalStatDiffName.str().c_str(), totalStatDiffName.str().c_str());
    }
    if (variable=="Zpt"){
      h_totalSyst[fill]=UnfoldingHistogramFactory::createZPtHistogram(totalSysName.str().c_str(), totalSysName.str().c_str());
      h_totalSyst_diff[fill]=UnfoldingHistogramFactory::createZPtHistogram(totalSysDiffName.str().c_str(), totalSysDiffName.str().c_str());
      h_totalStat[fill]=UnfoldingHistogramFactory::createZPtHistogram(totalStatName.str().c_str(), totalStatName.str().c_str());
      h_totalStat_diff[fill]=UnfoldingHistogramFactory::createZPtHistogram(totalStatDiffName.str().c_str(), totalStatDiffName.str().c_str());
    }
  }

  if (variable=="Njets"){
    h_crossSection_combination= UnfoldingHistogramFactory::createNjetsHistogram("h_xs_comb", "h_xs_comb");
    h_crossSection_comb_diff= UnfoldingHistogramFactory::createNjetsHistogram("h_xs_comb_diff", "h_xs_comb_diff");
    h_combStat = UnfoldingHistogramFactory::createNjetsHistogram("h_combStat", "h_combStat"); 
    h_combSyst = UnfoldingHistogramFactory::createNjetsHistogram("h_combSyst", "h_combSyst");
    h_combLumi = UnfoldingHistogramFactory::createNjetsHistogram("h_combLumi", "h_combLumi");
    h_combStat_diff = UnfoldingHistogramFactory::createNjetsHistogram("h_combStat_diff", "h_combStat_diff"); 
    h_combSyst_diff = UnfoldingHistogramFactory::createNjetsHistogram("h_combSyst_diff", "h_combSyst_diff");
  }
  else if (variable=="LeadingJetPt"){
    h_crossSection_combination= UnfoldingHistogramFactory::createLeadingJetHistogram("h_xs_comb", "h_xs_comb");
    h_crossSection_comb_diff= UnfoldingHistogramFactory::createLeadingJetHistogram("h_xs_comb_diff", "h_xs_comb_diff");
    h_combStat = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combStat", "h_combStat"); 
    h_combSyst = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combSyst", "h_combSyst");
    h_combLumi = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combLumi", "h_combLumi");
    h_combStat_diff = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combStat_diff", "h_combStat_diff"); 
    h_combSyst_diff = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combSyst_diff", "h_combSyst_diff");
  }
  else if (variable=="Zpt"){
    h_crossSection_combination= UnfoldingHistogramFactory::createZPtHistogram("h_xs_comb", "h_xs_comb");
    h_crossSection_comb_diff= UnfoldingHistogramFactory::createZPtHistogram("h_xs_comb_diff", "h_xs_comb_diff");
    h_combStat = UnfoldingHistogramFactory::createZPtHistogram("h_combStat", "h_combStat"); 
    h_combSyst = UnfoldingHistogramFactory::createZPtHistogram("h_combSyst", "h_combSyst");
    h_combLumi = UnfoldingHistogramFactory::createZPtHistogram("h_combLumi", "h_combLumi");
    h_combStat_diff = UnfoldingHistogramFactory::createZPtHistogram("h_combStat_diff", "h_combStat_diff"); 
    h_combSyst_diff = UnfoldingHistogramFactory::createZPtHistogram("h_combSyst_diff", "h_combSyst_diff");
  }
  else
    std::cout<<"UNKNOWN VARIABLE!!!"<<std::endl;


  for (int hist=0; hist<nChannels; hist++){
    std::ostringstream qcdScaleName,PDFsysName,leptTrgEffName_el, leptTrgEffName_mu, EtsysName, 
      muMomScaleName, elEnScaleName, pileupSysName, ZZxsName, ZgammaxsName, 
      dataDrivenName, bckgSysName, xsName, xsNameFinal, xsNameDiff, JESsysName, JERsysName, lumiName,
      dataDrivenNameEl, dataDrivenNameMu, unfSystName;
    qcdScaleName<<"h_qcdScale_"<<hist;
    PDFsysName<<"h_PDFsys_"<<hist;
    leptTrgEffName_el<<"ele_SF_"<<hist;
    leptTrgEffName_mu<<"mu_SF_"<<hist;
    EtsysName<<"h_Etsys_"<<hist;
    muMomScaleName<<"h_muMomScale_"<<hist;
    elEnScaleName<<"h_elEnScale_"<<hist;
    pileupSysName<<"h_pileupSys_"<<hist;
    ZZxsName<<"h_ZZxs_"<<hist;
    ZgammaxsName<<"h_Zgammaxs_"<<hist;
    dataDrivenName<<"h_dataDrivensys_"<<hist;
    dataDrivenNameEl<<"h_dataDrivenElsys_"<<hist;
    dataDrivenNameMu<<"h_dataDrivenMusys_"<<hist;
    bckgSysName<<"h_bckgSys_"<<hist;
    xsName<<"h_crossSection_"<<hist;
    xsNameFinal<<"h_crossSection_inclusive"<<hist;
    xsNameDiff<<"h_crossSection_incl_diff"<<hist;
    JESsysName<<"h_JESsys_"<<hist;
    JERsysName<<"h_JERsys_"<<hist;
    lumiName<<"h_lumi_"<<hist;
    unfSystName<<"h_unfSyst_"<<hist;

    h_qcdScale[hist] = (TH1D*) (finput ->Get(qcdScaleName.str().c_str())->Clone(qcdScaleName.str().c_str()));
    h_lumi[hist] = (TH1D*) (finput ->Get(lumiName.str().c_str())->Clone(lumiName.str().c_str()));
    h_PDFsys[hist] = (TH1D*) (finput ->Get(PDFsysName.str().c_str())->Clone(PDFsysName.str().c_str()));
    h_leptTrgEff_el[hist] = (TH1D*) (finput ->Get(leptTrgEffName_el.str().c_str())->Clone(leptTrgEffName_el.str().c_str()));
    h_leptTrgEff_mu[hist] = (TH1D*) (finput ->Get(leptTrgEffName_mu.str().c_str())->Clone(leptTrgEffName_mu.str().c_str()));
    h_Etsys[hist] = (TH1D*) (finput ->Get(EtsysName.str().c_str())->Clone(EtsysName.str().c_str()));
    h_muMomScale[hist] = (TH1D*) (finput ->Get(muMomScaleName.str().c_str())->Clone(muMomScaleName.str().c_str()));
    h_elEnScale[hist] = (TH1D*) (finput ->Get(elEnScaleName.str().c_str())->Clone(elEnScaleName.str().c_str()));
    h_pileupSys[hist] = (TH1D*) (finput ->Get(pileupSysName.str().c_str())->Clone(pileupSysName.str().c_str()));
    h_ZZxs[hist] = (TH1D*) (finput ->Get(ZZxsName.str().c_str())->Clone(ZZxsName.str().c_str()));
    h_Zgammaxs[hist] = (TH1D*) (finput ->Get(ZgammaxsName.str().c_str())->Clone(ZgammaxsName.str().c_str()));
    h_dataDrivensys[hist] = (TH1D*) (finput ->Get(dataDrivenName.str().c_str())->Clone(dataDrivenName.str().c_str()));
    h_dataDrivensys_el[hist] = (TH1D*) (finput ->Get(dataDrivenNameEl.str().c_str())->Clone(dataDrivenNameEl.str().c_str()));
    h_dataDrivensys_mu[hist] = (TH1D*) (finput ->Get(dataDrivenNameMu.str().c_str())->Clone(dataDrivenNameMu.str().c_str()));
    h_bckgSys[hist] = (TH1D*) (finput ->Get(bckgSysName.str().c_str())->Clone(bckgSysName.str().c_str()));
    h_crossSection[hist] = (TH1D*) (finput ->Get(xsName.str().c_str())->Clone(xsName.str().c_str()));
    h_crossSection_final[hist] = (TH1D*) (finput ->Get(xsNameFinal.str().c_str())->Clone(xsNameFinal.str().c_str()));
    h_crossSection_diff[hist] = (TH1D*) (finput ->Get(xsNameFinal.str().c_str())->Clone(xsNameDiff.str().c_str()));
    h_unfSyst[hist]= (TH1D*)(finput->Get(unfSystName.str().c_str())->Clone(unfSystName.str().c_str()));

    if (variable!="Zpt"){
      std::cout<<variable<<"  JET SYSTEMATICS :D"<<std::endl;
      h_JESsys[hist] = (TH1D*) (finput->Get(JESsysName.str().c_str())->Clone(JESsysName.str().c_str()));
      h_JERsys[hist] = (TH1D*) (finput->Get(JERsysName.str().c_str())->Clone(JERsysName.str().c_str()));
    }
    
  }

  //looop over each bin 
  
  for (int bin=1; bin< (h_crossSection[0]->GetNbinsX()+1); bin++){
    double elements[16];
    for (int el=0; el<16; el++) elements[el]=0;
    double statisticError[nChannels];
    double systematicError2[nChannels];
    double systematicError[nChannels];
    double lumiError[nChannels];
    //statistic and systematic errors:
    for (int nCh=0; nCh<nChannels; nCh++){
      //      statisticError[nCh]=h_crossSection[nCh]->GetBinError(bin);
      statisticError[nCh]=h_crossSection_final[nCh]->GetBinError(bin);
      systematicError2[nCh]=(pow(h_qcdScale[nCh]->GetBinContent(bin),2)+
			     pow(h_PDFsys[nCh]->GetBinContent(bin), 2)+
			     pow(h_leptTrgEff_el[nCh]->GetBinContent(bin),2)+
			     pow(h_leptTrgEff_mu[nCh]->GetBinContent(bin),2)+
			     pow(h_Etsys[nCh]->GetBinContent(bin), 2)+
			     pow(h_muMomScale[nCh]->GetBinContent(bin),2)+
			     pow(h_elEnScale[nCh]->GetBinContent(bin),2)+
			     pow(h_pileupSys[nCh]->GetBinContent(bin),2)+
			     pow(h_ZZxs[nCh]->GetBinContent(bin),2)+
			     pow(h_Zgammaxs[nCh]->GetBinContent(bin),2)+
			     pow(h_dataDrivensys_el[nCh]->GetBinContent(bin),2)+
			     pow(h_dataDrivensys_mu[nCh]->GetBinContent(bin),2)+
			     pow(h_unfSyst[nCh]->GetBinContent(bin),2)+
			     pow(h_bckgSys[nCh]->GetBinContent(bin),2));
      if (variable!="Zpt"){
	systematicError2[nCh]+= pow (h_JESsys[nCh]->GetBinContent(bin),2)+
	  pow (h_JERsys[nCh]->GetBinContent(bin),2);
      }

      h_totalStat[nCh]->SetBinContent(bin, statisticError[nCh]);
      //      systematicError[nCh]= (sqrt(systematicError2[nCh]))*(h_crossSection[nCh]->GetBinContent(bin));
      systematicError[nCh]= (sqrt(systematicError2[nCh]))*(h_crossSection_final[nCh]->GetBinContent(bin));
      lumiError[nCh]=(h_lumi[nCh]->GetBinContent(bin))*(h_crossSection_final[nCh]->GetBinContent(bin));
      h_totalSyst[nCh]->SetBinContent(bin, (systematicError[nCh]));
      std::cout<<bin<<" , "<<systematicError[nCh]<<" , "<<h_crossSection_final[nCh]->GetBinContent(bin)<<std::endl;
    }
    //common elements
    double commonSys[4][4];
    double lumiSys[4][4];

    for (int cha=0; cha<4; cha++){
      for (int chb=0; chb<4; chb++){
	commonSys[cha][chb]=(h_Etsys[cha]->GetBinContent(bin))*(h_Etsys[chb]->GetBinContent(bin))
	  + (h_pileupSys[cha]->GetBinContent(bin))*(h_pileupSys[chb]->GetBinContent(bin))
	  + (h_PDFsys[cha]->GetBinContent(bin))*(h_PDFsys[chb]->GetBinContent(bin))
	  + (h_qcdScale[cha]->GetBinContent(bin))*(h_qcdScale[chb]->GetBinContent(bin))
	  + (h_bckgSys[cha]->GetBinContent(bin))*(h_bckgSys[chb]->GetBinContent(bin))
	  + (h_ZZxs[cha]->GetBinContent(bin))*(h_ZZxs[chb]->GetBinContent(bin))
	  + (h_Zgammaxs[cha]->GetBinContent(bin))*(h_Zgammaxs[chb]->GetBinContent(bin))
	  + (h_unfSyst[cha]->GetBinContent(bin))*(h_unfSyst[chb]->GetBinContent(bin))
	  + (h_lumi[cha]->GetBinContent(bin))*(h_lumi[chb]->GetBinContent(bin));
	lumiSys[cha][chb]=(h_lumi[cha]->GetBinContent(bin))*(h_lumi[chb]->GetBinContent(bin))*(h_crossSection_final[cha]->GetBinContent(bin))*
	  (h_crossSection_final[chb]->GetBinContent(bin));
	if (variable!="Zpt"){
 	  commonSys[cha][chb]+=(h_JESsys[cha]->GetBinContent(bin))*(h_JESsys[chb]->GetBinContent(bin));
	}
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
    h_crossSection_final[0]->SetBinError(bin, sqrt(wo_lumi_1));
    h_crossSection_final[1]->SetBinError(bin, sqrt(wo_lumi_2));
    h_crossSection_final[2]->SetBinError(bin, sqrt(wo_lumi_3));
    h_crossSection_final[3]->SetBinError(bin, sqrt(wo_lumi_4));

    double corr_matrix_dd[4][4]={
      {1, 0, 1, 0},
      {0, 1, 0, 1}, 
      {1, 0, 1, 0},
      {0, 1, 0, 1}
    };
    //matrix is symetric
    //channels 0 and 1
    elements[4]= elements[1]= (h_crossSection_final[0]->GetBinContent(bin))*(h_crossSection_final[1]->GetBinContent(bin))*
      (commonSys[0][1] + (h_elEnScale[0]->GetBinContent(bin))*(h_elEnScale[1]->GetBinContent(bin))+
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[1]->GetBinContent(bin))+
       (h_leptTrgEff_el[0]->GetBinContent(bin))*(h_leptTrgEff_el[1]->GetBinContent(bin))+
       (h_leptTrgEff_mu[0]->GetBinContent(bin))*(h_leptTrgEff_mu[1]->GetBinContent(bin))
       //+(h_dataDrivensys[0]->GetBinContent(bin))*(h_dataDrivensys[1]->GetBinContent(bin))*corr_matrix_dd[0][1]);
       +(h_dataDrivensys_el[0]->GetBinContent(bin))*(h_dataDrivensys_el[1]->GetBinContent(bin))
       +(h_dataDrivensys_mu[0]->GetBinContent(bin))*(h_dataDrivensys_mu[1]->GetBinContent(bin)));
    
    
    //channels 0 and 2
    elements[8]=elements[2]=(h_crossSection_final[0]->GetBinContent(bin))*(h_crossSection_final[2]->GetBinContent(bin))*
      (commonSys[0][2] + (h_elEnScale[0]->GetBinContent(bin))*(h_elEnScale[2]->GetBinContent(bin)) +
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[2]->GetBinContent(bin))+
       (h_leptTrgEff_el[0]->GetBinContent(bin))*(h_leptTrgEff_el[2]->GetBinContent(2))+
       (h_leptTrgEff_mu[0]->GetBinContent(bin))*(h_leptTrgEff_mu[2]->GetBinContent(2))
       //+(h_dataDrivensys[0]->GetBinContent(bin))*(h_dataDrivensys[2]->GetBinContent(bin))*corr_matrix_dd[0][2]);
       +(h_dataDrivensys_el[0]->GetBinContent(bin))*(h_dataDrivensys_el[2]->GetBinContent(bin))
       +(h_dataDrivensys_mu[0]->GetBinContent(bin))*(h_dataDrivensys_mu[2]->GetBinContent(bin)));
    
    //channels 0 and 3
    elements[12]=elements[3]=(h_crossSection_final[0]->GetBinContent(bin))*(h_crossSection_final[3]->GetBinContent(bin))*
      (commonSys[0][3]+ (h_elEnScale[0]->GetBinContent(bin))*(h_elEnScale[3]->GetBinContent(bin))+
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[3]->GetBinContent(bin))+
       (h_leptTrgEff_el[0]->GetBinContent(bin))*(h_leptTrgEff_el[3]->GetBinContent(2))+
       (h_leptTrgEff_mu[0]->GetBinContent(bin))*(h_leptTrgEff_mu[3]->GetBinContent(2))
       //+(h_dataDrivensys[0]->GetBinContent(bin))*(h_dataDrivensys[3]->GetBinContent(bin))*corr_matrix_dd[0][3]);
       +(h_dataDrivensys_el[0]->GetBinContent(bin))*(h_dataDrivensys_el[3]->GetBinContent(bin))
       +(h_dataDrivensys_mu[0]->GetBinContent(bin))*(h_dataDrivensys_mu[3]->GetBinContent(bin)));
    
    //channels 1 and 2
    elements[9]=elements[6]= (h_crossSection_final[1]->GetBinContent(bin))*(h_crossSection_final[2]->GetBinContent(bin))*
      (commonSys[1][2]+ (h_elEnScale[1]->GetBinContent(bin))*(h_elEnScale[2]->GetBinContent(bin))
       + (h_muMomScale[1]->GetBinContent(bin))*(h_muMomScale[2]->GetBinContent(bin)) 
       + (h_leptTrgEff_el[1]->GetBinContent(bin))*(h_leptTrgEff_el[2]->GetBinContent(bin))
       + (h_leptTrgEff_mu[1]->GetBinContent(bin))*(h_leptTrgEff_mu[2]->GetBinContent(bin))
       //+ (h_dataDrivensys[1]->GetBinContent(bin))*(h_dataDrivensys[2]->GetBinContent(bin))*corr_matrix_dd[1][2]);
       + (h_dataDrivensys_el[1]->GetBinContent(bin))*(h_dataDrivensys_el[2]->GetBinContent(bin))
       + (h_dataDrivensys_mu[1]->GetBinContent(bin))*(h_dataDrivensys_mu[2]->GetBinContent(bin)));
    
    //channels 1 and 3
    elements[13]=elements[7]= (h_crossSection_final[1]->GetBinContent(bin))*(h_crossSection_final[3]->GetBinContent(bin))*
      (commonSys[1][3]+ (h_elEnScale[1]->GetBinContent(bin))*(h_elEnScale[3]->GetBinContent(bin)) 
       + (h_muMomScale[1]->GetBinContent(bin))*(h_muMomScale[3]->GetBinContent(bin))
       + (h_leptTrgEff_el[1]->GetBinContent(bin))*(h_leptTrgEff_el[3]->GetBinContent(bin))
       + (h_leptTrgEff_mu[1]->GetBinContent(bin))*(h_leptTrgEff_mu[3]->GetBinContent(bin))
       //       + (h_dataDrivensys[1]->GetBinContent(bin))*(h_dataDrivensys[3]->GetBinContent(bin))*corr_matrix_dd[1][3]);
       + (h_dataDrivensys_el[1]->GetBinContent(bin))*(h_dataDrivensys_el[3]->GetBinContent(bin))
       + (h_dataDrivensys_mu[1]->GetBinContent(bin))*(h_dataDrivensys_mu[3]->GetBinContent(bin)));
    
    //channels 2 and 3
    elements[14]=elements[11]= (h_crossSection_final[2]->GetBinContent(bin))*(h_crossSection_final[3]->GetBinContent(bin))*
      (commonSys[2][3]+ (h_elEnScale[2]->GetBinContent(bin))*(h_elEnScale[3]->GetBinContent(bin)) 
       + (h_muMomScale[2]->GetBinContent(bin))*(h_muMomScale[3]->GetBinContent(bin))
       + (h_leptTrgEff_el[2]->GetBinContent(bin))*(h_leptTrgEff_el[3]->GetBinContent(bin))
       + (h_leptTrgEff_mu[2]->GetBinContent(bin))*(h_leptTrgEff_mu[3]->GetBinContent(bin))
       //       + (h_dataDrivensys[2]->GetBinContent(bin))*(h_dataDrivensys[3]->GetBinContent(bin))*corr_matrix_dd[2][3]);
       + (h_dataDrivensys_el[2]->GetBinContent(bin))*(h_dataDrivensys_el[3]->GetBinContent(bin))
       + (h_dataDrivensys_mu[2]->GetBinContent(bin))*(h_dataDrivensys_mu[3]->GetBinContent(bin)));
    
    if (printBLUEmatrix){
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
    
    //invert !
    errMatInv = errMat.Invert();
    Double_t *mRefTest = errMat.GetMatrixArray();
    if (printBLUEmatrix){
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
    double final_Xsec = alphaCH[0]*(h_crossSection_final[0]->GetBinContent(bin)) + alphaCH[1]*(h_crossSection_final[1]->GetBinContent(bin)) 
      + alphaCH[2]*(h_crossSection_final[2]->GetBinContent(bin)) + alphaCH[3]*(h_crossSection_final[3]->GetBinContent(bin));
    
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

    Double_t stat_err_tot2 =  pow(alphaCH[0]*statisticError[0],2) + pow(alphaCH[1]*statisticError[1],2)
      +pow(alphaCH[2]*statisticError[2],2) + pow(alphaCH[3]*statisticError[3],2);
    
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
  }

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
    for (int i=1; i<=h_crossSection_final[channels]->GetNbinsX(); i++) {
      double value2 = h_crossSection_final[channels]->GetBinContent(i);
      double error2 = h_crossSection_final[channels]->GetBinError(i);
      double errorStat2 = h_totalStat[channels]->GetBinContent(i);
      double errorSyst2 = h_totalSyst[channels]->GetBinContent(i);
      double width2 = h_crossSection_final[channels]->GetBinWidth(i);
      double dsdx2 = value2/width2;
      double dsdx_err2 = dsdx2*error2/value2;
      double dsdx_errorStat= dsdx2*errorStat2/value2;
      double dsdx_errorSyst= dsdx2*errorSyst2/value2;
      h_crossSection_diff[channels]->SetBinContent(i,dsdx2);
      h_crossSection_diff[channels]->SetBinError(i,dsdx_err2);
      h_totalSyst_diff[channels]->SetBinContent(i, dsdx_errorSyst);
      h_totalStat_diff[channels]->SetBinContent(i, dsdx_errorStat);
    }
    
  }



  //LATEX OUTPUT
  if (latexOutput){
    // TString ranges[9]={"","","","","","","","","",""};
    TString rangesZpt[9]={"0-20 GeV", "20-40 GeV", "40-60 GeV", "60-80 GeV", "80-100 GeV", "100-120 GeV", "120-140 GeV", "140-200 GeV", "200-300 GeV"};
    TString rangesLeadingJetPt[9]={"30-60 GeV", "60-100 GeV", "100-150 GeV", "150-250 GeV"};
    TString rangesNjets[9]={"0 jets", "1 jet", "2 jets", "3 jets", "4 jets"};
   
   
    std::cout<<"--------------------------------------------------------"<<std::endl;
    std::cout<<"Latex output: "<<std::endl;
    std::cout<<"bin & 3e & 2e1mu & 1e2mu & 3mu & combination \\\\"<<std::endl;
    std::cout<<"\\hline"<<std::endl;

    
  
  for (int output=1; output<=h_crossSection_combination->GetNbinsX(); output++){
    outMain<<output<<" "<<h_crossSection_diff[0]->GetBinContent(output)<<" "<<h_crossSection_diff[1]->GetBinContent(output)<<" "<<h_crossSection_diff[2]->GetBinContent(output)<<" "<<h_crossSection_diff[3]->GetBinContent(output)<<" "<<h_crossSection_comb_diff->GetBinContent(output)<<std::endl; 
    outError1<<output<<" "<<h_totalStat_diff[0]->GetBinContent(output)<<" "<<h_totalStat_diff[1]->GetBinContent(output)<<" "<<h_totalStat_diff[2]->GetBinContent(output)<<" "<<h_totalStat_diff[3]->GetBinContent(output)<<" "<<h_combStat->GetBinContent(output)<<std::endl;
    outError2<<output<<" "<<h_totalSyst_diff[0]->GetBinContent(output)<<" "<<h_totalSyst_diff[1]->GetBinContent(output)<<" "<<h_totalSyst_diff[2]->GetBinContent(output)<<" "<<h_totalSyst_diff[3]->GetBinContent(output)<<" "<<h_combSyst->GetBinContent(output)<<std::endl; 
    outError3<<output<<" "<<0.026*h_crossSection_diff[0]->GetBinContent(output)<<" "<<0.026*h_crossSection_diff[1]->GetBinContent(output)<<" "<<0.026*h_crossSection_diff[2]->GetBinContent(output)<<" "<<0.026*h_crossSection_diff[3]->GetBinContent(output)<<" "<<h_combLumi->GetBinContent(output)<<std::endl; 

    if (variable=="Zpt"){
      std::cout<<scientific<<setprecision(4)<<rangesZpt[output-1]<<" & "<<h_crossSection_diff[0]->GetBinContent(output)<<"$ \\pm$ "<<h_totalStat_diff[0]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[0]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[0]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[1]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[1]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[1]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[1]->GetBinContent(output) <<" & "<<
	h_crossSection_diff[2]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[2]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[2]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[2]->GetBinContent(output) <<" & "<<
	h_crossSection_diff[3]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[3]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[3]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[3]->GetBinContent(output)<<" & "<<
      	h_crossSection_comb_diff->GetBinContent(output)<<" $\\pm$ "<<h_combStat->GetBinContent(output)<<" $\\pm$ "<<h_combSyst->GetBinContent(output)<<"$\\pm$ "<<h_combLumi->GetBinContent(output)<<"\\\\"<<std::endl;
    }
    if (variable=="LeadingJetPt"){
      std::cout<<fixed<<setprecision(4)<<rangesLeadingJetPt[output-1]<<" & "<<h_crossSection_diff[0]->GetBinContent(output)<<"$ \\pm$ "<<h_totalStat_diff[0]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[0]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[0]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[1]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[1]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[1]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[2]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[2]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[2]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[2]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[2]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[3]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[3]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[3]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[3]->GetBinContent(output)<<" & "<<
      	h_crossSection_comb_diff->GetBinContent(output)<<" $\\pm$ "<<h_combStat->GetBinContent(output)<<" $\\pm$ "<<h_combSyst->GetBinContent(output)<<"$\\pm\
$ "<<h_combLumi->GetBinContent(output)<<"\\\\"<<std::endl;
    }
    if (variable=="Njets"){
      std::cout<<fixed<<setprecision(3)<<rangesNjets[output-1]<<" & "<<h_crossSection_diff[0]->GetBinContent(output)<<"$ \\pm$ "<<h_totalStat_diff[0]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[0]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[0]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[1]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[1]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[1]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[1]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[2]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[2]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[2]->GetBinContent(output)<<" $\\pm$ "<<0.026*h_crossSection_diff[2]->GetBinContent(output)<<" & "<<
	h_crossSection_diff[3]->GetBinContent(output)<<" $\\pm$ "<<h_totalStat_diff[3]->GetBinContent(output)<<" $\\pm$ "<<h_totalSyst_diff[3]->GetBinContent(output)<<" & "<<
      	h_crossSection_comb_diff->GetBinContent(output)<<" $\\pm$ "<<h_combStat->GetBinContent(output)<<" $\\pm$ "<<h_combSyst->GetBinContent(output)<<"$\\pm\
$ "<<h_combLumi->GetBinContent(output)<<"\\\\"<<std::endl;
    }

  }
  }
  outMain.close();
  //END OF LATEX OUTPUT
  fout->cd();
  h_crossSection_final[0]->Write();
  h_crossSection_final[1]->Write();
  h_crossSection_final[2]->Write();
  h_crossSection_final[3]->Write();
  h_crossSection_diff[0]->Write();
  h_crossSection_diff[1]->Write();
  h_crossSection_diff[2]->Write();
  h_crossSection_diff[3]->Write();
  h_crossSection_combination->Write();
  h_crossSection_comb_diff->Write();
  h_combStat->Write();
  h_combSyst->Write();
 
}
