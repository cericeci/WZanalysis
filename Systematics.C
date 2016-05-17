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


#define  USEFLATQCD  false
#define  USEFLATPDF  false


using namespace std;

TH1D * GetTheorySystematics(std::string errorType,
			    std::string variable, 
			    std::string name, 
			    std::string title) {


  bool    useConstantUncertainty = true;
  double  qcdUncertainty        = 0.03;
  double  pdfUncertainty        = 0.014;
  double flatUncertainty;

  if (errorType == "qcd") {
    flatUncertainty = qcdUncertainty;
    useConstantUncertainty = USEFLATQCD;
  } else if (errorType == "pdf") {
    flatUncertainty = pdfUncertainty;   
    useConstantUncertainty = USEFLATPDF;
  }



  TH1D * hqcd = UnfoldingHistogramFactory::createHistogramForVar(variable,
								 name, title); 
  std::vector<double> values;

  // See table at
  // https://docs.google.com/spreadsheets/d/1NU-O7KRW2muRuzbYHBXY-UCr5idTlTBEZppDt5QCWvA/edit?usp=sharing

  if (errorType == "qcd") {

    if (variable == "Njets") {
      values.push_back(0.03);
      values.push_back(0.03);
      values.push_back(0.03);
      values.push_back(0.03);
    } else if (variable == "LeadingJetPt") {
      values.push_back(0.04);
      values.push_back(0.02);
      values.push_back(0.05);
      values.push_back(0.05);    
    } else if ( variable =="Zpt"){
      values.push_back(0.02);
      values.push_back(0.06);
      values.push_back(0.06);
      values.push_back(0.04);
      values.push_back(0.015);
      values.push_back(0.015);
      values.push_back(0.018);
      values.push_back(0.04);
      values.push_back(0.04);
    } 
  } else if (errorType == "pdf") {
    if (variable == "Njets") {
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
    } else if (variable == "LeadingJetPt") {
      values.push_back(0.03);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
    } else if ( variable =="Zpt"){
      values.push_back(0.03);
      values.push_back(0.03);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
      values.push_back(pdfUncertainty);
    } 
  }
  //
  int nbins = hqcd->GetNbinsX();

  if (nbins != values.size())  {
    std::cout << "QCD scale: nr of values does not correspond to nr of bins!!! \n";
    exit;
  }


  hqcd->SetBinContent(0, flatUncertainty);
  for (int ibin=1; ibin<=nbins; ibin++) {
    if (useConstantUncertainty) {
      hqcd->SetBinContent(ibin, flatUncertainty);
    } else {
      hqcd->SetBinContent(ibin, values[ibin-1]);
    }
  }

}



int main(int argc, char **argv)
{
  bool gotHistoBinning(false);
  char * binningFileName(0);
  char * variableName(0);
  bool gotVarName  = false;
  char c;
  
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
  string syst_type;

  if (gotVarName) {
    variable = variableName;
  } else {
    variable = "LeadingJetPt";
  }

  std::ostringstream outName;
  outName<<"sysResults/systematics_"<<variable<<".root";
  //  TFile * fout = new TFile("systematics.root", "RECREATE");
  TFile * fout = new TFile(outName.str().c_str(), "RECREATE");

  if (gotHistoBinning) {
    
    UnfoldingHistogramFactory * histoFac = UnfoldingHistogramFactory::GetInstance();
    histoFac->SetBinning(binningFileName);
    
  }

  double WZbr[4]={0.03363*0.1075,0.03363*0.1057,0.03366*0.1075,0.03366*0.1057};
  const int nChannels(4);

  TH1D * h_qcdScale[nChannels];
  TH1D * h_PDFSys[nChannels];
  TH1D * h_DataDrivenSys[nChannels];
  TH1D * h_DataDrivenSys_el[nChannels];
  TH1D * h_DataDrivenSys_mu[nChannels];
  TH1D * h_unf_syst[nChannels];
  TH1D * h_leptTrgEff_el[nChannels];
  TH1D * h_leptTrgEff_mu[nChannels];
  TH1D * h_Etsys[nChannels];
  TH1D * h_lumi[nChannels];
  TH1D * h_bckgSys[nChannels];
  TH1D * h_crossSection[nChannels];
  TH1D * h_crossSection_diff[nChannels];
  TH1D * h_crossSection_inclusive[nChannels];
  TH1D * h_crossSection_incl_diff[nChannels];
  TH1D * h_met_nominal[nChannels];
  TH1D * h_kterm[nChannels];

  TH1D* h_sys[nChannels][26];

  TMatrixD* covariance_unfolding[nChannels];

  std::vector<TString>types;
  types.push_back("pu_syst");
  types.push_back("ele_scale_syst");
  types.push_back("mu_scale_syst");
  types.push_back("ele_SF");
  types.push_back("mu_SF");
  types.push_back("JER");
  types.push_back("JES");
  types.push_back("ZZ");
  types.push_back("Zgamma");
  types.push_back("WV");
  types.push_back("WZZJets");
  types.push_back("ZZZJets");
  types.push_back("WWZJets");
  types.push_back("WWWJets");
  types.push_back("TTWJets");
  types.push_back("TTZJets");
  types.push_back("TTWWJets");
  types.push_back("TTGJets");
  types.push_back("WWGJets");
  types.push_back("met_elEn");
  types.push_back("met_jetEn");
  types.push_back("met_jetRes");
  types.push_back("met_muEn");
  types.push_back("met_tauEn");
  types.push_back("met_unEn");
  types.push_back("UnfoldingStat");
  //  types.push_back("met_nominal");

  for (int hist=0; hist<nChannels; hist++){
    std::ostringstream qcdScaleName,PDFsysName,leptTrgEffName_el, leptTrgEffName_mu,EtsysName, 
      muMomScaleName, elEnScaleName, pileupSysName, ZZxsName, ZgammaxsName, 
      dataDrivenName, bckgSysName, xsNameIncl, xsNameInclDiff, JESsysName, JERsysName, eleSFname, muSFname,
      unfStatName;
    std::ostringstream WVname, WZZname, ZZZname, WWZname, WWWname, TTWname, 
      TTZname, TTWWname, TTGname, WWGname, systName, lumiName;
    std::ostringstream metElEnName, metJetEnName, metJetResName, metMuEnName, metTauEnName, metUnEnName, metNominalName;
    std::ostringstream dataDriven_elName, dataDriven_muName, unfSystName;
    qcdScaleName<<"h_qcdScale_"<<hist;
    PDFsysName<<"h_PDFsys_"<<hist;
    leptTrgEffName_el<<"h_leptTrgEff_el_"<<hist;
    leptTrgEffName_mu<<"h_leptTrgEff_mu_"<<hist;
    EtsysName<<"h_Etsys_"<<hist;
    lumiName<<"h_lumi_"<<hist;
    muMomScaleName<<"h_muMomScale_"<<hist;
    elEnScaleName<<"h_elEnScale_"<<hist;
    pileupSysName<<"h_pileupSys_"<<hist;
    eleSFname<<"ele_SF_"<<hist;
    muSFname<<"mu_SF_"<<hist;
    unfStatName << "h_unfoldingStat_" << hist;
    ZZxsName<<"h_ZZxs_"<<hist;
    ZgammaxsName<<"h_Zgammaxs_"<<hist;
    dataDrivenName<<"h_dataDrivensys_"<<hist;
    dataDriven_elName<<"h_dataDrivenElsys_"<<hist;
    dataDriven_muName<<"h_dataDrivenMusys_"<<hist;
    bckgSysName<<"h_bckgSys_"<<hist;
    xsNameIncl<<"h_crossSection_inclusive"<<hist;
    xsNameInclDiff<<"h_crossSection_inclusive_diff"<<hist;
    JESsysName<<"h_JESsys_"<<hist;
    JERsysName<<"h_JERsys_"<<hist;
    WVname<<"h_WVJets_"<<hist;
    WZZname<<"h_WZZJets_"<<hist;
    ZZZname<<"h_ZZZJets_"<<hist;
    WWZname<<"h_WWZJets_"<<hist;
    WWWname<<"h_WWWJets_"<<hist;
    TTWname<<"h_TTWJets_"<<hist;
    TTZname<<"h_TTZJets_"<<hist;
    TTWWname<<"h_TTWWJets_"<<hist;
    TTGname<<"h_TTGJets_"<<hist;
    WWGname<<"h_WWGJets_"<<hist;
    systName<<"h_syst"<<hist;
    metElEnName<<"h_met_elEn_"<<hist;
    metJetEnName<<"h_met_jetEn_"<<hist;
    metJetResName<<"h_met_jetRes_"<<hist;
    metMuEnName<<"h_met_muEn_"<<hist; 
    metTauEnName<<"h_met_tauEn_"<<hist;
    metUnEnName<<"h_met_unEn_"<<hist;
    metNominalName<<"h_metNominal_"<<hist;
    unfSystName<<"h_unfSyst_"<<hist;

    h_qcdScale[hist] = GetTheorySystematics("qcd",variable, 
					    qcdScaleName.str().c_str(), 
					    qcdScaleName.str().c_str());
    h_PDFSys[hist]  = GetTheorySystematics("pdf",variable, 
					   PDFsysName.str().c_str(), 
					   PDFsysName.str().c_str());

    h_sys[hist][25] = UnfoldingHistogramFactory::createHistogramForVar(variable,
								       unfStatName.str().c_str(),
								       unfStatName.str().c_str());
    std::cout << "Unfolding stat histo : " << unfStatName.str().c_str() << std::endl;


    if (variable=="Njets"){
      //      h_qcdScale[hist]=UnfoldingHistogramFactory::createNjetsHistogram(qcdScaleName.str().c_str(), qcdScaleName.str().c_str());
      //      h_PDFSys[hist]= UnfoldingHistogramFactory::createNjetsHistogram(PDFsysName.str().c_str(), PDFsysName.str().c_str());
      h_leptTrgEff_el[hist]=UnfoldingHistogramFactory::createNjetsHistogram(leptTrgEffName_el.str().c_str(), leptTrgEffName_el.str().c_str());
      h_leptTrgEff_mu[hist]=UnfoldingHistogramFactory::createNjetsHistogram(leptTrgEffName_mu.str().c_str(), leptTrgEffName_mu.str().c_str());
      h_Etsys[hist]=UnfoldingHistogramFactory::createNjetsHistogram(EtsysName.str().c_str(), EtsysName.str().c_str());
      h_lumi[hist]=UnfoldingHistogramFactory::createNjetsHistogram(lumiName.str().c_str(), lumiName.str().c_str());
      h_DataDrivenSys[hist]=UnfoldingHistogramFactory::createNjetsHistogram(dataDrivenName.str().c_str(), dataDrivenName.str().c_str());
      h_DataDrivenSys_el[hist]=UnfoldingHistogramFactory::createNjetsHistogram(dataDriven_elName.str().c_str(), dataDriven_elName.str().c_str());
      h_DataDrivenSys_mu[hist]=UnfoldingHistogramFactory::createNjetsHistogram(dataDriven_muName.str().c_str(), dataDriven_muName.str().c_str());
      h_unf_syst[hist]=UnfoldingHistogramFactory::createNjetsHistogram(unfSystName.str().c_str(), unfSystName.str().c_str());
      h_bckgSys[hist]=UnfoldingHistogramFactory::createNjetsHistogram(bckgSysName.str().c_str(), bckgSysName.str().c_str());
      h_crossSection_inclusive[hist]=UnfoldingHistogramFactory::createNjetsHistogram(xsNameIncl.str().c_str(), xsNameIncl.str().c_str());     
      h_crossSection_incl_diff[hist]=UnfoldingHistogramFactory::createNjetsHistogram(xsNameInclDiff.str().c_str(), xsNameInclDiff.str().c_str());     
      h_met_nominal[hist]=UnfoldingHistogramFactory::createNjetsHistogram(metNominalName.str().c_str(), metNominalName.str().c_str());

      h_sys[hist][0] =UnfoldingHistogramFactory::createNjetsHistogram(pileupSysName.str().c_str(), pileupSysName.str().c_str());
      h_sys[hist][1] =UnfoldingHistogramFactory::createNjetsHistogram(elEnScaleName.str().c_str(), elEnScaleName.str().c_str());
      h_sys[hist][2] =UnfoldingHistogramFactory::createNjetsHistogram(muMomScaleName.str().c_str(), muMomScaleName.str().c_str());
      h_sys[hist][3] =UnfoldingHistogramFactory::createNjetsHistogram(eleSFname.str().c_str(), eleSFname.str().c_str());
      h_sys[hist][4] =UnfoldingHistogramFactory::createNjetsHistogram(muSFname.str().c_str(), muSFname.str().c_str());
      h_sys[hist][5] =UnfoldingHistogramFactory::createNjetsHistogram(JERsysName.str().c_str(), JERsysName.str().c_str());
      h_sys[hist][6] =UnfoldingHistogramFactory::createNjetsHistogram(JESsysName.str().c_str(), JESsysName.str().c_str());
      h_sys[hist][7] =UnfoldingHistogramFactory::createNjetsHistogram(ZZxsName.str().c_str(), ZZxsName.str().c_str());
      h_sys[hist][8] =UnfoldingHistogramFactory::createNjetsHistogram(ZgammaxsName.str().c_str(), ZgammaxsName.str().c_str());
      h_sys[hist][9] =UnfoldingHistogramFactory::createNjetsHistogram(WVname.str().c_str(), WVname.str().c_str());
      h_sys[hist][10]=UnfoldingHistogramFactory::createNjetsHistogram(WZZname.str().c_str(), WZZname.str().c_str());
      h_sys[hist][11]=UnfoldingHistogramFactory::createNjetsHistogram(ZZZname.str().c_str(), ZZZname.str().c_str());
      h_sys[hist][12]=UnfoldingHistogramFactory::createNjetsHistogram(WWZname.str().c_str(), WWZname.str().c_str());
      h_sys[hist][13]=UnfoldingHistogramFactory::createNjetsHistogram(WWWname.str().c_str(), WWWname.str().c_str());
      h_sys[hist][14]=UnfoldingHistogramFactory::createNjetsHistogram(TTWname.str().c_str(), TTWname.str().c_str());
      h_sys[hist][15]=UnfoldingHistogramFactory::createNjetsHistogram(TTZname.str().c_str(), TTZname.str().c_str());
      h_sys[hist][16]=UnfoldingHistogramFactory::createNjetsHistogram(TTWWname.str().c_str(), TTWWname.str().c_str());
      h_sys[hist][17]=UnfoldingHistogramFactory::createNjetsHistogram(TTGname.str().c_str(), TTGname.str().c_str());
      h_sys[hist][18]=UnfoldingHistogramFactory::createNjetsHistogram(WWGname.str().c_str(), WWGname.str().c_str());
      h_sys[hist][19]=UnfoldingHistogramFactory::createNjetsHistogram(metElEnName.str().c_str(), metElEnName.str().c_str());
      h_sys[hist][20]=UnfoldingHistogramFactory::createNjetsHistogram(metJetEnName.str().c_str(), metJetEnName.str().c_str());
      h_sys[hist][21]=UnfoldingHistogramFactory::createNjetsHistogram(metJetResName.str().c_str(), metJetResName.str().c_str());
      h_sys[hist][22]=UnfoldingHistogramFactory::createNjetsHistogram(metMuEnName.str().c_str(), metMuEnName.str().c_str());
      h_sys[hist][23]=UnfoldingHistogramFactory::createNjetsHistogram(metTauEnName.str().c_str(), metTauEnName.str().c_str());
      h_sys[hist][24]=UnfoldingHistogramFactory::createNjetsHistogram(metUnEnName.str().c_str(), metUnEnName.str().c_str());
    }

    if (variable=="Zpt"){
      //      h_qcdScale[hist]=UnfoldingHistogramFactory::createZPtHistogram(qcdScaleName.str().c_str(), qcdScaleName.str().c_str());
      //      h_PDFSys[hist]= UnfoldingHistogramFactory::createZPtHistogram(PDFsysName.str().c_str(), PDFsysName.str().c_str());
      h_leptTrgEff_el[hist]=UnfoldingHistogramFactory::createZPtHistogram(leptTrgEffName_el.str().c_str(), leptTrgEffName_el.str().c_str());
      h_leptTrgEff_mu[hist]=UnfoldingHistogramFactory::createZPtHistogram(leptTrgEffName_mu.str().c_str(), leptTrgEffName_mu.str().c_str());
      h_Etsys[hist]=UnfoldingHistogramFactory::createZPtHistogram(EtsysName.str().c_str(), EtsysName.str().c_str());
      h_lumi[hist]=UnfoldingHistogramFactory::createZPtHistogram(lumiName.str().c_str(), lumiName.str().c_str());
      h_DataDrivenSys[hist]=UnfoldingHistogramFactory::createZPtHistogram(dataDrivenName.str().c_str(), dataDrivenName.str().c_str());
      h_DataDrivenSys_el[hist]=UnfoldingHistogramFactory::createZPtHistogram(dataDriven_elName.str().c_str(), dataDriven_elName.str().c_str());
      h_DataDrivenSys_mu[hist]=UnfoldingHistogramFactory::createZPtHistogram(dataDriven_muName.str().c_str(), dataDriven_muName.str().c_str());
      h_unf_syst[hist]=UnfoldingHistogramFactory::createZPtHistogram(unfSystName.str().c_str(), unfSystName.str().c_str());
      h_bckgSys[hist]=UnfoldingHistogramFactory::createZPtHistogram(bckgSysName.str().c_str(), bckgSysName.str().c_str());
      h_crossSection_inclusive[hist]=UnfoldingHistogramFactory::createZPtHistogram(xsNameIncl.str().c_str(), xsNameIncl.str().c_str());     
      h_crossSection_incl_diff[hist]=UnfoldingHistogramFactory::createZPtHistogram(xsNameInclDiff.str().c_str(), xsNameInclDiff.str().c_str());     
      h_met_nominal[hist]=UnfoldingHistogramFactory::createZPtHistogram(metNominalName.str().c_str(), metNominalName.str().c_str());

      h_sys[hist][0] =UnfoldingHistogramFactory::createZPtHistogram(pileupSysName.str().c_str(), pileupSysName.str().c_str());
      h_sys[hist][1] =UnfoldingHistogramFactory::createZPtHistogram(elEnScaleName.str().c_str(), elEnScaleName.str().c_str());
      h_sys[hist][2] =UnfoldingHistogramFactory::createZPtHistogram(muMomScaleName.str().c_str(), muMomScaleName.str().c_str());
      h_sys[hist][3] =UnfoldingHistogramFactory::createZPtHistogram(eleSFname.str().c_str(), eleSFname.str().c_str());
      h_sys[hist][4] =UnfoldingHistogramFactory::createZPtHistogram(muSFname.str().c_str(), muSFname.str().c_str());
      h_sys[hist][5] =UnfoldingHistogramFactory::createZPtHistogram(JERsysName.str().c_str(), JERsysName.str().c_str());
      h_sys[hist][6] =UnfoldingHistogramFactory::createZPtHistogram(JESsysName.str().c_str(), JESsysName.str().c_str());
      h_sys[hist][7] =UnfoldingHistogramFactory::createZPtHistogram(ZZxsName.str().c_str(), ZZxsName.str().c_str());
      h_sys[hist][8] =UnfoldingHistogramFactory::createZPtHistogram(ZgammaxsName.str().c_str(), ZgammaxsName.str().c_str());
      h_sys[hist][9] =UnfoldingHistogramFactory::createZPtHistogram(WVname.str().c_str(), WVname.str().c_str());
      h_sys[hist][10]=UnfoldingHistogramFactory::createZPtHistogram(WZZname.str().c_str(), WZZname.str().c_str());
      h_sys[hist][11]=UnfoldingHistogramFactory::createZPtHistogram(ZZZname.str().c_str(), ZZZname.str().c_str());
      h_sys[hist][12]=UnfoldingHistogramFactory::createZPtHistogram(WWZname.str().c_str(), WWZname.str().c_str());
      h_sys[hist][13]=UnfoldingHistogramFactory::createZPtHistogram(WWWname.str().c_str(), WWWname.str().c_str());
      h_sys[hist][14]=UnfoldingHistogramFactory::createZPtHistogram(TTWname.str().c_str(), TTWname.str().c_str());
      h_sys[hist][15]=UnfoldingHistogramFactory::createZPtHistogram(TTZname.str().c_str(), TTZname.str().c_str());
      h_sys[hist][16]=UnfoldingHistogramFactory::createZPtHistogram(TTWWname.str().c_str(), TTWWname.str().c_str());
      h_sys[hist][17]=UnfoldingHistogramFactory::createZPtHistogram(TTGname.str().c_str(), TTGname.str().c_str());
      h_sys[hist][18]=UnfoldingHistogramFactory::createZPtHistogram(WWGname.str().c_str(), WWGname.str().c_str());
      h_sys[hist][19]=UnfoldingHistogramFactory::createZPtHistogram(metElEnName.str().c_str(), metElEnName.str().c_str());
      h_sys[hist][20]=UnfoldingHistogramFactory::createZPtHistogram(metJetEnName.str().c_str(), metJetEnName.str().c_str());
      h_sys[hist][21]=UnfoldingHistogramFactory::createZPtHistogram(metJetResName.str().c_str(), metJetResName.str().c_str());
      h_sys[hist][22]=UnfoldingHistogramFactory::createZPtHistogram(metMuEnName.str().c_str(), metMuEnName.str().c_str());
      h_sys[hist][23]=UnfoldingHistogramFactory::createZPtHistogram(metTauEnName.str().c_str(), metTauEnName.str().c_str());
      h_sys[hist][24]=UnfoldingHistogramFactory::createZPtHistogram(metUnEnName.str().c_str(), metUnEnName.str().c_str());
    }


    if (variable=="LeadingJetPt"){
      //      h_qcdScale[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(qcdScaleName.str().c_str(), qcdScaleName.str().c_str());
      //      h_PDFSys[hist]= UnfoldingHistogramFactory::createLeadingJetHistogram(PDFsysName.str().c_str(), PDFsysName.str().c_str());
      h_leptTrgEff_el[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(leptTrgEffName_el.str().c_str(), leptTrgEffName_el.str().c_str());
      h_leptTrgEff_mu[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(leptTrgEffName_mu.str().c_str(), leptTrgEffName_mu.str().c_str());
      h_Etsys[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(EtsysName.str().c_str(), EtsysName.str().c_str());
      h_lumi[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(lumiName.str().c_str(), lumiName.str().c_str());
      h_DataDrivenSys[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(dataDrivenName.str().c_str(), dataDrivenName.str().c_str());
      h_DataDrivenSys_el[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(dataDriven_elName.str().c_str(), dataDriven_elName.str().c_str());
      h_DataDrivenSys_mu[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(dataDriven_muName.str().c_str(), dataDriven_muName.str().c_str());
      h_unf_syst[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(unfSystName.str().c_str(), unfSystName.str().c_str());
      h_bckgSys[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(bckgSysName.str().c_str(), bckgSysName.str().c_str());
      h_crossSection_inclusive[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(xsNameIncl.str().c_str(), xsNameIncl.str().c_str());     
      h_crossSection_incl_diff[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(xsNameInclDiff.str().c_str(), xsNameInclDiff.str().c_str());     
      h_met_nominal[hist]=UnfoldingHistogramFactory::createLeadingJetHistogram(metNominalName.str().c_str(), metNominalName.str().c_str());

      h_sys[hist][0] =UnfoldingHistogramFactory::createLeadingJetHistogram(pileupSysName.str().c_str(), pileupSysName.str().c_str());
      h_sys[hist][1] =UnfoldingHistogramFactory::createLeadingJetHistogram(elEnScaleName.str().c_str(), elEnScaleName.str().c_str());
      h_sys[hist][2] =UnfoldingHistogramFactory::createLeadingJetHistogram(muMomScaleName.str().c_str(), muMomScaleName.str().c_str());
      h_sys[hist][3] =UnfoldingHistogramFactory::createLeadingJetHistogram(eleSFname.str().c_str(), eleSFname.str().c_str());
      h_sys[hist][4] =UnfoldingHistogramFactory::createLeadingJetHistogram(muSFname.str().c_str(), muSFname.str().c_str());
      h_sys[hist][5] =UnfoldingHistogramFactory::createLeadingJetHistogram(JERsysName.str().c_str(), JERsysName.str().c_str());
      h_sys[hist][6] =UnfoldingHistogramFactory::createLeadingJetHistogram(JESsysName.str().c_str(), JESsysName.str().c_str());
      h_sys[hist][7] =UnfoldingHistogramFactory::createLeadingJetHistogram(ZZxsName.str().c_str(), ZZxsName.str().c_str());
      h_sys[hist][8] =UnfoldingHistogramFactory::createLeadingJetHistogram(ZgammaxsName.str().c_str(), ZgammaxsName.str().c_str());
      h_sys[hist][9] =UnfoldingHistogramFactory::createLeadingJetHistogram(WVname.str().c_str(), WVname.str().c_str());
      h_sys[hist][10]=UnfoldingHistogramFactory::createLeadingJetHistogram(WZZname.str().c_str(), WZZname.str().c_str());
      h_sys[hist][11]=UnfoldingHistogramFactory::createLeadingJetHistogram(ZZZname.str().c_str(), ZZZname.str().c_str());
      h_sys[hist][12]=UnfoldingHistogramFactory::createLeadingJetHistogram(WWZname.str().c_str(), WWZname.str().c_str());
      h_sys[hist][13]=UnfoldingHistogramFactory::createLeadingJetHistogram(WWWname.str().c_str(), WWWname.str().c_str());
      h_sys[hist][14]=UnfoldingHistogramFactory::createLeadingJetHistogram(TTWname.str().c_str(), TTWname.str().c_str());
      h_sys[hist][15]=UnfoldingHistogramFactory::createLeadingJetHistogram(TTZname.str().c_str(), TTZname.str().c_str());
      h_sys[hist][16]=UnfoldingHistogramFactory::createLeadingJetHistogram(TTWWname.str().c_str(), TTWWname.str().c_str());
      h_sys[hist][17]=UnfoldingHistogramFactory::createLeadingJetHistogram(TTGname.str().c_str(), TTGname.str().c_str());
      h_sys[hist][18]=UnfoldingHistogramFactory::createLeadingJetHistogram(WWGname.str().c_str(), WWGname.str().c_str());
      h_sys[hist][19]=UnfoldingHistogramFactory::createLeadingJetHistogram(metElEnName.str().c_str(), metElEnName.str().c_str());
      h_sys[hist][20]=UnfoldingHistogramFactory::createLeadingJetHistogram(metJetEnName.str().c_str(), metJetEnName.str().c_str());
      h_sys[hist][21]=UnfoldingHistogramFactory::createLeadingJetHistogram(metJetResName.str().c_str(), metJetResName.str().c_str());
      h_sys[hist][22]=UnfoldingHistogramFactory::createLeadingJetHistogram(metMuEnName.str().c_str(), metMuEnName.str().c_str());
      h_sys[hist][23]=UnfoldingHistogramFactory::createLeadingJetHistogram(metTauEnName.str().c_str(), metTauEnName.str().c_str());
      h_sys[hist][24]=UnfoldingHistogramFactory::createLeadingJetHistogram(metUnEnName.str().c_str(), metUnEnName.str().c_str());
    }
    
   
  }

  std::ostringstream nominalName, nominalNameMet;
  nominalName<<"sysResults/unfolding_nominal_"<<variable<<".root";
  nominalNameMet<<"sysResults/unfolding_met_nominal_"<<variable<<".root";

  std::cout<<nominalName.str().c_str()<<std::endl;

  for (int i=0; i<types.size(); i++){

    std::ostringstream fileNameUp, fileNameDown;
    fileNameUp<<"sysResults/unfolding_"<<types[i]<<"_"<<variable<<"_UP.root";
    fileNameDown<<"sysResults/unfolding_"<<types[i]<<"_"<<variable<<"_DOWN.root";

    TFile * fUP= TFile::Open(fileNameUp.str().c_str());
    TFile * fDOWN= TFile::Open(fileNameDown.str().c_str());

    TFile* fnominal;
    std::cout<<"Nominal name MET: "<<nominalNameMet.str().c_str()<<std::endl;;
    if (i>=19 && i<25){
      fnominal= TFile::Open(nominalNameMet.str().c_str());
    }
    else if (i==25) {
      fnominal = fUP;
    } else 
      fnominal= TFile::Open(nominalName.str().c_str());

    
    for (int compute=0; compute<nChannels; compute++){
      std::ostringstream histName, histNameNewUp, histNameNewDown, fileNameNominal;
      histName<<"hdsigma"<<variable<<"_"<<(compute+1);
      histNameNewUp<<"hdsigma"<<variable<<"_UP_"<<compute;
      histNameNewDown<<"hdsigma"<<variable<<"_DOWN_"<<compute;
      TH1D * h_up= (TH1D*) (fUP->Get(histName.str().c_str())->Clone(histNameNewUp.str().c_str()));
      TH1D * h_down= (TH1D*) (fDOWN->Get(histName.str().c_str())->Clone(histNameNewDown.str().c_str()));
      TH1D * h_nominal= (TH1D*) (fnominal->Get(histName.str().c_str())->Clone(histName.str().c_str()));    
      for (int bin=0; bin< (h_nominal->GetNbinsX()+1); bin++){
	double up_value=h_up->GetBinContent(bin);
	double down_value=h_down->GetBinContent(bin);
	double nominal_value=h_nominal->GetBinContent(bin);
	double maxDiff=std::max(fabs(up_value-nominal_value), fabs(down_value-nominal_value));
	double syst_value=0;
	//double syst_value= maxDiff/nominal_value;
	
	if (nominal_value!=0.0){
	  syst_value= maxDiff/nominal_value;
	}
	if ((i==24) && (compute==2)){
	  std::cout<<types[i]<<std::endl;
	  std::cout<<"up: "<<up_value<<" down: "<<down_value<<" nominal: "<<nominal_value<<std::endl;
	  std::cout<<"sys: "<<syst_value<<std::endl;
	}
	h_sys[compute][i]->SetBinContent(bin, syst_value);
      }
      delete h_nominal;
      delete h_up;
      delete h_down;
    }
  }

  // Unfolding Statistics: make an average of the variation for each bin and 
  // apply it to all

  TH1D * h_unfoldStatSys_all =  
    UnfoldingHistogramFactory::createHistogramForVar(variable
						     ,"h_unfoldingStat_mean"
						     ,"Mean unfolding stat. sys");


  for (int ichan=0; ichan<nChannels; ichan++) {
    h_unfoldStatSys_all->Add(h_sys[ichan][25]);    
  }
  // factors:
  // 4 for average
  // sqrt(2) for going from difference between elments to sigma
  // sqrt(2) as stat. of full sample two times larger
  h_unfoldStatSys_all->Scale(1/8.); 

  // Copy this for all 4 channels
  TH1D * unfStatHistos[4];
  for (int ichan=0; ichan<nChannels; ichan++) {
    std::ostringstream key;
    key << "h_unfStat_" << ichan;
    unfStatHistos[ichan] = (TH1D*) h_unfoldStatSys_all->Clone(key.str().c_str());
  }
  
  std::ostringstream dataDrivenName, dataDrivenName_el, dataDrivenName_mu, ktermName_up, ktermName_down;
  dataDrivenName_el<<"sysResults/unfolding_dataDriven_el_"<<variable<<".root";
  dataDrivenName_mu<<"sysResults/unfolding_dataDriven_mu_"<<variable<<".root";
  ktermName_up<<"sysResults/unfolding_kterm_up_"<<variable<<".root";
  ktermName_down<<"sysResults/unfolding_kterm_down_"<<variable<<".root";

  TFile * fDataDriven_el=TFile::Open(dataDrivenName_el.str().c_str());
  TFile * fDataDriven_mu=TFile::Open(dataDrivenName_mu.str().c_str());
  TFile * fkTerm_up= TFile::Open(ktermName_up.str().c_str());
  TFile * fkTerm_down=TFile::Open(ktermName_down.str().c_str());
  
  
  for (int other=0; other<nChannels; other++){
    TFile * fnominal= TFile::Open(nominalName.str().c_str());  
    //cross section
    std::ostringstream crossSect, crossSect2, dataDrivSys, xsName, xsName2;
    crossSect<<"hdsigma"<<variable<<"_"<<(other+1);
    xsName<<"h_crossSection_"<<other;
    h_crossSection[other]= (TH1D*) (fnominal)->Get(crossSect.str().c_str())->Clone(xsName.str().c_str());

    // Adding Covariance matrix reading
    std::ostringstream covarianceInName;
    std::ostringstream covarianceOutName;
    covarianceInName << "Covariance_hresult" << variable << "_" << other+1;
    covariance_unfolding[other] = (TMatrixD*) (fnominal)->Get(covarianceInName.str().c_str())->Clone(covarianceOutName.str().c_str());
    //    h_crossSection[other]= (TH1D*) (fnominal)->Get(crossSect.str().c_str())->Clone(crossSect.str().c_str());
    std::cout<<"Integral: "<<h_crossSection[other]->Integral()<<std::endl;
    TH1D * h_dataDrivenUp_el=(TH1D*)(fDataDriven_el)->Get(crossSect.str().c_str())->Clone();
    TH1D * h_dataDrivenUp_mu=(TH1D*)(fDataDriven_mu)->Get(crossSect.str().c_str())->Clone();
    TH1D * h_kterm_down=(TH1D*)(fkTerm_up)->Get(crossSect.str().c_str())->Clone();
    TH1D * h_kterm_up=  (TH1D*)(fkTerm_down)->Get(crossSect.str().c_str())->Clone();
    for (int bin1=0; bin1<(h_crossSection[other]->GetNbinsX()+1);bin1++){
      double xs=(h_crossSection[other]->GetBinContent(bin1))/WZbr[other];
      h_crossSection_inclusive[other]->SetBinContent(bin1, xs);
      double back2=0;
      for (int b=9; b<19; b++){
	back2+=pow(h_sys[other][b]->GetBinContent(bin1),2);
      }
      h_bckgSys[other]->SetBinContent(bin1, sqrt(back2));
      double met2=0;

      for (int m=19; m<=24; m++){
	met2+=pow(h_sys[other][m]->GetBinContent(bin1),2);
      }
      h_Etsys[other]->SetBinContent(bin1,sqrt(met2));
      double ltrig=sqrt(pow(h_sys[other][3]->GetBinContent(bin1),2)+pow(h_sys[other][4]->GetBinContent(bin1),2));
      //      h_leptTrgEff_el[other]->SetBinContent(bin1, ltrig);
      // ORIGINAL VALUE      
      //      h_qcdScale[other]->SetBinContent(bin1, 0.03);
      // NEW 
      //      h_qcdScale[other]->SetBinContent(bin1, 0.05);
      h_lumi[other]->SetBinContent(bin1, 0.026);
      //      h_PDFSys[other]->SetBinContent(bin1, 0.014);
      double dataDrivenSyst_el=fabs(h_crossSection[other]->GetBinContent(bin1)-h_dataDrivenUp_el->GetBinContent(bin1))/(h_crossSection[other]->GetBinContent(bin1));
      double dataDrivenSyst_mu=fabs(h_crossSection[other]->GetBinContent(bin1)-h_dataDrivenUp_mu->GetBinContent(bin1))/(h_crossSection[other]->GetBinContent(bin1));
      h_DataDrivenSys_el[other]->SetBinContent(bin1, dataDrivenSyst_el);
      h_DataDrivenSys_mu[other]->SetBinContent(bin1, dataDrivenSyst_mu);

      double unf_syst_max=std::max(fabs(h_crossSection[other]->GetBinContent(bin1)-h_kterm_up->GetBinContent(bin1)),
				   fabs(h_crossSection[other]->GetBinContent(bin1)-h_kterm_down->GetBinContent(bin1)));
      std::cout<<h_kterm_up->GetBinContent(bin1)<<" , "<<h_kterm_down->GetBinContent(bin1)<<std::endl;
      double unf_syst=unf_syst_max/(h_crossSection[other]->GetBinContent(bin1));
      h_unf_syst[other]->SetBinContent(bin1, unf_syst);
    }

  }
  
  fout->cd();
  h_crossSection[0]->Write();
  h_crossSection[1]->Write();
  h_crossSection[2]->Write();
  h_crossSection[3]->Write();

  
  h_crossSection_inclusive[0]->Write();
  h_crossSection_inclusive[1]->Write();
  h_crossSection_inclusive[2]->Write();
  h_crossSection_inclusive[3]->Write();

  for (int ich=0; ich<nChannels; ich++) {

    std::ostringstream matrixName;
    matrixName << "unfolding_covariance_" << variable << "_ch" << ich;
    covariance_unfolding[ich]->Write(matrixName.str().c_str());
  // covariance_unfolding[0]->Write("covariance_unfolding_0");
  // covariance_unfolding[1]->Write("covariance_unfolding_1");
  // covariance_unfolding[2]->Write("covariance_unfolding_2");
  // covariance_unfolding[3]->Write("covariance_unfolding_3");

  }

  h_unfoldStatSys_all->Write();
  for (int i=0; i<nChannels; i++) 
    unfStatHistos[i]->Write();
  
  fout->Write();
  fout->Close();

}
