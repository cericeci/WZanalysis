#include "TFile.h"
#include "TString.h"
#include "/users/ltikvica/CMSSW_4_2_9_HLT1/src/WZanalysis/xsections.h"
#include <iostream>
#include <sstream>
#include "TH1F.h"
#include "TPad.h" 
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <iomanip>

//#include "CMS_lumi.h"

TH1F *  DrawOverflow(int rebin, TH1F *h)
{
  // This function paint the histogram h with an extra bin for overflows

  const char* name  = h->GetName();
  //  const char* title = h->GetTitle();
  std::ostringstream title;
  title<<h->GetTitle()<<"_test";
  Int_t nxold    = h->GetNbinsX();
  Int_t nx    = h->GetNbinsX()+ rebin ;
  std::cout << "nx "<< nx << std::endl;
  Double_t x1 = h->GetBinLowEdge(1);
  //Double_t x1 = h->GetXaxis()->GetXmin();
  std::cout << "x1 "<< x1 << std::endl;
  //std::cout << "x2old "<<  h->GetXaxis()->GetXmax()<<std::endl;
  Double_t bw = h->GetBinWidth(1);
  // std::cout << "bw "<< bw << std::endl;
  Double_t x2 = h->GetBinLowEdge(nxold + 1)+  bw  * rebin;
  //Double_t x2 = h->GetXaxis()->GetXmax() +bw * rebin;
  std::cout << "x2 "<< x2 << std::endl;
  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(name, title.str().c_str(), nx, x1, x2);

  // Fill the new hitogram including the extra bin for overflows
  for (Int_t i=1; i<= nxold + 1; i++) {
    // htmp->Fill(htmp->GetBinCenter(i) - bw  , h->GetBinContent(i));
    htmp->Fill( i , h->GetBinContent(i+1));
    //        std::cout << "bin low edge" << htmp->GetBinLowEdge(i) << " get bin content " << htmp->GetBinContent(i) << std::endl; 
  }



  // Fill the underflows, already done
  htmp->Fill(x1-1, h->GetBinContent(0));


  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());
  return htmp;

}

void MoveOverflowBins(TH1*     h,
		      Double_t xmin,
		      Double_t xmax)
{
  UInt_t nbins = h->GetNbinsX();

  TAxis* axis = (TAxis*)h->GetXaxis();
  
  UInt_t firstBin = (xmin != -999) ? axis->FindBin(xmin) : 1;
  UInt_t lastBin  = (xmax != -999) ? axis->FindBin(xmax) : nbins;

  Double_t firstVal = 0;
  Double_t firstErr = 0;

  Double_t lastVal = 0;
  Double_t lastErr = 0;

  for (UInt_t i=0; i<=nbins+1; i++) {

    if (i <= firstBin) {
      firstVal += h->GetBinContent(i);
      firstErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i >= lastBin) {
      lastVal += h->GetBinContent(i);
      lastErr += (h->GetBinError(i)*h->GetBinError(i));
    }

    if (i < firstBin || i > lastBin) {
      h->SetBinContent(i, 0);
      h->SetBinError  (i, 0);
    }
  }

  firstErr = sqrt(firstErr);
  lastErr  = sqrt(lastErr);

  h->SetBinContent(firstBin, firstVal);
  h->SetBinError  (firstBin, firstErr);

  h->SetBinContent(lastBin, lastVal);
  h->SetBinError  (lastBin, lastErr);
}


void plotLatex_B(double lumi){
  TLatex latexLabel;
  latexLabel.SetNDC();
  latexLabel.SetTextSize(0.05);
  latexLabel.DrawLatex(0.18, 0.96, "#font[132]{CMS Preliminary 2012}");
  latexLabel.DrawLatex(0.77, 0.955, "#font[132]{#sqrt{s} = 8 TeV}");
  latexLabel.SetTextSize(0.04);
  latexLabel.DrawLatex(0.71, 0.87, Form("#font[132]{#intL dt= %5.2f fb^{-1}}", lumi));
  // latexLabel.DrawLatex(0.3, 0.87, "|#eta^{#gamma}| < 1.4442");
}

void styleHisto1D(TH1F* histo, const char* titleX, double binWidth, bool unit){
  //  histo->GetXaxis()->SetTitle(titleX);
  //  std::cout<<histo->GetXaxis()<<std::endl;
  if (unit)
    histo->GetYaxis()->SetTitle(Form("Number of Events / %1.2f GeV", binWidth));
  else
    histo->GetYaxis()->SetTitle(Form("Number of Events / %1.2f", binWidth));
  histo->GetXaxis()->SetLabelFont(132);
  histo->GetYaxis()->SetLabelFont(132);
  histo->GetXaxis()->SetLabelOffset(0.007);
  histo->GetYaxis()->SetLabelOffset(0.007);
  histo->GetXaxis()->SetLabelSize(0.);
  histo->GetYaxis()->SetLabelSize(0.03);
  histo->GetXaxis()->SetTitleFont(132);
  histo->GetYaxis()->SetTitleFont(132);
  histo->GetXaxis()->SetTitleSize(0.06);
  histo->GetYaxis()->SetTitleSize(0.03);
  histo->GetXaxis()->SetTitleOffset(1.);
  //  histo->GetYaxis()->SetTitleOffset(1.);
histo->GetYaxis()->SetTitleOffset(1.5);
  histo->GetXaxis()->SetNdivisions(510);
  histo->GetYaxis()->SetNdivisions(510);

  //  std::cout<<"Ovdje..."<<std::endl;
  histo->SetMarkerSize(0.7);

  histo->SetMaximum(histo->GetMaximum()*1.5);
  //  histo->SetMaximum(histo->GetMaximum()*2.0);

  //histo->Draw("SAMEPE1");
  
  histo->Draw("PE1");
}

double skimEfficiency(TFile * f) {

  
  TH1F * numEvents = (TH1F* ) f->Get("numEvents");
  TH1F * nrEvents = (TH1F* ) f->Get("hNrEvents");
  
  
  int nevents =   numEvents->GetBinContent(1);
  int nevents_pat    =   numEvents->GetBinContent(2);
  int nevents_filter =   numEvents->GetBinContent(3);

  int  nskim= nrEvents->GetEntries();

  double skimEfficiency = (float) nskim / (float) nevents;
  return skimEfficiency;
  
}



int  plotDataVsMC(TFile * fDataDriven,
		  TFile * fZz,
		  TFile * fZgamma,
		  TFile * fWv,
		  TFile * fVvv,
		  TFile * fWz,
		  TFile * fdata,
		  TString histKey ,
		  int channel,
		  TCanvas * canv,
		  TString xAxisTitle="",
		  double binWidth=0,
		  bool unit=true)
		  
{

  //luminosity:
  double luminosity(18258);
  //cross sections:



  TH1F * h_DataDriven   = (TH1F*) (fDataDriven ->Get(histKey))       ->Clone(histKey+"_dataDriven");
  TH1F * h_Zz           = (TH1F*) (fZz ->Get(histKey))    ->Clone(histKey+"_zz");
  TH1F * h_Zgamma       = (TH1F*) (fZgamma ->Get(histKey))    ->Clone(histKey+"_zgamma");
  TH1F * h_Wv           = (TH1F*) (fWv ->Get(histKey))    ->Clone(histKey+"_wv");
  TH1F * h_Vvv          = (TH1F*) (fVvv ->Get(histKey))    ->Clone(histKey+"_vvv");

  TH1F * h_Wz           = (TH1F*) (fWz    ->Get(histKey))       ->Clone(histKey+"_wz");

  TH1F * h_data         = (TH1F*) (fdata    ->Get(histKey))     ->Clone("Data");
  
  TH1F * h_All=(TH1F*) h_DataDriven->Clone("h_All");
  h_All->Add(h_Zz);
  h_All->Add(h_Zgamma);
  h_All->Add(h_Wv);
  h_All->Add(h_Vvv);
  h_All->Add(h_Wz);

  double all_error[5]={0.054, 0.0696, 0.07972, 0.072, 0.037};

  //setting error
  for (int bins=0; bins<(h_All->GetNbinsX()+1); bins++){
    double binError=all_error[channel]*(h_All->GetBinContent(bins));
    h_All->SetBinError(bins, binError);
  }

  THStack *  th = new THStack(histKey,histKey);
  th->Add(h_Zgamma);
  th->Add(h_Wv);  
  th->Add(h_Vvv);    
  th->Add(h_Zz);
  th->Add(h_DataDriven);
  th->Add(h_Wz);   


  h_Wz         ->SetFillColor(kOrange-2);    
  h_Zgamma     ->SetFillColor(kRed+2);
  h_Zz         ->SetFillColor(kRed+1);    
  h_Wv         ->SetFillColor(kAzure);    
  h_DataDriven ->SetFillColor(kGray+1);    
  h_Vvv        ->SetFillColor(kBlack);    

  h_Wz         ->SetLineColor(kOrange-2);    
  h_Zgamma     ->SetLineColor(kRed+2);
  h_DataDriven ->SetLineColor(kGray+1);    
  h_Zz         ->SetLineColor(kRed+1);    
  h_Wv         ->SetLineColor(kAzure);    
  h_Vvv        ->SetLineColor(kBlack);    


  h_All->SetFillColor(kBlack);
  h_All->SetFillStyle(3345);
  //  h_All->SetLineColor(kWhite);
  // h_All->SetLineWidth(0);
  //h_All->SetMarkerColor(kOrange-2);
  //h_All->SetMarkerSize(0);
  
  
  float dataMax = h_data->GetMaximum();
  float mcMax =   th->GetMaximum();
  /*  
  std::cout<<"data max:"<<dataMax<<std::endl;
  std::cout<<"MC max:"<< mcMax<<std::endl;
  */
  
  if (dataMax > mcMax) {
    th->SetMaximum(dataMax);
  }
  th->SetMinimum(0);
  h_data->SetMinimum(0);
  //else {th->SetMaximum(mcMax);}
  float hmin = th->GetMinimum();
  
  if (hmin>0.) th->SetMinimum(0.5*hmin);
  
  TAxis * a =  th->GetXaxis(); // ->SetTitle(xAxisTitle);

  h_data->SetMarkerStyle(20);
  h_data->SetLabelFont(62);
  //  gPad->SetBottomMargin(0.28);
  gPad->SetBottomMargin(0.29);
  
  styleHisto1D(h_data, xAxisTitle,binWidth, unit);
  
  //th->Draw();
  //th->Draw("SAME");  
  h_All->Draw("e2");
  
  //  h_data->Draw();  
  //  h_data->Draw("SAMEPE1");  


  //plotLatex_B(19.602);


  TH1F *pull = new TH1F("Pull_"+histKey,"",h_data->GetNbinsX(),h_data->GetBinLowEdge(1),h_data->GetBinLowEdge(h_data->GetNbinsX())+h_data->GetBinWidth(1));
  pull->Sumw2();
  //pull->SetBinErrorOption(TH1::kPoisson);
  
  TH1F *pulltmp = new TH1F("Pulltmp"+histKey,"",h_data->GetNbinsX(),h_data->GetBinLowEdge(1),
			   h_data->GetBinLowEdge(h_data->GetNbinsX())+h_data->GetBinWidth(1));
  pulltmp->Sumw2();

  pulltmp->Add(h_Wz);
  pulltmp->Add(h_Vvv);
  pulltmp->Add(h_Wv);
  pulltmp->Add(h_Zgamma);
  pulltmp->Add(h_Zz);
  pulltmp->Add(h_DataDriven);



  for (int j=1;j<=(h_data->GetNbinsX());j++) {
    if ((pulltmp->GetBinContent(j) != 0.) && (h_data->GetBinContent(j)!=0)) {
  //  std::cout<<":)"<<std::endl;
      pull->SetBinContent(j,h_data->GetBinContent(j)/pulltmp->GetBinContent(j)-1.);
    
      pull->SetBinError(j,h_data->GetBinError(j)/pulltmp->GetBinContent(j));
    }
  }
  
  pull->SetMinimum(0-1.);
  pull->SetMaximum(3.5-1.);
  pull->SetMarkerSize(0.7*1.5);	
  
  TPad * pad = new TPad("pad"+histKey, "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.73);//0.72
  pad->SetBottomMargin(0.11);//0.72
  pad->SetFillColor(0);
  pad->SetGridy(1);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd();
  
//pad->SetRightMargin(0.06);
  //pad->SetLeftMargin(0.12);
  pad->SetObjectStat(0);
  pull->SetMarkerStyle(20);//20
	  
  pull->Draw("E");
  pull->Draw("SAME E");



  pull->GetYaxis()->SetTitle("data / prediction - 1");
  pull->GetXaxis()->SetTitle(xAxisTitle);

  //  pull->GetYaxis()->SetLabelOffset(0.035);
  pull->GetXaxis()->SetTitleSize(0.03);
  pull->GetXaxis()->SetLabelSize(0.03);
  pull->GetYaxis()->SetTitleSize(0.03);
  pull->GetYaxis()->SetLabelSize(0.03);
  pull->GetYaxis()->SetTitleOffset(1.5);
  pull->GetXaxis()->SetTitleOffset(1.5);
  pull->GetYaxis()->CenterTitle();

  /*  
  TPad * pad = new TPad("pad", "pad", 0.0, 0.0, 1.0, 1.0);
  pad->SetTopMargin(0.73);//0.72
  pad->SetBottomMargin(0.11);//0.72
  pad->SetFillColor(0);
  pad->SetGridy(1);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd();
  //pad->SetRightMargin(0.06);
  //pad->SetLeftMargin(0.12);
  pad->SetObjectStat(0);
  */
  
  pad->Draw();
  
  TLegend* leg1 = new TLegend(0.2, 0.7, 0.4, 0.9);
  
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetShadowColor(0);
  leg1->SetBorderSize(0);
  leg1->SetTextSize(3);
  
  
  //  TLegend* leg2 = new TLegend(0.7, 0.6, 0.9, 0.8);
  TLegend* leg2 = new TLegend(0.65, 0.72, 0.75, 0.92);
  
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetShadowColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  
  //this might be useful

  std::ostringstream nameDataDriven,nameWV,nameVVV, nameWZ, nameData, nameZZ, nameZgamma;
  int nbins_dd= h_DataDriven->GetNbinsX();

  nameDataDriven<<"data driven"; //("<<fixed<<setprecision(2)<< h_DataDriven->Integral(0, 1000)<<")";
  nameZZ<<"ZZ";// ("<<fixed<<setprecision(2)<< h_Zz->Integral(0,1000)<<")";
  nameZgamma<<"Zgamma"; //("<<  fixed<<setprecision(2)<< h_Zgamma->Integral(0,1000)<<")";
  nameWV<<"WV ";//("<<  fixed<<setprecision(2)<< h_Wv->Integral(0,1000)<<")";
  nameVVV<<"VVV ";//("<<  fixed<<setprecision(2)<< h_Vvv->Integral(0,1000)<<")";
  nameWZ<<"WZ ";//("<<  fixed<<setprecision(2)<< h_Wz->Integral(0,1000)<<")";

 if (channel==0)
   nameData<<"Data eee ";//("<<  fixed<<setprecision(2)<< h_data->Integral(0,1000)<<")";
  if (channel==1)
    nameData<<"Data ee#mu ";//("<<   fixed<<setprecision(2)<<h_data->Integral(0,1000)<<")";
  if (channel==2)
    nameData<<"Data e#mu#mu ";//("<<   fixed<<setprecision(2)<<h_data->Integral(0,1000)<<")";
  if (channel==3)
    nameData<<"Data #mu#mu#mu ";//("<<   fixed<<setprecision(2)<<h_data->Integral(0,1000)<<")";
  if (channel==4)
    nameData<<"Data ";//("<<  fixed<<setprecision(2)<< h_data->Integral(0,1000)<<")";
  /*
  nameDataDriven<<"data driven ("<< h_DataDriven->Integral(0, 500)<<")";
  nameZZ<<"ZZ ("<< h_Zz->Integral(0,500)<<")";
  nameZgamma<<"Zgamma ("<< h_Zgamma->Integral(0,500)<<")";
  nameWV<<"WV ("<< h_Wv->Integral(0,500)<<")";
  nameVVV<<"VVV ("<< h_Vvv->Integral(0,500)<<")";
  nameWZ<<"WZ ("<< h_Wz->Integral(0,500)<<")";
  */
  /* 
 if (channel==0)
    nameData<<"Data eee ("<< h_data->Integral(0,500)<<")";
  if (channel==1)
    nameData<<"Data ee#mu ("<< h_data->Integral(0,500)<<")";
  if (channel==2)
    nameData<<"Data e#mu#mu ("<< h_data->Integral(0,500)<<")";
  if (channel==3)
    nameData<<"Data #mu#mu#mu ("<< h_data->Integral(0,500)<<")";
  if (channel==4)
    nameData<<"Data ("<< h_data->Integral(0,500)<<")";
  */



  leg2->AddEntry(h_DataDriven, nameDataDriven.str().c_str());
  leg2->AddEntry(h_Zz, nameZZ.str().c_str());
  leg2->AddEntry(h_Zgamma, nameZgamma.str().c_str());
  leg2->AddEntry(h_Wv, nameWV.str().c_str());
  leg2->AddEntry(h_Vvv, nameVVV.str().c_str());
  leg2->AddEntry(h_Wz, nameWZ.str().c_str());
  leg2->AddEntry(h_data, nameData.str().c_str());
  
  /*
   leg2->AddEntry(h_Zjets, "zjets");
   leg2->AddEntry(h_Gvjets, "GVjets");
   leg2->AddEntry(h_Tt, "tt");
   leg2->AddEntry(h_Zzto2e2mu, "ZZ");
   leg2->AddEntry(h_Wz, "WZ");
   leg2->AddEntry(h_data, "DATA");
   */
  leg1->Draw();
  leg2->Draw();
  int iPos=11;
  int iPeriod=2; 
  writeExtraText = true; 
  lumi_8TeV  = "19.6 fb^{-1}";
  
  //  CMS_lumi(pad, iPeriod, iPos); 
  CMS_lumi(canv, iPeriod, iPos); 

   
   return 1;
   
}


void plotHistogramsNew()
{


  // Use CMS style: white background, no stat, ...
  /*
  gROOT->ProcessLine(".L ./CMSStyle.C");
  gROOT->ProcessLine("CMSstyle()");
  gStyle->SetOptStat(0);
  */
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  int iPos=11;
  int iPeriod=2; 
  writeExtraText = true; 
  lumi_8TeV  = "19.6 fb^{-1}";
/*    
  TFile *f1   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zjets.root");
  TFile *f2   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/top.root");
  TFile *f3   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/ZZ.root");
  TFile *f4   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zgamma.root");
  TFile *f5   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WV.root");
  TFile *f6   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/VVV.root");
  TFile *f7   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WZ.root");
  TFile *f8   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data.root");
  */
  TFile *f1   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data_driven.root");
  TFile *f2   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/ZZ.root");
  TFile *f3   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/Zgamma.root");
  TFile *f4   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WV.root");
  TFile *f5   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/VVV.root");
  TFile *f6   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/WZ.root");
  TFile *f7   = TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/data.root");
  
  const int n=5;

  TCanvas * c1[n];//, c2[4], c3[4], c4[4];
  TCanvas * c2[n];
  
  TCanvas * c3[n];
  TCanvas * c4[n];
  TCanvas * c5[n];
  TCanvas * c6[n];
  TCanvas * c7[n];
  TCanvas * c8[n];
  TCanvas * c9[n];
  TCanvas * c10[n];
  TCanvas * c11[n];
  TCanvas * c12[n];
  TCanvas * c13[n];
  TCanvas * c14[n];
  TCanvas * c15[n];
  TCanvas * c16[n];
  TCanvas * c17[n];
  TCanvas * c18[n];
  TCanvas * c19[n];
  TCanvas * c20[n];
  TCanvas * c21[n];
  TCanvas * c22[n];
  
  for (int canvas=0; canvas<n; canvas++){
    std::ostringstream name1, name2, name3, name4, name5, name6, name7, name8, name9, name10, name11, name12, name13, name14, name15, name16, name17, name18, name19, name20, name21, name22;
    name1<<"Zmass1_"<<canvas<<std::endl;
    name2<<"Zmass2_"<<canvas<<std::endl;
    name3<<"MET1_"<<canvas<<std::endl;
    name4<<"MET2_"<<canvas<<std::endl;
    name5<<"hZpt1_"<<canvas<<std::endl;
    name6<<"hZpt2_"<<canvas<<std::endl;
    name7<<"hLeadingJelPt1_"<<canvas<<std::endl;
    name8<<"hLeadingJelPt2_"<<canvas<<std::endl;
    name9<<"hNjets1_"<<canvas<<std::endl;
    name10<<"hNjets2_"<<canvas<<std::endl;
    name11<<"hDeltaPhi1_"<<canvas<<std::endl;
    name12<<"hDeltaPhi2_"<<canvas<<std::endl;
    name13<<"hZlepton1Pt1_"<<canvas<<std::endl;
    name14<<"hZlepton1Pt2_"<<canvas<<std::endl;
    name15<<"hWleptonPt1_"<<canvas<<std::endl;
    name16<<"hWleptonPt2_"<<canvas<<std::endl;
    name17<<"hMTW1_"<<canvas<<std::endl;
    name18<<"hMTW2_"<<canvas<<std::endl;
    name19<<"hZlepton2Pt1_"<<canvas<<std::endl;
    name20<<"hZlepton2Pt2_"<<canvas<<std::endl;
    name21<<"h3Lmass1_"<<canvas<<std::endl;
    name22<<"h3Lmass2_"<<canvas<<std::endl;

    c1[canvas]  =new TCanvas(name1.str().c_str(), name1.str().c_str());
    c2[canvas]  =new TCanvas(name2.str().c_str(), name2.str().c_str());
    c3[canvas]  =new TCanvas(name3.str().c_str(), name3.str().c_str());
    c4[canvas]  =new TCanvas(name4.str().c_str(), name4.str().c_str());
    c5[canvas]  =new TCanvas(name5.str().c_str(), name5.str().c_str());
    c6[canvas]  =new TCanvas(name6.str().c_str(), name6.str().c_str());
    c7[canvas]  =new TCanvas(name7.str().c_str(), name7.str().c_str());
    c8[canvas]  =new TCanvas(name8.str().c_str(), name8.str().c_str());
    c9[canvas]  =new TCanvas(name9.str().c_str(), name9.str().c_str());
    c10[canvas] =new TCanvas(name10.str().c_str(), name10.str().c_str());
    c11[canvas] =new TCanvas(name11.str().c_str(), name11.str().c_str());
    c12[canvas] =new TCanvas(name12.str().c_str(), name12.str().c_str());
    c13[canvas] =new TCanvas(name13.str().c_str(), name13.str().c_str());
    c14[canvas] =new TCanvas(name14.str().c_str(), name14.str().c_str());
    c15[canvas] =new TCanvas(name15.str().c_str(), name15.str().c_str());
    c16[canvas] =new TCanvas(name16.str().c_str(), name16.str().c_str());
    c17[canvas] =new TCanvas(name17.str().c_str(), name17.str().c_str());
    c18[canvas] =new TCanvas(name18.str().c_str(), name18.str().c_str());
    c19[canvas] =new TCanvas(name19.str().c_str(), name19.str().c_str());
    c20[canvas] =new TCanvas(name20.str().c_str(), name20.str().c_str());
    c21[canvas] =new TCanvas(name21.str().c_str(), name21.str().c_str());
    c22[canvas] =new TCanvas(name22.str().c_str(), name22.str().c_str());
    
    }
  
  std::ostringstream type;
  type<<"pdf";

  for (int histo=0; histo<n; histo++){
    std::ostringstream Zmass1, Zmass2, MET1, MET2, Zpt1, Zpt2, LeadingJetPt1, LeadingJetPt2, Njets1, Njets2, DeltaPhi1, DeltaPhi2, Zlepton1pt1, Zlepton1pt2, Zlepton2pt1, Zlepton2pt2, Wleptonpt1, Wleptonpt2, MTW1, MTW2, n3Lepmass1, n3Lepmass2;;
    std::ostringstream Zmass1save, Zmass2save, MET1save, MET2save, Zpt1save, Zpt2save, LeadingJetPt1save, LeadingJetPt2save, Njets1save, Njets2save, DeltaPhi1save, DeltaPhi2save, Zlepton1pt1save, Zlepton1pt2save, Zlepton2pt1save, Zlepton2pt2save, Wleptonpt1save, Wleptonpt2save, MTW1save, MTW2save, n3Lepmass1save, n3Lepmass2save;
    Zmass1<<"hZmass1_"<<histo;
    Zmass2<<"hZmass2_"<<histo;
    MET1<<"hMET1_"<<histo;
    MET2<<"hMET2_"<<histo;
    Zpt1<<"hZpt1_"<<histo;
    Zpt2<<"hZpt2_"<<histo;
    LeadingJetPt1<<"hLeadingJetPt1_"<<histo;
    LeadingJetPt2<<"hLeadingJetPt2_"<<histo;
    Njets1<<"hNjetsBigger1_"<<histo;
    Njets2<<"hNjetsBigger2_"<<histo;
    //    Njets1<<"hNjets1_"<<histo;
    //Njets2<<"hNjets2_"<<histo;
    DeltaPhi1<<"hDeltaPhi1_"<<histo;
    DeltaPhi2<<"hDeltaPhi2_"<<histo;
    Zlepton1pt1<<"hZlepton1pt1_"<<histo;
    Zlepton2pt1<<"hZlepton2pt1_"<<histo;
    Zlepton1pt2<<"hZlepton1pt2_"<<histo;
    Zlepton2pt2<<"hZlepton2pt2_"<<histo;
    Wleptonpt1<<"hWleptonpt1_"<<histo;
    Wleptonpt2<<"hWleptonpt2_"<<histo;
    MTW1<<"hMTW1_"<<histo;
    MTW2<<"hMTW2_"<<histo;
    n3Lepmass1<<"h3Lmass1_"<<histo;
    n3Lepmass2<<"h3Lmass2_"<<histo;
        
    c1[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zmass1.str().c_str(), histo, c1[histo], "M_{Z}(GeV)", 1);
    Zmass1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zmass1_"<<histo<<".pdf";
    //    CMS_lumi(c1[histo], iPeriod, iPos); 
    //    c1[histo]->Update();
    //c1[histo]->GetFrame()->Draw(); 
    //c1[histo]->Print(Zmass1save.str().c_str());
    c1[histo]->SaveAs(Zmass1save.str().c_str());
    
    
    c2[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zmass2.str().c_str(), histo, c2[histo], "M_{Z}(GeV)", 1);
    Zmass2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zmass2_"<<histo<<".pdf";
    //    CMS_lumi(c2[histo], iPeriod, iPos); 
    c2[histo]->SaveAs(Zmass2save.str().c_str());
    
    /*
    c3[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, MET1.str().c_str(), histo, c3[histo], "MET(GeV)", 5);
    MET1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/MET1_"<<histo<<".pdf";
//    CMS_lumi(c3[histo], iPeriod, iPos); 
    c3[histo]->SaveAs(MET1save.str().c_str());

    c4[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, MET2.str().c_str(), histo, c4[histo], "MET(GeV)", 5);
    MET2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/MET2_"<<histo<<".pdf";
    //CMS_lumi(c4[histo], iPeriod, iPos);
    c4[histo]->SaveAs(MET2save.str().c_str());
 
    c5[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zpt1.str().c_str(), histo, c5[histo], "Z_{pt}(GeV)", 10);
    Zpt1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zpt1_"<<histo<<".pdf";
    //CMS_lumi(c5[histo], iPeriod, iPos); 
    c5[histo]->SaveAs(Zpt1save.str().c_str());

    c6[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zpt2.str().c_str(),histo, c6[histo], "Z_{pt}(GeV)", 10);
    Zpt2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zpt2_"<<histo<<".pdf";
    //CMS_lumi(c6[histo], iPeriod, iPos); 
    c6[histo]->SaveAs(Zpt2save.str().c_str());
    
    c7[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, LeadingJetPt1.str().c_str(), histo, c7[histo], "LeadingJet_{pt}(GeV)", 5);
    LeadingJetPt1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/LeadingJetPt1_"<<histo<<".pdf";
    //CMS_lumi(c7[histo], iPeriod, iPos); 
    c7[histo]->SaveAs(LeadingJetPt1save.str().c_str());

    c8[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, LeadingJetPt2.str().c_str(), histo, c8[histo],"LeadingJet_{pt}(GeV)", 5);
    LeadingJetPt2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/LeadingJetPt2_"<<histo<<".pdf";
    //CMS_lumi(c8[histo], iPeriod, iPos); 
    c8[histo]->SaveAs(LeadingJetPt2save.str().c_str());
    
   
    c9[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Njets1.str().c_str(), histo, c9[histo], "N_{jets}", 1, false);
    Njets1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Njets1_"<<histo<<".pdf";
    //CMS_lumi(c9[histo], iPeriod, iPos); 
    c9[histo]->SaveAs(Njets1save.str().c_str());

    c10[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Njets2.str().c_str(),histo, c10[histo], "N_{jets}", 1, false);
    Njets2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Njets2_"<<histo<<".pdf";
    //CMS_lumi(c10[histo], iPeriod, iPos); 
    c10[histo]->SaveAs(Njets2save.str().c_str());
    */
    /*
    c11[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, DeltaPhi1.str().c_str(), histo, c11[histo], "#Delta #Phi", 0.1, false);
    DeltaPhi1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/DeltaPhi1_"<<histo<<".pdf";
    //CMS_lumi(c11[histo], iPeriod, iPos); 
    c11[histo]->SaveAs(DeltaPhi1save.str().c_str());

    c12[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, DeltaPhi2.str().c_str(), histo, c12[histo], "#Delta #Phi", 0.1, false);
    DeltaPhi2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/DeltaPhi2_"<<histo<<".pdf";
    //CMS_lumi(c12[histo], iPeriod, iPos); 
    c12[histo]->SaveAs(DeltaPhi2save.str().c_str());
    */
    /*
    c13[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zlepton1pt1.str().c_str(), histo, c13[histo], "p_{T}(GeV)", 10);
    Zlepton1pt1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zlepton1pt1_"<<histo<<".pdf";
    //CMS_lumi(c13[histo], iPeriod, iPos); 
    c13[histo]->SaveAs(Zlepton1pt1save.str().c_str());

    c14[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zlepton1pt2.str().c_str(), histo, c14[histo], "p_{T}(GeV)", 10);
    Zlepton1pt2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zlepton1pt2_"<<histo<<".pdf";
    //CMS_lumi(c14[histo], iPeriod, iPos); 
    c14[histo]->SaveAs(Zlepton1pt2save.str().c_str());

    c15[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Wleptonpt1.str().c_str(), histo, c15[histo], "p_{T}(GeV)", 10);
    Wleptonpt1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Wleptonpt1_"<<histo<<".pdf";
    //CMS_lumi(c15[histo], iPeriod, iPos); 
    c15[histo]->SaveAs(Wleptonpt1save.str().c_str());

    c16[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Wleptonpt2.str().c_str(), histo, c16[histo], "p_{T}(GeV)", 10);
    Wleptonpt2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Wleptonpt2_"<<histo<<".pdf";
    //CMS_lumi(c16[histo], iPeriod, iPos); 
    c16[histo]->SaveAs(Wleptonpt2save.str().c_str());
    
    
    c17[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, MTW1.str().c_str(), histo, c17[histo], "M_{T}^{W}(GeV)", 5);
    MTW1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/MTW1_"<<histo<<".pdf";
    //    CMS_lumi(c17[histo], iPeriod, iPos); 
    c17[histo]->SaveAs(MTW1save.str().c_str());
    
    c18[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, MTW2.str().c_str(), histo, c18[histo], "M_{T}^{W}(GeV)", 5);
    MTW2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/MTW2_"<<histo<<".pdf";
    //    CMS_lumi(c18[histo], iPeriod, iPos); 
    c18[histo]->SaveAs(MTW2save.str().c_str());

    c19[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zlepton2pt1.str().c_str(), histo, c19[histo], "p_{T}(GeV)", 10);
    Zlepton2pt1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zlepton2pt1_"<<histo<<".pdf";
    //CMS_lumi(c19[histo], iPeriod, iPos); 
    c19[histo]->SaveAs(Zlepton2pt1save.str().c_str());

    c20[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zlepton2pt2.str().c_str(), histo, c20[histo], "p_{T}(GeV)", 10);
    Zlepton2pt2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/Zlepton2pt2_"<<histo<<".pdf";
    // CMS_lumi(c20[histo], iPeriod, iPos); 
    c20[histo]->SaveAs(Zlepton2pt2save.str().c_str());
    
    
    c21[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, n3Lepmass1.str().c_str(), histo, c21[histo], "M(GeV)", 10);
    n3Lepmass1save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/invMassLep1_"<<histo<<".pdf";
    // CMS_lumi(c20[histo], iPeriod, iPos); 
    c21[histo]->SaveAs(n3Lepmass1save.str().c_str());

    c22[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, n3Lepmass2.str().c_str(), histo, c22[histo], "M(GeV)", 10);
    n3Lepmass2save<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/plotoviSrpanj/invMassLep2_"<<histo<<".pdf";
    // CMS_lumi(c20[histo], iPeriod, iPos); 
    c22[histo]->SaveAs(n3Lepmass2save.str().c_str());
    */
  }
 
 return;
}
