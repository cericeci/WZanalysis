#include "TFile.h"
#include "TString.h"
#include "/users/ltikvica/CMSSW_4_2_9_HLT1/src/WZanalysis/xsections.h"
#include <iostream>
#include <sstream>


TH1F *  DrawOverflow(int rebin, TH1F *h)
{
  // This function paint the histogram h with an extra bin for overflows

  const char* name  = h->GetName();
  const char* title = h->GetTitle();
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
  TH1F *htmp = new TH1F(name, title, nx, x1, x2);

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

void styleHisto1D(TH1F* histo, const char* titleX, double binWidth){
  histo->GetXaxis()->SetTitle(titleX);
  //  std::cout<<histo->GetXaxis()<<std::endl;
  histo->GetYaxis()->SetTitle(Form("Number of Events / %1.2f GeV", binWidth));

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
  histo->GetYaxis()->SetTitleOffset(1.);
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
		  TString xAxisTitle="",
		  double binWidth=0)
		  
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
 
  

  THStack *  th = new THStack(histKey,histKey);
  th->Add(h_Wv); 
  th->Add(h_Vvv);    
  th->Add(h_Zz);
  th->Add(h_DataDriven);
  th->Add(h_Wz);   
  th->Add(h_Zgamma);
 

  
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


  //c01->RedrawAxis();
  
  float dataMax = h_data->GetMaximum();
  float mcMax =   th->GetMaximum();
  /*  
  std::cout<<"data max:"<<dataMax<<std::endl;
  std::cout<<"MC max:"<< mcMax<<std::endl;
  */
  
  if (dataMax > mcMax) {
    th->SetMaximum(dataMax);
  }
  
  //else {th->SetMaximum(mcMax);}
  float hmin = th->GetMinimum();
  
  if (hmin>0.) th->SetMinimum(0.5*hmin);
  
  TAxis * a =  th->GetXaxis(); // ->SetTitle(xAxisTitle);

  h_data->SetMarkerStyle(20);
  h_data->SetLabelFont(62);
  //  gPad->SetBottomMargin(0.28);
  gPad->SetBottomMargin(0.29);
  
  styleHisto1D(h_data, xAxisTitle,binWidth);
  
  //th->Draw();
  th->Draw("SAME");  
  
  //  h_data->Draw();  
  h_data->Draw("SAMEPE1");  
  
  plotLatex_B(19.602);

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
  pull->GetYaxis()->SetTitleOffset(1.0);
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
  
  
  TLegend* leg2 = new TLegend(0.7, 0.6, 0.9, 0.8);
  
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetShadowColor(0);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  
  //this might be useful

  std::ostringstream nameDataDriven,nameWV,nameVVV, nameWZ, nameData, nameZZ, nameZgamma;
  nameDataDriven<<"data driven ("<< h_DataDriven->Integral(0, 500)<<")";
  nameZZ<<"ZZ ("<< h_Zz->Integral(0,500)<<")";
  nameZgamma<<"Zgamma ("<< h_Zgamma->Integral(0,500)<<")";
  nameWV<<"WV ("<< h_Wv->Integral(0,500)<<")";
  nameVVV<<"VVV ("<< h_Vvv->Integral(0,500)<<")";
  nameWZ<<"WZ ("<< h_Wz->Integral(0,500)<<")";
  nameData<<"Data ("<< h_data->Integral(0,500)<<")";
  

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
   
   return 1;
   
}


void plotHistogramsNew()
{


  // Use CMS style: white background, no stat, ...
  gROOT->ProcessLine(".L ./CMSStyle.C");
  gROOT->ProcessLine("CMSstyle()");
  gStyle->SetOptStat(0);
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
  
  TCanvas * c1[4];
  
  for (int canvas=0; canvas<4; canvas++){
    std::ostringstream name;
    name<<"c1_"<<canvas<<std::endl;
    c1[canvas]=new TCanvas(name.str().c_str(), name.str().c_str());
  }
  
  for (int histo=0; histo<4; histo++){
    std::ostringstream Zmass1, Zmass2;
    Zmass1<<"hZmass1_"<<histo;
    c1[histo]->cd();
    plotDataVsMC(f1,f2,f3,f4,f5,f6,f7, Zmass1.str().c_str(),"M_{Z}(GeV)", 1.2);
  }
 


 return;
}
