#include "TFile.h"
#include "TString.h"
//#include "xsections.h"
#include <iostream>
#include <sstream>
#include "TH1F.h"
#include "TPad.h" 
#include "TLatex.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TROOT.h"

void plotFinal(TString variable, float adjustMin=-999)
{
  //gROOT->LoadMacro("tdrstyle.C");
  //  setTDRStyle();

  gROOT->ProcessLine(".L CMSStyle.C");
  CMSstyle();
	
  gROOT->LoadMacro("CMS_lumi.C");	
  
  int iPos=22;
  int iPeriod=2;
  writeExtraText = true;
  lumi_8TeV  = "19.6 fb^{-1}";

  bool plotAllChannels = true;
  bool plotMadgraph = true;
  bool plotMCFM     = true;
  // No point in plotting jet multiplicity for MCFM
  if (variable == "Njets") plotMCFM = false;

  TFile * fmg;
  TFile * fmcfm_wp;
  TFile * fmcfm_wm;


  std::ostringstream fileName, outputName;
  fileName<<"unfoldingFinalResults/combination_"<<variable<<".root";
  outputName<<"unfoldingFinalResults/all_unfolding_"<<variable<<".root";
  TString mg_fileName = "unfoldingFinalResults/mg_full.root";

  TFile *f1  = TFile::Open(fileName.str().c_str());
  TFile * fout= new TFile(outputName.str().c_str(), "RECREATE");

  const int nChannels;
  std::ostringstream histName;
  TH1F * crossSection[nChannels];
  TH1F * hComb= (TH1F*) (f1->Get("h_xs_comb") ->Clone("hComb"));
  TH1F * hComb_diff= (TH1F*) (f1->Get("h_xs_comb_diff") ->Clone("hComb_diff"));
  TH1F * hcrossSection_0= (TH1F*) (f1->Get("h_crossSection_inclusive0") ->Clone("h_crossSection_inclusive0"));
  TH1F * hcrossSection_1= (TH1F*) (f1->Get("h_crossSection_inclusive1") ->Clone("h_crossSection_inclusive1"));
  TH1F * hcrossSection_2= (TH1F*) (f1->Get("h_crossSection_inclusive2") ->Clone("h_crossSection_inclusive2"));
  TH1F * hcrossSection_3= (TH1F*) (f1->Get("h_crossSection_inclusive3") ->Clone("h_crossSection_inclusive3"));

  TH1F * hcrossSection_diff_0= (TH1F*) (f1->Get("h_crossSection_incl_diff0") ->Clone("h_crossSection_incl_diff0"));
  TH1F * hcrossSection_diff_1= (TH1F*) (f1->Get("h_crossSection_incl_diff1") ->Clone("h_crossSection_incl_diff1"));
  TH1F * hcrossSection_diff_2= (TH1F*) (f1->Get("h_crossSection_incl_diff2") ->Clone("h_crossSection_incl_diff2"));
  TH1F * hcrossSection_diff_3= (TH1F*) (f1->Get("h_crossSection_incl_diff3") ->Clone("h_crossSection_incl_diff3"));


  // Madgraph plots

  TH1D * mg_histos[4];

  if (plotMadgraph) {
    fmg = new TFile(mg_fileName);
    for (int i=0; i<4; i++) {
      TString mg_hname = "hGenXs";
      mg_hname += variable;
      mg_hname += "_";
      mg_hname += i+1;
      
      mg_histos[i] = (TH1D*) (fmg->Get(mg_hname));
      mg_histos[i]->SetLineWidth(3);
      mg_histos[i]->Print();
    }
  }

  // MCFM plots

  TH1D * mcfmPlots[2];
  if (plotMCFM) {

    TFile * fmcfm_wp = new TFile("unfoldingFinalResults/mcfm-plots-71.root");
    TFile * fmcfm_wm = new TFile("unfoldingFinalResults/mcfm-plots-76.root");
    
    TString mcfmPlotName = "h"+variable;
    std::cout << "MCFM histo : " << mcfmPlotName;
    
    mcfmPlots[0] = (TH1D*) fmcfm_wp->Get(mcfmPlotName)->Clone("mcfm_tot");
    mcfmPlots[1] = (TH1D*) fmcfm_wm->Get(mcfmPlotName)->Clone("mcfm_tot");
    mcfmPlots[0]->Add(mcfmPlots[1]);
    
    mcfmPlots[0]->Print();
    mcfmPlots[0]->SetLineColor(4);
    mcfmPlots[0]->SetLineWidth(3);
  }

  TString yAxisName;
  if (variable=="Njets") yAxisName="d#sigma/dN";
  else yAxisName="d#sigma/dP_{T}";
  TLegend* leg1= new TLegend(0.7, 0.7, 0.9, 0.9);
  //  TLegend* leg1= new TLegend(0.6, 0.7, 0.6, 0.9);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->AddEntry(hComb, "Combination");
  if (plotAllChannels) {
    leg1->AddEntry(hcrossSection_0, "3e");
    leg1->AddEntry(hcrossSection_1, "2e1mu");
    leg1->AddEntry(hcrossSection_2, "1e2mu");
    leg1->AddEntry(hcrossSection_3, "3mu");
  }
  if (plotMadgraph) leg1->AddEntry(mg_histos[0], "MADGRAPH (NLO norm.)");
  if (plotMCFM)     leg1->AddEntry(mcfmPlots[0], "MCFM");

  
  TCanvas * c1= new TCanvas("crossSection", "crossSection"); 
  hComb->SetLineColor(kRed);
  //  hComb->SetMarkerSize(1000);
  hComb->SetMarkerStyle(21);
  hComb->SetMarkerColor(kRed);
  hComb->SetFillColor(2);
  hComb->SetFillStyle(3001);
  hComb->GetXaxis()->SetTitle(variable);
  hComb->GetYaxis()->SetTitle("#sigma");
  hcrossSection_0->SetLineColor(kBlue);
  hcrossSection_0->SetMarkerColor(kBlue);
  hcrossSection_0->SetMarkerStyle(20);
  
  hcrossSection_1->SetLineColor(kYellow);
  hcrossSection_1->SetMarkerColor(kYellow);
  hcrossSection_1->SetMarkerStyle(20);
  
  hcrossSection_2->SetLineColor(kGreen);
  hcrossSection_2->SetMarkerColor(kGreen);
  hcrossSection_2->SetMarkerStyle(20);

  hcrossSection_3->SetLineColor(kViolet);
  hcrossSection_3->SetMarkerColor(kViolet);
  hcrossSection_3->SetMarkerStyle(20);
  
  hComb->Draw("E2");
  hcrossSection_0->Draw("SAMEPE");
  hcrossSection_1->Draw("SAMEPE");
  hcrossSection_2->Draw("SAMEPE");
  hcrossSection_3->Draw("SAMEPE");
  
  leg1->Draw();
  std::ostringstream saveName;
  saveName<<"unfoldingFinalResults/combinationDsigma_"<<variable<<".pdf";
  //  c1->SetLogy();
  CMS_lumi(c1, iPeriod, iPos); 
  c1->SaveAs(saveName.str().c_str());
  

  TCanvas * c2= new TCanvas("crossSection_diff", "crossSection_diff"); 
  hComb_diff->SetLineColor(kRed);
  //  hComb->SetMarkerSize(1000);
  hComb_diff->SetMarkerStyle(21);
  hComb_diff->SetMarkerColor(kRed);
  hComb_diff->SetMarkerColor(kRed);
  hComb_diff->SetFillColor(2);  
  hComb_diff->SetFillStyle(3001);
  hComb_diff->GetXaxis()->SetTitle(variable);
  hComb_diff->GetYaxis()->SetTitle(yAxisName);

  if (adjustMin > 0.) {
    hComb_diff->SetMinimum(adjustMin*hComb_diff->GetMinimum());
  }

  
  hcrossSection_diff_0->SetLineColor(kBlue);
  hcrossSection_diff_0->SetMarkerColor(kBlue);
  hcrossSection_diff_0->SetMarkerStyle(20);
  
  hcrossSection_diff_1->SetLineColor(kYellow);
  hcrossSection_diff_1->SetMarkerColor(kYellow);
  hcrossSection_diff_1->SetMarkerStyle(20);
  
  hcrossSection_diff_2->SetLineColor(kGreen);
  hcrossSection_diff_2->SetMarkerColor(kGreen);
  hcrossSection_diff_2->SetMarkerStyle(20);

  hcrossSection_diff_3->SetLineColor(kViolet);
  hcrossSection_diff_3->SetMarkerColor(kViolet);
  hcrossSection_diff_3->SetMarkerStyle(20);
  
  hComb_diff->Draw("E2");
  if (plotAllChannels) {
    hcrossSection_diff_0->Draw("SAMEPE");
    hcrossSection_diff_1->Draw("SAMEPE");
    hcrossSection_diff_2->Draw("SAMEPE");
    hcrossSection_diff_3->Draw("SAMEPE");
  }

  if (plotMadgraph) mg_histos[0]->Draw("SAMEH");
  if (plotMCFM)     mcfmPlots[0]->Draw("SAMEH");

  leg1->Draw();
  std::ostringstream saveNameDiff, saveNameDiffLog;
  saveNameDiff<<"unfoldingFinalResults/combinationDsigma_diff_"<<variable<<".pdf";
  CMS_lumi(c2, iPeriod, iPos);    
  c2->SaveAs(saveNameDiff.str().c_str());
  c2->SetLogy();
  saveNameDiffLog<<"unfoldingFinalResults/combinationDsigma_diff_"<<variable<<"_LOG.pdf";
  c2->SaveAs(saveNameDiffLog.str().c_str());

  //ovo se treba zakomentirati pa onda nece crtati gluposti tu prije
  fout->cd();
  hComb_diff->Write();
  mg_histos[0]->Write();
  if (variable!="Njets")
    mcfmPlots[0]->Write();
  fout->Close();
}

