#include "TROOT.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPad.h"
#include "TMath.h"
#include "TColor.h"


const Font_t _cmsTextFont   =  61;
const Font_t _extraTextFont =  52;
const Font_t _lumiTextFont  =  42;
const Size_t _msize         = 1.1;


void drawFigure4bc(TString var="Njets") 
{
  gInterpreter->ExecuteMacro("GoodStyle.C");
  gROOT->LoadMacro("tdrstyle.C");                                                                                                                           
  setTDRStyle();    

  TString variable = "p_{T,max}^{#font[12]{l}}";
  std::ostringstream filename, madname;
  filename<<"rootfiles/all_unfolding_"<<var<<".root";
  TFile* file = new TFile(filename.str().c_str(), "read");
  madname<<"hGenXs"<<var<<"_1";
  std::cout<<madname.str().c_str()<<std::endl;
  TH1F* xsValue          = (TH1F*)(file->Get("hComb_diff")->Clone("xsValue"));
  TH1F* xsValue_Madgraph = (TH1F*)(file->Get(madname.str().c_str())->Clone("xsValue_Madgraph"));
  //  TH1F* xsValue_MCnlo    = (TH1F*)(file->Get("mcfm_tot")->Clone("xsValue_MCnlo"));

  // Set the data errors- I don't need this because I already have complete error in my plot
  //----------------------------------------------------------------------------
  
  // Data cosmetics
  //----------------------------------------------------------------------------
  xsValue->SetLineWidth(1);
  xsValue->SetMarkerSize(_msize);
  xsValue->SetMarkerStyle(kFullCircle);
  xsValue->SetMarkerColor(kBlack);
  xsValue->SetLineColor(kBlack);
  xsValue->SetFillStyle(1001);
  xsValue->SetFillColor(kWhite);
  // Madgraph cosmetics
  //----------------------------------------------------------------------------
  xsValue_Madgraph->SetFillColor(kOrange);
  //xsValue_Madgraph->SetFillColor(kWhite);
  xsValue_Madgraph->SetFillStyle(1001);
  xsValue_Madgraph->SetLineColor(kOrange+7);
  xsValue_Madgraph->SetLineWidth(1);
  xsValue_Madgraph->SetMarkerColor(kOrange+7);
  xsValue_Madgraph->SetMarkerSize(_msize);
  xsValue_Madgraph->SetMarkerStyle(21);


  // MCNLO cosmetics
  //----------------------------------------------------------------------------
  //xsValue_MCnlo->SetFillColor(kAzure-9);
   //xsValue_MCnlo->SetFillColor(kWhite);
  //xsValue_MCnlo->SetFillStyle(1001);
  //xsValue_MCnlo->SetLineColor(kAzure);
  //xsValue_MCnlo->SetLineWidth(1);
  //xsValue_MCnlo->SetMarkerColor(kAzure);
  //xsValue_MCnlo->SetMarkerSize(_msize);
  //  xsValue_MCnlo->SetMarkerStyle(21);
  //xsValue_MCnlo->SetMarkerStyle(24);

  //  TCanvas * c1=new TCanvas("c1", "c1");
  
//xsValue_MCnlo->Draw("pey0");
   

  // Set the canvas and pads
  //----------------------------------------------------------------------------
  //  TCanvas* canvas = new TCanvas("wwxs", "wwxs", 600, 600);
  TCanvas* canvas = new TCanvas("wwxs", "wwxs"); 

  //defalut
  //TCanvas* canvas = new TCanvas("wwxs", "wwxs", 600, 850);
  //  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.55, 1, 1.000);
  //TPad* pad2 = new TPad("pad2", "pad2", 0, 0.39, 1, 0.552);
  //TPad* pad3 = new TPad("pad3", "pad3", 0, 0.23, 1, 0.392);

  TPad* pad1 = new TPad("pad1", "pad1", 0, 0.332, 1, 1.000);
  TPad* pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.33);
  //TPad* pad3 = new TPad("pad3", "pad3", 0, 0, 1, 0.332);
  
  pad1->SetTopMargin(0.09);
  pad2->SetTopMargin(0);
  //  pad3->SetTopMargin(0);

  pad1->SetBottomMargin(0);
  pad2->SetBottomMargin(0.45);
  //  pad3->SetBottomMargin(0.15);
  //pad3->SetBottomMargin(0.45);

  pad1->SetLeftMargin(0.16);
  pad2->SetLeftMargin(0.16);
  //pad3->SetLeftMargin(0.16);

  pad1->SetRightMargin(0.06);
  pad2->SetRightMargin(0.06);
  //pad3->SetRightMargin(0.06);



  // pad1
  //----------------------------------------------------------------------------
  pad1->Draw();
  pad1->cd();
  pad1->SetLogy();


  // Draw
  //----------------------------------------------------------------------------
  AxisFonts(xsValue->GetXaxis(), variable + " (GeV)");


  if (var=="LeadingJetPt")
    AxisFonts(xsValue->GetYaxis(), "d#sigma(WZ#rightarrow3l#nu)/dp_{T}^{LeadingJet}");
  if (var=="Njets")
    AxisFonts(xsValue->GetYaxis(), "d#sigma(WZ#rightarrow3l#nu)/dN_{jets}");

  
  //TH1F* hpowError = (TH1F*)xsValue_Powheg->Clone();
  //  TH1F* hmadError = (TH1F*)xsValue_Madgraph->Clone();
  //TH1F* hmcError  = (TH1F*)xsValue_MCnlo->Clone();

  xsValue         ->Draw("pe");
  xsValue_Madgraph->Draw("p,same");
  //xsValue_Madgraph->Draw("p");
  
  //  xsValue_MCnlo   ->Draw("p,same");



  if (var=="LeadingJetPt")
    xsValue->SetMinimum(5e-4);
  if (var=="Njets")
    xsValue->SetMinimum(0.02);

 
  // Legend
  //----------------------------------------------------------------------------
  DrawLegend(0.718, 0.80, xsValue,   " Data",     "lp");
  DrawLegend(0.718, 0.74, xsValue_Madgraph, " Madgraph", "flp");  
  //  DrawLegend(0.718, 0.68, xsValue_MCnlo,  " MC@NLO",   "flp");


  // Draw text 
  //----------------------------------------------------------------------------
  DrawLatex(_cmsTextFont,   0.173, 0.935, 0.065, 11, "CMS");
  DrawLatex(_extraTextFont, 0.268, 0.935, 0.035, 11, "Preliminary");
  DrawLatex(_lumiTextFont,  0.940, 0.935, 0.050, 31, "19.6 fb^{-1} (8 TeV)");


  // Prepare the ratios
  //----------------------------------------------------------------------------
  TH1F* ratio_mad    = xsValue_Madgraph->Clone("ratio");
  //  TH1F* ratio_mcnlo  = xsValue_Madgraph->Clone("ratio");
  //  TH1F* ratio_mcnlo  = xsValue_MCnlo->Clone("ratio");
  TH1F* hratio_mad   = xsValue_Madgraph->Clone("ratio");
  //TH1F* hratio_mcnlo = xsValue_MCnlo->Clone("ratio");
  TH1F* ratioErr     = xsValue->Clone("ratio");


  ratioErr->SetFillColor  (kGray+2);
  ratioErr->SetFillStyle  (   3004);
  ratioErr->SetLineColor  (kGray+2);
  ratioErr->SetMarkerColor(kGray+2);
  ratioErr->SetMarkerSize (      0);


  // Set the bin content
  //----------------------------------------------------------------------------
  for (UInt_t ibin=1; ibin<=ratio->GetNbinsX(); ibin++) {
   
    Double_t madValue = xsValue_Madgraph->GetBinContent(ibin);
    //Double_t madError = xsValue_Madgraph->GetBinError  (ibin);
   
    //Double_t mcnloValue = xsValue_MCnlo->GetBinContent(ibin);
    //Double_t mcnloError = xsValue_MCnlo->GetBinError  (ibin);
   
    Double_t dataValue = xsValue->GetBinContent(ibin);
    
    Double_t dataError = xsValue->GetBinError(ibin);
   
    Double_t ratioValue_mad = (madValue > 0) ? madValue / dataValue : 0.0;
    //Double_t ratioError_mad = (madValue > 0) ? madError / dataValue : 0.0;
    Double_t ratioError_mad = madValue/pow(dataValue,2)*dataError;
    //Double_t ratioValue_mcnlo = (mcnloValue > 0) ? mcnloValue / dataValue : 0.0;
    // Double_t ratioError_mcnlo = (mcnloValue > 0) ? mcnloError / dataValue : 0.0;
    //Double_t ratioError_mcnlo = mcnloValue/pow(dataValue,2)*dataError;
   
    Double_t uncertaintyError = (dataValue > 0) ? dataError / dataValue : 0.0;
   
    ratio_mad ->SetBinContent(ibin, ratioValue_mad);
    hratio_mad->SetBinContent(ibin, ratioValue_mad);
    hratio_mad->SetBinError  (ibin, ratioError_mad);
   
    //ratio_mcnlo ->SetBinContent(ibin, ratioValue_mcnlo);
    //hratio_mcnlo->SetBinContent(ibin, ratioValue_mcnlo);
    //hratio_mcnlo->SetBinError  (ibin, ratioError_mcnlo);
   
    ratioErr->SetBinContent(ibin, 1.0);
    ratioErr->SetBinError  (ibin, uncertaintyError);  //??? ovo nije bas jasno sto je
  }


  //  AxisFontsRatio(ratioErr->GetYaxis(), "y", "Theory / Data");
  //AxisFontsRatio(ratioErr->GetXaxis(), "x", variable + " (GeV)");
  AxisFontsRatio(ratio_mad->GetYaxis(), "y", "Theory / Data");

  if (var=="LeadingJetPt")
    AxisFontsRatio(ratio_mad->GetXaxis(), "x", "p_{T}^{LeadingJet} (GeV)");
  if (var=="Njets")
    AxisFontsRatio(ratio_mad->GetXaxis(), "x", "N_{jets}");

  //  AxisFontsRatio(ratio_mcnlo->GetYaxis(), "y", "Theory / Data");
  //AxisFontsRatio(ratio_mcnlo->GetXaxis(), "x", "p_{T}^{Z} (GeV)");


  //  ratio_mcnlo->SetFillColor(kAzure-9);
   //xsValue_MCnlo->SetFillColor(kWhite);
  //ratio_mcnlo->SetFillStyle(1001);
  //ratio_mcnlo->SetLineColor(kAzure);
  //ratio_mcnlo->SetLineWidth(1);
  //ratio_mcnlo->SetMarkerColor(kAzure);
  //ratio_mcnlo->SetMarkerSize(_msize);
  //  xsValue_MCnlo->SetMarkerStyle(21);
  //ratio_mcnlo->SetMarkerStyle(24);

  // Draw pad2
  //----------------------------------------------------------------------------
  canvas->cd();
  pad2->Draw();
  pad2->cd();
  
  //ratioErr  ->Draw("e2");
  
  ratio_mad ->Draw("pz");
  hratio_mad->Draw("e2,same");
  ratio_mad ->Draw("p, same");
  ratio_mad->GetYaxis()->SetRangeUser(0.01, 2.3);  

  pad2->Modified();
  
  DrawLatex(43, 0.2, 0.85, 15.0, 11, "Madgraph+Pythia normalized to #sigma_{NNLO}");
  

  // Draw pad3
  //----------------------------------------------------------------------------
  //canvas->cd();
  //pad3->Draw();
  //pad3->cd();
  //std::cout<<"Default option: "<<ratio_mcnlo->GetOption()<<std::endl;
  
  //ratioErr    ->Draw("e2");
  //ratio_mcnlo->Draw("pz");
  //hratio_mcnlo->Draw("e2,same");

  //  ratio_mcnlo ->Draw("pz, same");
  //ratio_mcnlo->GetYaxis()->SetRangeUser(0.1, 2.3);  
  //ratio_mcnlo->GetYaxis()->SetRangeUser(0.0, 0.5);  
  
  //pad3->Modified();

  //DrawLatex(43, 0.2, 0.89, 15.0, 11, "MC@NLO+Herwig normalized to #sigma_{NNLO}");
  
    

  // Save
  //----------------------------------------------------------------------------
  pad1->cd(); pad1->GetFrame()->DrawClone();
  pad2->cd(); pad2->GetFrame()->DrawClone();
  //pad3->cd(); pad3->GetFrame()->DrawClone();

  canvas->cd();

  std::ostringstream saveName, saveName2;
  saveName <<"pdf/unfolded_" << var << ".pdf"; 
  saveName2<<"png/unfolded_" << var << ".png";
  canvas->SaveAs(saveName.str().c_str());
  canvas->SaveAs(saveName2.str().c_str());
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFonts(TAxis* axis, TString title)
{
  axis->SetLabelFont  (   43);
  axis->SetLabelOffset(0.007);
  axis->SetLabelSize  (   22);
  axis->SetNdivisions (  110);
  axis->SetTitle      (title);
  axis->SetTitleFont  (   43);
  //  axis->SetTitleOffset(  2.5);
  axis->SetTitleOffset(  1.5);
  axis->SetTitleSize  (   22);
}


//------------------------------------------------------------------------------
// AxisFonts
//------------------------------------------------------------------------------
void AxisFontsRatio(TAxis*  axis,
		    TString which,
		    TString title)
{
  if (which.Contains("x"))
    {
      axis->SetLabelFont  (   43);
      axis->SetLabelOffset(0.022);
      axis->SetLabelSize  (   22);
      axis->SetNdivisions (  510);
      axis->SetTitle      (title);
      axis->SetTitleFont  (   43);
      axis->SetTitleOffset(5.0);
      axis->SetTitleSize  (   22);
    }
  else
    {
      axis->CenterTitle   ( true);
      axis->SetLabelFont  (   43);
      axis->SetLabelOffset(0.012);
      axis->SetLabelSize  (   15);
      axis->SetNdivisions (  505);
      axis->SetTitle      (title);
      axis->SetTitleFont  (   43);
      axis->SetTitleOffset(  1.5);
      axis->SetTitleSize  (   15);
    }
}


//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
void DrawLegend(Float_t x1,
		Float_t y1,
		TH1F*   hist,
		TString label,
		TString option)
{
  TLegend* legend = new TLegend(x1,
				y1,
				x1 + 0.22,
				y1 + 0.05);

  legend->SetBorderSize(    0);
  legend->SetFillColor (    0);
  legend->SetTextAlign (   12);
  legend->SetTextFont  (   42);
  legend->SetTextSize  (0.045);

  legend->AddEntry (hist, label.Data(), option.Data());
  //  legend->AddEntry (hist, label.Data());
  legend->Draw();
}


//------------------------------------------------------------------------------
// DrawLatex
//------------------------------------------------------------------------------
void DrawLatex(Font_t      tfont,
	       Double_t    x,
	       Double_t    y,
	       Double_t    tsize,
	       Short_t     align,
	       const char* text,
	       Bool_t      setndc = true)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC      (setndc);
  tl->SetTextAlign( align);
  tl->SetTextFont ( tfont);
  tl->SetTextSize ( tsize);

  tl->Draw("same");
}
