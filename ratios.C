const UInt_t nchannels = 5;


Double_t xs7tev         [nchannels] = {18.8907, 19.474, 22.436, 16.974, 18.737};
Double_t xs7tevErrorStat[nchannels] = { 1.2036,  2.188,  2.991,  2.122,  2.584};
Double_t xs7tevErrorSyst[nchannels] = { 0.6947,  0.732,  0.966,  0.661,  0.729};
Double_t xs7tevErrorLumi[nchannels] = { 0.4396,  0.457,  0.512,  0.394,  0.438};

Double_t xs8tev         [nchannels] = {23.8343, 24.6729, 22.2613  , 23.8032 , 24.8118};
Double_t xs8tevErrorStat[nchannels] = { 0.7792,  1.9209,  1.61596 ,  1.52062,  1.29358};
Double_t xs8tevErrorSyst[nchannels] = { 1.0259,  1.36,  1.15759,  1.42819,  1.38946};
Double_t xs8tevErrorLumi[nchannels] = { 0.62,  0.64 , 0.58,  0.62,  0.65};

TString label[nchannels] = {"inclusive", "#mu#mu#mu", "eee", "#mu#mue", "ee#mu"};


// Settings
//------------------------------------------------------------------------------
TString _format = "pdf";


//------------------------------------------------------------------------------
// ratios
//------------------------------------------------------------------------------
void ratios(Int_t ecm = 8)
{
  Double_t luminosity = (ecm == 8) ? 19602 : 4920;  // pb-1
  Double_t xsWplusZ   = (ecm == 8) ?  13.9 : 11.4;  // pb
  Double_t xsWminusZ  = (ecm == 8) ?   8.1 :  6.4;  // pb

  gInterpreter->ExecuteMacro("HiggsPaperStyle.C");
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  gROOT->LoadMacro("CMS_lumi.C");

  gStyle->SetEndErrorSize(5);

  gSystem->mkdir(_format, kTRUE);

  int iPos=22;
  int iPeriod=2; 
  writeExtraText = true; 
  lumi_8TeV  = "19.6 fb^{-1}";

  // Loop
  //----------------------------------------------------------------------------
  TGraphErrors* gStat = new TGraphErrors(nchannels);
  TGraphErrors* gSyst = new TGraphErrors(nchannels);
  TGraphErrors* gLumi = new TGraphErrors(nchannels);

  for (UInt_t i=0; i<nchannels; i++) {

    Double_t xs          = (ecm == 8) ? xs8tev         [i] : xs7tev         [i];
    Double_t xsErrorStat = (ecm == 8) ? xs8tevErrorStat[i] : xs7tevErrorStat[i];
    Double_t xsErrorSyst = (ecm == 8) ? xs8tevErrorSyst[i] : xs7tevErrorSyst[i];
    Double_t xsErrorLumi = (ecm == 8) ? xs8tevErrorLumi[i] : xs7tevErrorLumi[i];

    Double_t f = xs / (xsWplusZ + xsWminusZ);

    Double_t errorSquared = (xsErrorStat * xsErrorStat);

    gStat->SetPointError(i, f * sqrt(errorSquared) / xs, 0.0);

    errorSquared += (xsErrorSyst * xsErrorSyst);

    gSyst->SetPointError(i, f * sqrt(errorSquared) / xs, 0.0);

    errorSquared += (xsErrorLumi * xsErrorLumi);

    gLumi->SetPointError(i, f * sqrt(errorSquared) / xs, 0.0);

    gStat->SetPoint(i, f, i+1);
    gSyst->SetPoint(i, f, i+1);
    gLumi->SetPoint(i, f, i+1);
  }


  // Cosmetics
  //----------------------------------------------------------------------------
  gStat->SetLineWidth  (2);
  gStat->SetMarkerSize (1.3);
  //  gStat->SetMarkerSize (2.3);
  gStat->SetMarkerStyle(kFullCircle);

  gSyst->SetLineColor  (kRed);
  gSyst->SetLineWidth  (2);
  gSyst->SetMarkerSize (1.3);
  gSyst->SetMarkerStyle(kFullCircle);

  gLumi->SetLineColor  (kBlue);
  gLumi->SetLineWidth  (2);
  gLumi->SetMarkerSize (1.3);
  gLumi->SetMarkerStyle(kFullCircle);


  // Draw
  //----------------------------------------------------------------------------
  TCanvas* canvas = new TCanvas();

  canvas->SetLeftMargin(canvas->GetRightMargin());

  Double_t xmin = 0.6;
  Double_t xmax = 2.2;
  Double_t ymin = 0.50;
  Double_t ymax = nchannels + ymin;
  
  TH2F* dummy = new TH2F("dummy", "", 100, xmin, xmax, 100, ymin, ymax);

  dummy->Draw();


  // Vertical line at 1
  //----------------------------------------------------------------------------
  TLine* line = new TLine(1.0, ymin, 1.0, ymax);

  line->SetLineWidth(2);

  line->Draw("same");


  // Ratios
  //----------------------------------------------------------------------------
  gLumi->Draw("p||,same");
  gSyst->Draw("p||,same");
  gStat->Draw("p,same");


  // Labels
  //----------------------------------------------------------------------------
  for (UInt_t i=0; i<nchannels; i++) {

    Double_t x = gStat->GetX()[i];
    Double_t y = gStat->GetY()[i];

    DrawTLatex(xmin+0.05, y, 0.035, 12, Form("%s", label[i].Data()));

    Double_t statError  = gStat->GetErrorX(i);
    Double_t systError  = gSyst->GetErrorX(i);
    Double_t lumiError  = gLumi->GetErrorX(i);

    lumiError = sqrt(lumiError*lumiError - systError*systError);

    systError = sqrt(systError*systError - statError*statError);

    DrawTLatex(xmax-0.05, y, 0.035, 32, Form("%.2f #pm %.2f #pm %.2f #pm %.2f",
					     x, statError, systError, lumiError));
  }

  //  DrawTLatex(0.940, 0.983, 0.05, 33,
	     //	     Form("#sqrt{s} = %d TeV, L = %.1f fb^{-1}", ecm, luminosity/1e3), true);

  dummy->GetXaxis()->CenterTitle();
  dummy->GetXaxis()->SetTitleOffset(0.9);
  dummy->GetXaxis()->SetTitleSize(0.040);
  dummy->GetXaxis()->SetLabelSize(0.030);
  dummy->GetXaxis()->SetTitle("#sigma_{WZ}^{exp} / #sigma_{WZ}^{theo}");
  dummy->GetYaxis()->SetTitle("");


  // Remove y-axis labels
  //----------------------------------------------------------------------------
  TAxis* yaxis = dummy->GetYaxis();
  
  for (UInt_t j=1; j<yaxis->GetNbins(); j++) yaxis->SetBinLabel(j, "");

  CMS_lumi(canvas, iPeriod, iPos);

  // Save
  //----------------------------------------------------------------------------
  canvas->Update();
  canvas->GetFrame()->DrawClone();
  //  canvas->RedrawAxi();

  //  canvas->SaveAs(Form("%s/ratios%dtev.%s",
  //	      _format.Data(),
  //		      ecm,
  //	      _format.Data()));
  canvas->SaveAs("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis2/WZanalysis/plotovi/crossSectionRatio.pdf");
}


//------------------------------------------------------------------------------
// DrawTLatex
//------------------------------------------------------------------------------
void DrawTLatex(Double_t    x,
		Double_t    y,
		Double_t    tsize,
		Short_t     align,
		const char* text,
		Bool_t      setndc = false)
{
  TLatex* tl = new TLatex(x, y, text);

  tl->SetNDC      (setndc);
  tl->SetTextAlign( align);
  tl->SetTextFont (    42);
  //  tl->SetTextSize ( 15);
  tl->SetTextSize ( tsize);
  tl->Draw("same");  

  
  
}
