void mergeHistograms()
{
  TFile* fData   = new TFile("rootfiles/data.root",        "read");
  TFile* fFakes  = new TFile("rootfiles/data_driven.root", "read");
  TFile* fZgamma = new TFile("rootfiles/Zgamma.root",      "read");
  TFile* fZZ     = new TFile("rootfiles/ZZ.root",          "read");
  TFile* fWZ     = new TFile("rootfiles/WZ.root",          "read");
  TFile* fVVV    = new TFile("rootfiles/VVV.root",         "read");
  TFile* fWV     = new TFile("rootfiles/WV.root",          "read");
  //  TFile* fSyst   = new TFile("rootfiles/syst.root",        "read");

  TH1F* data   = (TH1F*)fData  ->Get("hZmass2_4");
  TH1F* fakes  = (TH1F*)fFakes ->Get("hZmass2_4");
  TH1F* Zgamma = (TH1F*)fZgamma->Get("hZmass2_4");
  TH1F* ZZ     = (TH1F*)fZZ    ->Get("hZmass2_4");
  TH1F* WZ     = (TH1F*)fWZ    ->Get("hZmass2_4");
  TH1F* VVV    = (TH1F*)fVVV   ->Get("hZmass2_4");
  TH1F* WV     = (TH1F*)fWV    ->Get("hZmass2_4");
  //  TH1F* allmc  = (TH1F*)fSyst  ->Get("hZmass2_4");
  TH1F* allmc  = data->Clone();

  data  ->SetNameTitle("data",   "data");
  fakes ->SetNameTitle("fakes",  "fakes");
  Zgamma->SetNameTitle("Zgamma", "Zgamma");
  ZZ    ->SetNameTitle("ZZ",     "ZZ");
  WZ    ->SetNameTitle("WZ",     "WZ");
  VVV   ->SetNameTitle("VVV",    "VVV");
  WV    ->SetNameTitle("WV",     "WV");
  allmc ->SetNameTitle("allmc",  "allmc");

  TFile* output = new TFile("rootfiles/invMass2Lep_8TeV.root", "recreate");

  output->cd();

  data  ->Write();
  fakes ->Write();
  Zgamma->Write();
  ZZ    ->Write();
  WZ    ->Write();
  VVV   ->Write();
  WV    ->Write();
  allmc ->Write();

  output->Close();
}
