//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 11 23:00:24 2015 by ROOT version 5.32/00
// from TTree MET/0
// found on file: wzMetSystematics.root
//////////////////////////////////////////////////////////

#ifndef metsys_h
#define metsys_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class metsys {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Float_t         patPFMet;
   Float_t         patPFMetPhi;
   Float_t         patType1CorrectedPFMet;
   Float_t         patType1CorrectedPFMetPhi;
   Float_t         patType1CorrectedPFMetElectronEnDown;
   Float_t         patType1CorrectedPFMetElectronEnDownPhi;
   Float_t         patType1CorrectedPFMetElectronEnUp;
   Float_t         patType1CorrectedPFMetElectronEnUpPhi;
   Float_t         patType1CorrectedPFMetJetEnDown;
   Float_t         patType1CorrectedPFMetJetEnDownPhi;
   Float_t         patType1CorrectedPFMetJetEnUp;
   Float_t         patType1CorrectedPFMetJetEnUpPhi;
   Float_t         patType1CorrectedPFMetJetResDown;
   Float_t         patType1CorrectedPFMetJetResDownPhi;
   Float_t         patType1CorrectedPFMetJetResUp;
   Float_t         patType1CorrectedPFMetJetResUpPhi;
   Float_t         patType1CorrectedPFMetMuonEnDown;
   Float_t         patType1CorrectedPFMetMuonEnDownPhi;
   Float_t         patType1CorrectedPFMetMuonEnUp;
   Float_t         patType1CorrectedPFMetMuonEnUpPhi;
   Float_t         patType1CorrectedPFMetTauEnDown;
   Float_t         patType1CorrectedPFMetTauEnDownPhi;
   Float_t         patType1CorrectedPFMetTauEnUp;
   Float_t         patType1CorrectedPFMetTauEnUpPhi;
   Float_t         patType1CorrectedPFMetUnclusteredEnDown;
   Float_t         patType1CorrectedPFMetUnclusteredEnDownPhi;
   Float_t         patType1CorrectedPFMetUnclusteredEnUp;
   Float_t         patType1CorrectedPFMetUnclusteredEnUpPhi;
   Float_t         Phi;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_patPFMet;   //!
   TBranch        *b_patPFMetPhi;   //!
   TBranch        *b_patType1CorrectedPFMet;   //!
   TBranch        *b_patType1CorrectedPFMetPhi;   //!
   TBranch        *b_patType1CorrectedPFMetElectronEnDown;   //!
   TBranch        *b_patType1CorrectedPFMetElectronEnDownPhi;   //!
   TBranch        *b_patType1CorrectedPFMetElectronEnUp;   //!
   TBranch        *b_patType1CorrectedPFMetElectronEnUpPhi;   //!
   TBranch        *b_patType1CorrectedPFMetJetEnDown;   //!
   TBranch        *b_patType1CorrectedPFMetJetEnDownPhi;   //!
   TBranch        *b_patType1CorrectedPFMetJetEnUp;   //!
   TBranch        *b_patType1CorrectedPFMetJetEnUpPhi;   //!
   TBranch        *b_patType1CorrectedPFMetJetResDown;   //!
   TBranch        *b_patType1CorrectedPFMetJetResDownPhi;   //!
   TBranch        *b_patType1CorrectedPFMetJetResUp;   //!
   TBranch        *b_patType1CorrectedPFMetJetResUpPhi;   //!
   TBranch        *b_patType1CorrectedPFMetMuonEnDown;   //!
   TBranch        *b_patType1CorrectedPFMetMuonEnDownPhi;   //!
   TBranch        *b_patType1CorrectedPFMetMuonEnUp;   //!
   TBranch        *b_patType1CorrectedPFMetMuonEnUpPhi;   //!
   TBranch        *b_patType1CorrectedPFMetTauEnDown;   //!
   TBranch        *b_patType1CorrectedPFMetTauEnDownPhi;   //!
   TBranch        *b_patType1CorrectedPFMetTauEnUp;   //!
   TBranch        *b_patType1CorrectedPFMetTauEnUpPhi;   //!
   TBranch        *b_patType1CorrectedPFMetUnclusteredEnDown;   //!
   TBranch        *b_patType1CorrectedPFMetUnclusteredEnDownPhi;   //!
   TBranch        *b_patType1CorrectedPFMetUnclusteredEnUp;   //!
   TBranch        *b_patType1CorrectedPFMetUnclusteredEnUpPhi;   //!
   TBranch        *b_Phi;   //!

   metsys(TTree *tree=0);
   virtual ~metsys();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef metsys_cxx
metsys::metsys(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("wzMetSystematics.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("wzMetSystematics.root");
      }
      f->GetObject("MET",tree);

   }
   Init(tree);
}

metsys::~metsys()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t metsys::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t metsys::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void metsys::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("patPFMet", &patPFMet, &b_patPFMet);
   fChain->SetBranchAddress("patPFMetPhi", &patPFMetPhi, &b_patPFMetPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMet", &patType1CorrectedPFMet, &b_patType1CorrectedPFMet);
   fChain->SetBranchAddress("patType1CorrectedPFMetPhi", &patType1CorrectedPFMetPhi, &b_patType1CorrectedPFMetPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetElectronEnDown", &patType1CorrectedPFMetElectronEnDown, &b_patType1CorrectedPFMetElectronEnDown);
   fChain->SetBranchAddress("patType1CorrectedPFMetElectronEnDownPhi", &patType1CorrectedPFMetElectronEnDownPhi, &b_patType1CorrectedPFMetElectronEnDownPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetElectronEnUp", &patType1CorrectedPFMetElectronEnUp, &b_patType1CorrectedPFMetElectronEnUp);
   fChain->SetBranchAddress("patType1CorrectedPFMetElectronEnUpPhi", &patType1CorrectedPFMetElectronEnUpPhi, &b_patType1CorrectedPFMetElectronEnUpPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetEnDown", &patType1CorrectedPFMetJetEnDown, &b_patType1CorrectedPFMetJetEnDown);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetEnDownPhi", &patType1CorrectedPFMetJetEnDownPhi, &b_patType1CorrectedPFMetJetEnDownPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetEnUp", &patType1CorrectedPFMetJetEnUp, &b_patType1CorrectedPFMetJetEnUp);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetEnUpPhi", &patType1CorrectedPFMetJetEnUpPhi, &b_patType1CorrectedPFMetJetEnUpPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetResDown", &patType1CorrectedPFMetJetResDown, &b_patType1CorrectedPFMetJetResDown);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetResDownPhi", &patType1CorrectedPFMetJetResDownPhi, &b_patType1CorrectedPFMetJetResDownPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetResUp", &patType1CorrectedPFMetJetResUp, &b_patType1CorrectedPFMetJetResUp);
   fChain->SetBranchAddress("patType1CorrectedPFMetJetResUpPhi", &patType1CorrectedPFMetJetResUpPhi, &b_patType1CorrectedPFMetJetResUpPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetMuonEnDown", &patType1CorrectedPFMetMuonEnDown, &b_patType1CorrectedPFMetMuonEnDown);
   fChain->SetBranchAddress("patType1CorrectedPFMetMuonEnDownPhi", &patType1CorrectedPFMetMuonEnDownPhi, &b_patType1CorrectedPFMetMuonEnDownPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetMuonEnUp", &patType1CorrectedPFMetMuonEnUp, &b_patType1CorrectedPFMetMuonEnUp);
   fChain->SetBranchAddress("patType1CorrectedPFMetMuonEnUpPhi", &patType1CorrectedPFMetMuonEnUpPhi, &b_patType1CorrectedPFMetMuonEnUpPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetTauEnDown", &patType1CorrectedPFMetTauEnDown, &b_patType1CorrectedPFMetTauEnDown);
   fChain->SetBranchAddress("patType1CorrectedPFMetTauEnDownPhi", &patType1CorrectedPFMetTauEnDownPhi, &b_patType1CorrectedPFMetTauEnDownPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetTauEnUp", &patType1CorrectedPFMetTauEnUp, &b_patType1CorrectedPFMetTauEnUp);
   fChain->SetBranchAddress("patType1CorrectedPFMetTauEnUpPhi", &patType1CorrectedPFMetTauEnUpPhi, &b_patType1CorrectedPFMetTauEnUpPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetUnclusteredEnDown", &patType1CorrectedPFMetUnclusteredEnDown, &b_patType1CorrectedPFMetUnclusteredEnDown);
   fChain->SetBranchAddress("patType1CorrectedPFMetUnclusteredEnDownPhi", &patType1CorrectedPFMetUnclusteredEnDownPhi, &b_patType1CorrectedPFMetUnclusteredEnDownPhi);
   fChain->SetBranchAddress("patType1CorrectedPFMetUnclusteredEnUp", &patType1CorrectedPFMetUnclusteredEnUp, &b_patType1CorrectedPFMetUnclusteredEnUp);
   fChain->SetBranchAddress("patType1CorrectedPFMetUnclusteredEnUpPhi", &patType1CorrectedPFMetUnclusteredEnUpPhi, &b_patType1CorrectedPFMetUnclusteredEnUpPhi);
   fChain->SetBranchAddress("Phi", &Phi, &b_Phi);
   Notify();
}

Bool_t metsys::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void metsys::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t metsys::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef metsys_cxx
