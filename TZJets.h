//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 29 22:16:45 2016 by ROOT version 5.32/00
// from TTree latino/probe_tree
// found on file: /STORE/wz8tev/latino_999_TZJets.root
//////////////////////////////////////////////////////////

#ifndef TZJets_h
#define TZJets_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TZJets {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         baseW;
   Float_t         bdt1;
   Float_t         bdt2;
   Float_t         bdt3;
   Float_t         bdt4;
   Float_t         bdt5;
   Float_t         bdt6;
   Float_t         ch1;
   Float_t         ch2;
   Float_t         ch3;
   Float_t         ch4;
   Float_t         chSumEt;
   Float_t         channel;
   Float_t         chmet;
   Float_t         chmetphi;
   Float_t         cjeteta1;
   Float_t         cjeteta2;
   Float_t         cjetid1;
   Float_t         cjetid2;
   Float_t         cjetmva1;
   Float_t         cjetmva2;
   Float_t         cjetphi1;
   Float_t         cjetphi2;
   Float_t         cjetpt1;
   Float_t         cjetpt2;
   Float_t         dataset;
   Float_t         detajj;
   Float_t         dphill;
   Float_t         dphilljet;
   Float_t         dphilljetjet;
   Float_t         dphillmet;
   Float_t         dphilmet;
   Float_t         dphilmet1;
   Float_t         dphilmet2;
   Float_t         drll;
   Float_t         dymva0;
   Float_t         dymva1;
   Float_t         dz1;
   Float_t         dz2;
   Float_t         dz3;
   Float_t         dz4;
   Float_t         effAW;
   Float_t         effBW;
   Float_t         eta1;
   Float_t         eta2;
   Float_t         eta3;
   Float_t         eta4;
   Float_t         fakeAW;
   Float_t         fakeBW;
   Float_t         fakeW;
   Float_t         fermiW;
   Float_t         fourW;
   Float_t         gammaMRStar;
   Float_t         hardbjpb;
   Float_t         hardtche;
   Float_t         higggsGenpt;
   Float_t         hypo;
   Float_t         imet;
   Float_t         ip1;
   Float_t         ip2;
   Float_t         ip3;
   Float_t         ip4;
   Float_t         isSTA1;
   Float_t         isSTA2;
   Float_t         isSTA3;
   Float_t         isSTA4;
   Float_t         isSTA5;
   Float_t         isSTA6;
   Float_t         iso1;
   Float_t         iso2;
   Float_t         iso3;
   Float_t         iso4;
   Float_t         iso5;
   Float_t         iso6;
   Float_t         isomva1;
   Float_t         isomva2;
   Float_t         isomva3;
   Float_t         isomva4;
   Float_t         jetCHM1;
   Float_t         jetCHM2;
   Float_t         jetCHM3;
   Float_t         jetCHM4;
   Float_t         jetCHM5;
   Float_t         jetCHM6;
   Float_t         jetCHM7;
   Float_t         jetGenPartoneta1;
   Float_t         jetGenPartoneta2;
   Float_t         jetGenPartoneta3;
   Float_t         jetGenPartonphi1;
   Float_t         jetGenPartonphi2;
   Float_t         jetGenPartonphi3;
   Float_t         jetGenPartonpid1;
   Float_t         jetGenPartonpid2;
   Float_t         jetGenPartonpid3;
   Float_t         jetGenPartonpt1;
   Float_t         jetGenPartonpt2;
   Float_t         jetGenPartonpt3;
   Float_t         jetGeneta1;
   Float_t         jetGeneta2;
   Float_t         jetGeneta3;
   Float_t         jetGenphi1;
   Float_t         jetGenphi2;
   Float_t         jetGenphi3;
   Float_t         jetGenpt1;
   Float_t         jetGenpt2;
   Float_t         jetGenpt3;
   Float_t         jetNHM1;
   Float_t         jetNHM2;
   Float_t         jetNHM3;
   Float_t         jetNHM4;
   Float_t         jetNHM5;
   Float_t         jetNHM6;
   Float_t         jetNHM7;
   Float_t         jetPhM1;
   Float_t         jetPhM2;
   Float_t         jetPhM3;
   Float_t         jetPhM4;
   Float_t         jetPhM5;
   Float_t         jetPhM6;
   Float_t         jetPhM7;
   Float_t         jetRho;
   Float_t         jetbjpb1;
   Float_t         jetbjpb2;
   Float_t         jetbjpb3;
   Float_t         jetbjpb4;
   Float_t         jetbjpb5;
   Float_t         jetbjpb6;
   Float_t         jetbjpb7;
   Float_t         jetcsv1;
   Float_t         jetcsv2;
   Float_t         jetcsv3;
   Float_t         jetcsv4;
   Float_t         jetcsv5;
   Float_t         jetcsv6;
   Float_t         jetcsv7;
   Float_t         jeteta1;
   Float_t         jeteta2;
   Float_t         jeteta3;
   Float_t         jeteta4;
   Float_t         jeteta5;
   Float_t         jeteta6;
   Float_t         jeteta7;
   Float_t         jetid1;
   Float_t         jetid2;
   Float_t         jetid3;
   Float_t         jetid4;
   Float_t         jetid5;
   Float_t         jetid6;
   Float_t         jetid7;
   Float_t         jetmass1;
   Float_t         jetmass2;
   Float_t         jetmass3;
   Float_t         jetmass4;
   Float_t         jetmass5;
   Float_t         jetmass6;
   Float_t         jetmass7;
   Float_t         jetmva1;
   Float_t         jetmva2;
   Float_t         jetmva3;
   Float_t         jetmva4;
   Float_t         jetmva5;
   Float_t         jetmva6;
   Float_t         jetmva7;
   Float_t         jetphi1;
   Float_t         jetphi2;
   Float_t         jetphi3;
   Float_t         jetphi4;
   Float_t         jetphi5;
   Float_t         jetphi6;
   Float_t         jetphi7;
   Float_t         jetpt1;
   Float_t         jetpt2;
   Float_t         jetpt3;
   Float_t         jetpt4;
   Float_t         jetpt5;
   Float_t         jetpt6;
   Float_t         jetpt7;
   Float_t         jetptd1;
   Float_t         jetptd2;
   Float_t         jetptd3;
   Float_t         jetptd4;
   Float_t         jetptd5;
   Float_t         jetptd6;
   Float_t         jetptd7;
   Float_t         jettche1;
   Float_t         jettche2;
   Float_t         jettche3;
   Float_t         jettche4;
   Float_t         jettche5;
   Float_t         jettche6;
   Float_t         jettche7;
   Float_t         jettchp1;
   Float_t         jettchp2;
   Float_t         jettchp3;
   Float_t         jettchp4;
   Float_t         jettchp5;
   Float_t         jettchp6;
   Float_t         jettchp7;
   Float_t         kfW;
   Float_t         leptonGeneta1;
   Float_t         leptonGeneta2;
   Float_t         leptonGeneta3;
   Float_t         leptonGenphi1;
   Float_t         leptonGenphi2;
   Float_t         leptonGenphi3;
   Float_t         leptonGenpid1;
   Float_t         leptonGenpid2;
   Float_t         leptonGenpid3;
   Float_t         leptonGenpt1;
   Float_t         leptonGenpt2;
   Float_t         leptonGenpt3;
   Float_t         lh1;
   Float_t         lh2;
   Float_t         lh3;
   Float_t         lh4;
   Float_t         mctruth;
   Float_t         metGeneta;
   Float_t         metGenphi;
   Float_t         metGenpt;
   Float_t         mjj;
   Float_t         mll;
   Float_t         mpmet;
   Float_t         mth;
   Float_t         mtw1;
   Float_t         mtw2;
   Float_t         nbjet;
   Float_t         nbjetcsvloose;
   Float_t         nbjetcsvmedium;
   Float_t         nbjetcsvtight;
   Float_t         nbjettche;
   Float_t         nbrem1;
   Float_t         nbrem2;
   Float_t         nbrem3;
   Float_t         nbrem4;
   Float_t         neutrinoGeneta1;
   Float_t         neutrinoGeneta2;
   Float_t         neutrinoGeneta3;
   Float_t         neutrinoGenphi1;
   Float_t         neutrinoGenphi2;
   Float_t         neutrinoGenphi3;
   Float_t         neutrinoGenpid1;
   Float_t         neutrinoGenpid2;
   Float_t         neutrinoGenpid3;
   Float_t         neutrinoGenpt1;
   Float_t         neutrinoGenpt2;
   Float_t         neutrinoGenpt3;
   Float_t         nextra;
   Float_t         njet;
   Float_t         njetid;
   Float_t         njetvbf;
   Float_t         pchmet;
   Float_t         pdfid1;
   Float_t         pdfid2;
   Float_t         pdfscalePDF;
   Float_t         pdfx1;
   Float_t         pdfx1PDF;
   Float_t         pdfx2;
   Float_t         pdfx2PDF;
   Float_t         peaking;
   Float_t         pfSumEt;
   Float_t         pfmet;
   Float_t         pfmetMEtSig;
   Float_t         pfmetSignificance;
   Float_t         pfmetTypeI;
   Float_t         pfmetTypeIphi;
   Float_t         pfmetphi;
   Float_t         phi1;
   Float_t         phi2;
   Float_t         phi3;
   Float_t         phi4;
   Float_t         ppfmet;
   Float_t         predmet;
   Float_t         pt1;
   Float_t         pt2;
   Float_t         pt3;
   Float_t         pt4;
   Float_t         pt5;
   Float_t         pt6;
   Float_t         ptll;
   Float_t         redmet;
   Float_t         sceta1;
   Float_t         sceta2;
   Float_t         sceta3;
   Float_t         sceta4;
   Float_t         softbjpb;
   Float_t         softtche;
   Float_t         tightmu;
   Float_t         triggAW;
   Float_t         triggBW;
   Float_t         trigger;
   Float_t         worstJetLepPt;
   Float_t         yll;
   Float_t         itpu;
   Float_t         nvtx;
   Float_t         ootpum1;
   Float_t         ootpup1;
   Float_t         puAW;
   Float_t         puBW;
   Float_t         trpu;
   Int_t           bveto;
   Int_t           bveto_ip;
   Int_t           bveto_mu;
   Int_t           bveto_munj;
   Int_t           bveto_munj05;
   Int_t           bveto_munj30;
   Int_t           bveto_munj3005;
   Int_t           bveto_nj;
   Int_t           bveto_nj05;
   Int_t           bveto_nj30;
   Int_t           bveto_nj3005;
   Int_t           dphiveto;
   Int_t           pass2012ICHEP1;
   Int_t           pass2012ICHEP2;
   Int_t           pass2012ICHEP3;
   Int_t           pass2012ICHEP4;
   Int_t           pass2012ICHEP5;
   Int_t           pass2012ICHEP6;
   Int_t           pass2012ICHEPLOOSE1;
   Int_t           pass2012ICHEPLOOSE2;
   Int_t           pass2012ICHEPLOOSE3;
   Int_t           pass2012ICHEPLOOSE4;
   Int_t           pass2012ICHEPLOOSE5;
   Int_t           pass2012ICHEPLOOSE6;
   Int_t           sameflav;
   Int_t           zveto;
   UInt_t          run;
   UInt_t          lumi;
   UInt_t          event;
   Float_t         triggW;
   Float_t         effW;
   Float_t         effWUp;
   Float_t         effWDown;
   Float_t         puW;

   // List of branches
   TBranch        *b_baseW;   //!
   TBranch        *b_bdt1;   //!
   TBranch        *b_bdt2;   //!
   TBranch        *b_bdt3;   //!
   TBranch        *b_bdt4;   //!
   TBranch        *b_bdt5;   //!
   TBranch        *b_bdt6;   //!
   TBranch        *b_ch1;   //!
   TBranch        *b_ch2;   //!
   TBranch        *b_ch3;   //!
   TBranch        *b_ch4;   //!
   TBranch        *b_chSumEt;   //!
   TBranch        *b_channel;   //!
   TBranch        *b_chmet;   //!
   TBranch        *b_chmetphi;   //!
   TBranch        *b_cjeteta1;   //!
   TBranch        *b_cjeteta2;   //!
   TBranch        *b_cjetid1;   //!
   TBranch        *b_cjetid2;   //!
   TBranch        *b_cjetmva1;   //!
   TBranch        *b_cjetmva2;   //!
   TBranch        *b_cjetphi1;   //!
   TBranch        *b_cjetphi2;   //!
   TBranch        *b_cjetpt1;   //!
   TBranch        *b_cjetpt2;   //!
   TBranch        *b_dataset;   //!
   TBranch        *b_detajj;   //!
   TBranch        *b_dphill;   //!
   TBranch        *b_dphilljet;   //!
   TBranch        *b_dphilljetjet;   //!
   TBranch        *b_dphillmet;   //!
   TBranch        *b_dphilmet;   //!
   TBranch        *b_dphilmet1;   //!
   TBranch        *b_dphilmet2;   //!
   TBranch        *b_drll;   //!
   TBranch        *b_dymva0;   //!
   TBranch        *b_dymva1;   //!
   TBranch        *b_dz1;   //!
   TBranch        *b_dz2;   //!
   TBranch        *b_dz3;   //!
   TBranch        *b_dz4;   //!
   TBranch        *b_effAW;   //!
   TBranch        *b_effBW;   //!
   TBranch        *b_eta1;   //!
   TBranch        *b_eta2;   //!
   TBranch        *b_eta3;   //!
   TBranch        *b_eta4;   //!
   TBranch        *b_fakeAW;   //!
   TBranch        *b_fakeBW;   //!
   TBranch        *b_fakeW;   //!
   TBranch        *b_fermiW;   //!
   TBranch        *b_fourW;   //!
   TBranch        *b_gammaMRStar;   //!
   TBranch        *b_hardbjpb;   //!
   TBranch        *b_hardtche;   //!
   TBranch        *b_higggsGenpt;   //!
   TBranch        *b_hypo;   //!
   TBranch        *b_imet;   //!
   TBranch        *b_ip1;   //!
   TBranch        *b_ip2;   //!
   TBranch        *b_ip3;   //!
   TBranch        *b_ip4;   //!
   TBranch        *b_isSTA1;   //!
   TBranch        *b_isSTA2;   //!
   TBranch        *b_isSTA3;   //!
   TBranch        *b_isSTA4;   //!
   TBranch        *b_isSTA5;   //!
   TBranch        *b_isSTA6;   //!
   TBranch        *b_iso1;   //!
   TBranch        *b_iso2;   //!
   TBranch        *b_iso3;   //!
   TBranch        *b_iso4;   //!
   TBranch        *b_iso5;   //!
   TBranch        *b_iso6;   //!
   TBranch        *b_isomva1;   //!
   TBranch        *b_isomva2;   //!
   TBranch        *b_isomva3;   //!
   TBranch        *b_isomva4;   //!
   TBranch        *b_jetCHM1;   //!
   TBranch        *b_jetCHM2;   //!
   TBranch        *b_jetCHM3;   //!
   TBranch        *b_jetCHM4;   //!
   TBranch        *b_jetCHM5;   //!
   TBranch        *b_jetCHM6;   //!
   TBranch        *b_jetCHM7;   //!
   TBranch        *b_jetGenPartoneta1;   //!
   TBranch        *b_jetGenPartoneta2;   //!
   TBranch        *b_jetGenPartoneta3;   //!
   TBranch        *b_jetGenPartonphi1;   //!
   TBranch        *b_jetGenPartonphi2;   //!
   TBranch        *b_jetGenPartonphi3;   //!
   TBranch        *b_jetGenPartonpid1;   //!
   TBranch        *b_jetGenPartonpid2;   //!
   TBranch        *b_jetGenPartonpid3;   //!
   TBranch        *b_jetGenPartonpt1;   //!
   TBranch        *b_jetGenPartonpt2;   //!
   TBranch        *b_jetGenPartonpt3;   //!
   TBranch        *b_jetGeneta1;   //!
   TBranch        *b_jetGeneta2;   //!
   TBranch        *b_jetGeneta3;   //!
   TBranch        *b_jetGenphi1;   //!
   TBranch        *b_jetGenphi2;   //!
   TBranch        *b_jetGenphi3;   //!
   TBranch        *b_jetGenpt1;   //!
   TBranch        *b_jetGenpt2;   //!
   TBranch        *b_jetGenpt3;   //!
   TBranch        *b_jetNHM1;   //!
   TBranch        *b_jetNHM2;   //!
   TBranch        *b_jetNHM3;   //!
   TBranch        *b_jetNHM4;   //!
   TBranch        *b_jetNHM5;   //!
   TBranch        *b_jetNHM6;   //!
   TBranch        *b_jetNHM7;   //!
   TBranch        *b_jetPhM1;   //!
   TBranch        *b_jetPhM2;   //!
   TBranch        *b_jetPhM3;   //!
   TBranch        *b_jetPhM4;   //!
   TBranch        *b_jetPhM5;   //!
   TBranch        *b_jetPhM6;   //!
   TBranch        *b_jetPhM7;   //!
   TBranch        *b_jetRho;   //!
   TBranch        *b_jetbjpb1;   //!
   TBranch        *b_jetbjpb2;   //!
   TBranch        *b_jetbjpb3;   //!
   TBranch        *b_jetbjpb4;   //!
   TBranch        *b_jetbjpb5;   //!
   TBranch        *b_jetbjpb6;   //!
   TBranch        *b_jetbjpb7;   //!
   TBranch        *b_jetcsv1;   //!
   TBranch        *b_jetcsv2;   //!
   TBranch        *b_jetcsv3;   //!
   TBranch        *b_jetcsv4;   //!
   TBranch        *b_jetcsv5;   //!
   TBranch        *b_jetcsv6;   //!
   TBranch        *b_jetcsv7;   //!
   TBranch        *b_jeteta1;   //!
   TBranch        *b_jeteta2;   //!
   TBranch        *b_jeteta3;   //!
   TBranch        *b_jeteta4;   //!
   TBranch        *b_jeteta5;   //!
   TBranch        *b_jeteta6;   //!
   TBranch        *b_jeteta7;   //!
   TBranch        *b_jetid1;   //!
   TBranch        *b_jetid2;   //!
   TBranch        *b_jetid3;   //!
   TBranch        *b_jetid4;   //!
   TBranch        *b_jetid5;   //!
   TBranch        *b_jetid6;   //!
   TBranch        *b_jetid7;   //!
   TBranch        *b_jetmass1;   //!
   TBranch        *b_jetmass2;   //!
   TBranch        *b_jetmass3;   //!
   TBranch        *b_jetmass4;   //!
   TBranch        *b_jetmass5;   //!
   TBranch        *b_jetmass6;   //!
   TBranch        *b_jetmass7;   //!
   TBranch        *b_jetmva1;   //!
   TBranch        *b_jetmva2;   //!
   TBranch        *b_jetmva3;   //!
   TBranch        *b_jetmva4;   //!
   TBranch        *b_jetmva5;   //!
   TBranch        *b_jetmva6;   //!
   TBranch        *b_jetmva7;   //!
   TBranch        *b_jetphi1;   //!
   TBranch        *b_jetphi2;   //!
   TBranch        *b_jetphi3;   //!
   TBranch        *b_jetphi4;   //!
   TBranch        *b_jetphi5;   //!
   TBranch        *b_jetphi6;   //!
   TBranch        *b_jetphi7;   //!
   TBranch        *b_jetpt1;   //!
   TBranch        *b_jetpt2;   //!
   TBranch        *b_jetpt3;   //!
   TBranch        *b_jetpt4;   //!
   TBranch        *b_jetpt5;   //!
   TBranch        *b_jetpt6;   //!
   TBranch        *b_jetpt7;   //!
   TBranch        *b_jetptd1;   //!
   TBranch        *b_jetptd2;   //!
   TBranch        *b_jetptd3;   //!
   TBranch        *b_jetptd4;   //!
   TBranch        *b_jetptd5;   //!
   TBranch        *b_jetptd6;   //!
   TBranch        *b_jetptd7;   //!
   TBranch        *b_jettche1;   //!
   TBranch        *b_jettche2;   //!
   TBranch        *b_jettche3;   //!
   TBranch        *b_jettche4;   //!
   TBranch        *b_jettche5;   //!
   TBranch        *b_jettche6;   //!
   TBranch        *b_jettche7;   //!
   TBranch        *b_jettchp1;   //!
   TBranch        *b_jettchp2;   //!
   TBranch        *b_jettchp3;   //!
   TBranch        *b_jettchp4;   //!
   TBranch        *b_jettchp5;   //!
   TBranch        *b_jettchp6;   //!
   TBranch        *b_jettchp7;   //!
   TBranch        *b_kfW;   //!
   TBranch        *b_leptonGeneta1;   //!
   TBranch        *b_leptonGeneta2;   //!
   TBranch        *b_leptonGeneta3;   //!
   TBranch        *b_leptonGenphi1;   //!
   TBranch        *b_leptonGenphi2;   //!
   TBranch        *b_leptonGenphi3;   //!
   TBranch        *b_leptonGenpid1;   //!
   TBranch        *b_leptonGenpid2;   //!
   TBranch        *b_leptonGenpid3;   //!
   TBranch        *b_leptonGenpt1;   //!
   TBranch        *b_leptonGenpt2;   //!
   TBranch        *b_leptonGenpt3;   //!
   TBranch        *b_lh1;   //!
   TBranch        *b_lh2;   //!
   TBranch        *b_lh3;   //!
   TBranch        *b_lh4;   //!
   TBranch        *b_mctruth;   //!
   TBranch        *b_metGeneta;   //!
   TBranch        *b_metGenphi;   //!
   TBranch        *b_metGenpt;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_mll;   //!
   TBranch        *b_mpmet;   //!
   TBranch        *b_mth;   //!
   TBranch        *b_mtw1;   //!
   TBranch        *b_mtw2;   //!
   TBranch        *b_nbjet;   //!
   TBranch        *b_nbjetcsvloose;   //!
   TBranch        *b_nbjetcsvmedium;   //!
   TBranch        *b_nbjetcsvtight;   //!
   TBranch        *b_nbjettche;   //!
   TBranch        *b_nbrem1;   //!
   TBranch        *b_nbrem2;   //!
   TBranch        *b_nbrem3;   //!
   TBranch        *b_nbrem4;   //!
   TBranch        *b_neutrinoGeneta1;   //!
   TBranch        *b_neutrinoGeneta2;   //!
   TBranch        *b_neutrinoGeneta3;   //!
   TBranch        *b_neutrinoGenphi1;   //!
   TBranch        *b_neutrinoGenphi2;   //!
   TBranch        *b_neutrinoGenphi3;   //!
   TBranch        *b_neutrinoGenpid1;   //!
   TBranch        *b_neutrinoGenpid2;   //!
   TBranch        *b_neutrinoGenpid3;   //!
   TBranch        *b_neutrinoGenpt1;   //!
   TBranch        *b_neutrinoGenpt2;   //!
   TBranch        *b_neutrinoGenpt3;   //!
   TBranch        *b_nextra;   //!
   TBranch        *b_njet;   //!
   TBranch        *b_njetid;   //!
   TBranch        *b_njetvbf;   //!
   TBranch        *b_pchmet;   //!
   TBranch        *b_pdfid1;   //!
   TBranch        *b_pdfid2;   //!
   TBranch        *b_pdfscalePDF;   //!
   TBranch        *b_pdfx1;   //!
   TBranch        *b_pdfx1PDF;   //!
   TBranch        *b_pdfx2;   //!
   TBranch        *b_pdfx2PDF;   //!
   TBranch        *b_peaking;   //!
   TBranch        *b_pfSumEt;   //!
   TBranch        *b_pfmet;   //!
   TBranch        *b_pfmetMEtSig;   //!
   TBranch        *b_pfmetSignificance;   //!
   TBranch        *b_pfmetTypeI;   //!
   TBranch        *b_pfmetTypeIphi;   //!
   TBranch        *b_pfmetphi;   //!
   TBranch        *b_phi1;   //!
   TBranch        *b_phi2;   //!
   TBranch        *b_phi3;   //!
   TBranch        *b_phi4;   //!
   TBranch        *b_ppfmet;   //!
   TBranch        *b_predmet;   //!
   TBranch        *b_pt1;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_pt3;   //!
   TBranch        *b_pt4;   //!
   TBranch        *b_pt5;   //!
   TBranch        *b_pt6;   //!
   TBranch        *b_ptll;   //!
   TBranch        *b_redmet;   //!
   TBranch        *b_sceta1;   //!
   TBranch        *b_sceta2;   //!
   TBranch        *b_sceta3;   //!
   TBranch        *b_sceta4;   //!
   TBranch        *b_softbjpb;   //!
   TBranch        *b_softtche;   //!
   TBranch        *b_tightmu;   //!
   TBranch        *b_triggAW;   //!
   TBranch        *b_triggBW;   //!
   TBranch        *b_trigger;   //!
   TBranch        *b_worstJetLepPt;   //!
   TBranch        *b_yll;   //!
   TBranch        *b_itpu;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_ootpum1;   //!
   TBranch        *b_ootpup1;   //!
   TBranch        *b_puAW;   //!
   TBranch        *b_puBW;   //!
   TBranch        *b_trpu;   //!
   TBranch        *b_bveto;   //!
   TBranch        *b_bveto_ip;   //!
   TBranch        *b_bveto_mu;   //!
   TBranch        *b_bveto_munj;   //!
   TBranch        *b_bveto_munj05;   //!
   TBranch        *b_bveto_munj30;   //!
   TBranch        *b_bveto_munj3005;   //!
   TBranch        *b_bveto_nj;   //!
   TBranch        *b_bveto_nj05;   //!
   TBranch        *b_bveto_nj30;   //!
   TBranch        *b_bveto_nj3005;   //!
   TBranch        *b_dphiveto;   //!
   TBranch        *b_pass2012ICHEP1;   //!
   TBranch        *b_pass2012ICHEP2;   //!
   TBranch        *b_pass2012ICHEP3;   //!
   TBranch        *b_pass2012ICHEP4;   //!
   TBranch        *b_pass2012ICHEP5;   //!
   TBranch        *b_pass2012ICHEP6;   //!
   TBranch        *b_pass2012ICHEPLOOSE1;   //!
   TBranch        *b_pass2012ICHEPLOOSE2;   //!
   TBranch        *b_pass2012ICHEPLOOSE3;   //!
   TBranch        *b_pass2012ICHEPLOOSE4;   //!
   TBranch        *b_pass2012ICHEPLOOSE5;   //!
   TBranch        *b_pass2012ICHEPLOOSE6;   //!
   TBranch        *b_sameflav;   //!
   TBranch        *b_zveto;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_triggW;   //!
   TBranch        *b_effW;   //!
   TBranch        *b_effWUp;   //!
   TBranch        *b_effWDown;   //!
   TBranch        *b_puW;   //!

   TZJets(TTree *tree=0);
   virtual ~TZJets();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TZJets_cxx
TZJets::TZJets(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/STORE/wz8tev/latino_999_TZJets.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/STORE/wz8tev/latino_999_TZJets.root");
      }
      f->GetObject("latino",tree);

   }
   Init(tree);
}

TZJets::~TZJets()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TZJets::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TZJets::LoadTree(Long64_t entry)
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

void TZJets::Init(TTree *tree)
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

   fChain->SetBranchAddress("baseW", &baseW, &b_baseW);
   fChain->SetBranchAddress("bdt1", &bdt1, &b_bdt1);
   fChain->SetBranchAddress("bdt2", &bdt2, &b_bdt2);
   fChain->SetBranchAddress("bdt3", &bdt3, &b_bdt3);
   fChain->SetBranchAddress("bdt4", &bdt4, &b_bdt4);
   fChain->SetBranchAddress("bdt5", &bdt5, &b_bdt5);
   fChain->SetBranchAddress("bdt6", &bdt6, &b_bdt6);
   fChain->SetBranchAddress("ch1", &ch1, &b_ch1);
   fChain->SetBranchAddress("ch2", &ch2, &b_ch2);
   fChain->SetBranchAddress("ch3", &ch3, &b_ch3);
   fChain->SetBranchAddress("ch4", &ch4, &b_ch4);
   fChain->SetBranchAddress("chSumEt", &chSumEt, &b_chSumEt);
   fChain->SetBranchAddress("channel", &channel, &b_channel);
   fChain->SetBranchAddress("chmet", &chmet, &b_chmet);
   fChain->SetBranchAddress("chmetphi", &chmetphi, &b_chmetphi);
   fChain->SetBranchAddress("cjeteta1", &cjeteta1, &b_cjeteta1);
   fChain->SetBranchAddress("cjeteta2", &cjeteta2, &b_cjeteta2);
   fChain->SetBranchAddress("cjetid1", &cjetid1, &b_cjetid1);
   fChain->SetBranchAddress("cjetid2", &cjetid2, &b_cjetid2);
   fChain->SetBranchAddress("cjetmva1", &cjetmva1, &b_cjetmva1);
   fChain->SetBranchAddress("cjetmva2", &cjetmva2, &b_cjetmva2);
   fChain->SetBranchAddress("cjetphi1", &cjetphi1, &b_cjetphi1);
   fChain->SetBranchAddress("cjetphi2", &cjetphi2, &b_cjetphi2);
   fChain->SetBranchAddress("cjetpt1", &cjetpt1, &b_cjetpt1);
   fChain->SetBranchAddress("cjetpt2", &cjetpt2, &b_cjetpt2);
   fChain->SetBranchAddress("dataset", &dataset, &b_dataset);
   fChain->SetBranchAddress("detajj", &detajj, &b_detajj);
   fChain->SetBranchAddress("dphill", &dphill, &b_dphill);
   fChain->SetBranchAddress("dphilljet", &dphilljet, &b_dphilljet);
   fChain->SetBranchAddress("dphilljetjet", &dphilljetjet, &b_dphilljetjet);
   fChain->SetBranchAddress("dphillmet", &dphillmet, &b_dphillmet);
   fChain->SetBranchAddress("dphilmet", &dphilmet, &b_dphilmet);
   fChain->SetBranchAddress("dphilmet1", &dphilmet1, &b_dphilmet1);
   fChain->SetBranchAddress("dphilmet2", &dphilmet2, &b_dphilmet2);
   fChain->SetBranchAddress("drll", &drll, &b_drll);
   fChain->SetBranchAddress("dymva0", &dymva0, &b_dymva0);
   fChain->SetBranchAddress("dymva1", &dymva1, &b_dymva1);
   fChain->SetBranchAddress("dz1", &dz1, &b_dz1);
   fChain->SetBranchAddress("dz2", &dz2, &b_dz2);
   fChain->SetBranchAddress("dz3", &dz3, &b_dz3);
   fChain->SetBranchAddress("dz4", &dz4, &b_dz4);
   fChain->SetBranchAddress("effAW", &effAW, &b_effAW);
   fChain->SetBranchAddress("effBW", &effBW, &b_effBW);
   fChain->SetBranchAddress("eta1", &eta1, &b_eta1);
   fChain->SetBranchAddress("eta2", &eta2, &b_eta2);
   fChain->SetBranchAddress("eta3", &eta3, &b_eta3);
   fChain->SetBranchAddress("eta4", &eta4, &b_eta4);
   fChain->SetBranchAddress("fakeAW", &fakeAW, &b_fakeAW);
   fChain->SetBranchAddress("fakeBW", &fakeBW, &b_fakeBW);
   fChain->SetBranchAddress("fakeW", &fakeW, &b_fakeW);
   fChain->SetBranchAddress("fermiW", &fermiW, &b_fermiW);
   fChain->SetBranchAddress("fourW", &fourW, &b_fourW);
   fChain->SetBranchAddress("gammaMRStar", &gammaMRStar, &b_gammaMRStar);
   fChain->SetBranchAddress("hardbjpb", &hardbjpb, &b_hardbjpb);
   fChain->SetBranchAddress("hardtche", &hardtche, &b_hardtche);
   fChain->SetBranchAddress("higggsGenpt", &higggsGenpt, &b_higggsGenpt);
   fChain->SetBranchAddress("hypo", &hypo, &b_hypo);
   fChain->SetBranchAddress("imet", &imet, &b_imet);
   fChain->SetBranchAddress("ip1", &ip1, &b_ip1);
   fChain->SetBranchAddress("ip2", &ip2, &b_ip2);
   fChain->SetBranchAddress("ip3", &ip3, &b_ip3);
   fChain->SetBranchAddress("ip4", &ip4, &b_ip4);
   fChain->SetBranchAddress("isSTA1", &isSTA1, &b_isSTA1);
   fChain->SetBranchAddress("isSTA2", &isSTA2, &b_isSTA2);
   fChain->SetBranchAddress("isSTA3", &isSTA3, &b_isSTA3);
   fChain->SetBranchAddress("isSTA4", &isSTA4, &b_isSTA4);
   fChain->SetBranchAddress("isSTA5", &isSTA5, &b_isSTA5);
   fChain->SetBranchAddress("isSTA6", &isSTA6, &b_isSTA6);
   fChain->SetBranchAddress("iso1", &iso1, &b_iso1);
   fChain->SetBranchAddress("iso2", &iso2, &b_iso2);
   fChain->SetBranchAddress("iso3", &iso3, &b_iso3);
   fChain->SetBranchAddress("iso4", &iso4, &b_iso4);
   fChain->SetBranchAddress("iso5", &iso5, &b_iso5);
   fChain->SetBranchAddress("iso6", &iso6, &b_iso6);
   fChain->SetBranchAddress("isomva1", &isomva1, &b_isomva1);
   fChain->SetBranchAddress("isomva2", &isomva2, &b_isomva2);
   fChain->SetBranchAddress("isomva3", &isomva3, &b_isomva3);
   fChain->SetBranchAddress("isomva4", &isomva4, &b_isomva4);
   fChain->SetBranchAddress("jetCHM1", &jetCHM1, &b_jetCHM1);
   fChain->SetBranchAddress("jetCHM2", &jetCHM2, &b_jetCHM2);
   fChain->SetBranchAddress("jetCHM3", &jetCHM3, &b_jetCHM3);
   fChain->SetBranchAddress("jetCHM4", &jetCHM4, &b_jetCHM4);
   fChain->SetBranchAddress("jetCHM5", &jetCHM5, &b_jetCHM5);
   fChain->SetBranchAddress("jetCHM6", &jetCHM6, &b_jetCHM6);
   fChain->SetBranchAddress("jetCHM7", &jetCHM7, &b_jetCHM7);
   fChain->SetBranchAddress("jetGenPartoneta1", &jetGenPartoneta1, &b_jetGenPartoneta1);
   fChain->SetBranchAddress("jetGenPartoneta2", &jetGenPartoneta2, &b_jetGenPartoneta2);
   fChain->SetBranchAddress("jetGenPartoneta3", &jetGenPartoneta3, &b_jetGenPartoneta3);
   fChain->SetBranchAddress("jetGenPartonphi1", &jetGenPartonphi1, &b_jetGenPartonphi1);
   fChain->SetBranchAddress("jetGenPartonphi2", &jetGenPartonphi2, &b_jetGenPartonphi2);
   fChain->SetBranchAddress("jetGenPartonphi3", &jetGenPartonphi3, &b_jetGenPartonphi3);
   fChain->SetBranchAddress("jetGenPartonpid1", &jetGenPartonpid1, &b_jetGenPartonpid1);
   fChain->SetBranchAddress("jetGenPartonpid2", &jetGenPartonpid2, &b_jetGenPartonpid2);
   fChain->SetBranchAddress("jetGenPartonpid3", &jetGenPartonpid3, &b_jetGenPartonpid3);
   fChain->SetBranchAddress("jetGenPartonpt1", &jetGenPartonpt1, &b_jetGenPartonpt1);
   fChain->SetBranchAddress("jetGenPartonpt2", &jetGenPartonpt2, &b_jetGenPartonpt2);
   fChain->SetBranchAddress("jetGenPartonpt3", &jetGenPartonpt3, &b_jetGenPartonpt3);
   fChain->SetBranchAddress("jetGeneta1", &jetGeneta1, &b_jetGeneta1);
   fChain->SetBranchAddress("jetGeneta2", &jetGeneta2, &b_jetGeneta2);
   fChain->SetBranchAddress("jetGeneta3", &jetGeneta3, &b_jetGeneta3);
   fChain->SetBranchAddress("jetGenphi1", &jetGenphi1, &b_jetGenphi1);
   fChain->SetBranchAddress("jetGenphi2", &jetGenphi2, &b_jetGenphi2);
   fChain->SetBranchAddress("jetGenphi3", &jetGenphi3, &b_jetGenphi3);
   fChain->SetBranchAddress("jetGenpt1", &jetGenpt1, &b_jetGenpt1);
   fChain->SetBranchAddress("jetGenpt2", &jetGenpt2, &b_jetGenpt2);
   fChain->SetBranchAddress("jetGenpt3", &jetGenpt3, &b_jetGenpt3);
   fChain->SetBranchAddress("jetNHM1", &jetNHM1, &b_jetNHM1);
   fChain->SetBranchAddress("jetNHM2", &jetNHM2, &b_jetNHM2);
   fChain->SetBranchAddress("jetNHM3", &jetNHM3, &b_jetNHM3);
   fChain->SetBranchAddress("jetNHM4", &jetNHM4, &b_jetNHM4);
   fChain->SetBranchAddress("jetNHM5", &jetNHM5, &b_jetNHM5);
   fChain->SetBranchAddress("jetNHM6", &jetNHM6, &b_jetNHM6);
   fChain->SetBranchAddress("jetNHM7", &jetNHM7, &b_jetNHM7);
   fChain->SetBranchAddress("jetPhM1", &jetPhM1, &b_jetPhM1);
   fChain->SetBranchAddress("jetPhM2", &jetPhM2, &b_jetPhM2);
   fChain->SetBranchAddress("jetPhM3", &jetPhM3, &b_jetPhM3);
   fChain->SetBranchAddress("jetPhM4", &jetPhM4, &b_jetPhM4);
   fChain->SetBranchAddress("jetPhM5", &jetPhM5, &b_jetPhM5);
   fChain->SetBranchAddress("jetPhM6", &jetPhM6, &b_jetPhM6);
   fChain->SetBranchAddress("jetPhM7", &jetPhM7, &b_jetPhM7);
   fChain->SetBranchAddress("jetRho", &jetRho, &b_jetRho);
   fChain->SetBranchAddress("jetbjpb1", &jetbjpb1, &b_jetbjpb1);
   fChain->SetBranchAddress("jetbjpb2", &jetbjpb2, &b_jetbjpb2);
   fChain->SetBranchAddress("jetbjpb3", &jetbjpb3, &b_jetbjpb3);
   fChain->SetBranchAddress("jetbjpb4", &jetbjpb4, &b_jetbjpb4);
   fChain->SetBranchAddress("jetbjpb5", &jetbjpb5, &b_jetbjpb5);
   fChain->SetBranchAddress("jetbjpb6", &jetbjpb6, &b_jetbjpb6);
   fChain->SetBranchAddress("jetbjpb7", &jetbjpb7, &b_jetbjpb7);
   fChain->SetBranchAddress("jetcsv1", &jetcsv1, &b_jetcsv1);
   fChain->SetBranchAddress("jetcsv2", &jetcsv2, &b_jetcsv2);
   fChain->SetBranchAddress("jetcsv3", &jetcsv3, &b_jetcsv3);
   fChain->SetBranchAddress("jetcsv4", &jetcsv4, &b_jetcsv4);
   fChain->SetBranchAddress("jetcsv5", &jetcsv5, &b_jetcsv5);
   fChain->SetBranchAddress("jetcsv6", &jetcsv6, &b_jetcsv6);
   fChain->SetBranchAddress("jetcsv7", &jetcsv7, &b_jetcsv7);
   fChain->SetBranchAddress("jeteta1", &jeteta1, &b_jeteta1);
   fChain->SetBranchAddress("jeteta2", &jeteta2, &b_jeteta2);
   fChain->SetBranchAddress("jeteta3", &jeteta3, &b_jeteta3);
   fChain->SetBranchAddress("jeteta4", &jeteta4, &b_jeteta4);
   fChain->SetBranchAddress("jeteta5", &jeteta5, &b_jeteta5);
   fChain->SetBranchAddress("jeteta6", &jeteta6, &b_jeteta6);
   fChain->SetBranchAddress("jeteta7", &jeteta7, &b_jeteta7);
   fChain->SetBranchAddress("jetid1", &jetid1, &b_jetid1);
   fChain->SetBranchAddress("jetid2", &jetid2, &b_jetid2);
   fChain->SetBranchAddress("jetid3", &jetid3, &b_jetid3);
   fChain->SetBranchAddress("jetid4", &jetid4, &b_jetid4);
   fChain->SetBranchAddress("jetid5", &jetid5, &b_jetid5);
   fChain->SetBranchAddress("jetid6", &jetid6, &b_jetid6);
   fChain->SetBranchAddress("jetid7", &jetid7, &b_jetid7);
   fChain->SetBranchAddress("jetmass1", &jetmass1, &b_jetmass1);
   fChain->SetBranchAddress("jetmass2", &jetmass2, &b_jetmass2);
   fChain->SetBranchAddress("jetmass3", &jetmass3, &b_jetmass3);
   fChain->SetBranchAddress("jetmass4", &jetmass4, &b_jetmass4);
   fChain->SetBranchAddress("jetmass5", &jetmass5, &b_jetmass5);
   fChain->SetBranchAddress("jetmass6", &jetmass6, &b_jetmass6);
   fChain->SetBranchAddress("jetmass7", &jetmass7, &b_jetmass7);
   fChain->SetBranchAddress("jetmva1", &jetmva1, &b_jetmva1);
   fChain->SetBranchAddress("jetmva2", &jetmva2, &b_jetmva2);
   fChain->SetBranchAddress("jetmva3", &jetmva3, &b_jetmva3);
   fChain->SetBranchAddress("jetmva4", &jetmva4, &b_jetmva4);
   fChain->SetBranchAddress("jetmva5", &jetmva5, &b_jetmva5);
   fChain->SetBranchAddress("jetmva6", &jetmva6, &b_jetmva6);
   fChain->SetBranchAddress("jetmva7", &jetmva7, &b_jetmva7);
   fChain->SetBranchAddress("jetphi1", &jetphi1, &b_jetphi1);
   fChain->SetBranchAddress("jetphi2", &jetphi2, &b_jetphi2);
   fChain->SetBranchAddress("jetphi3", &jetphi3, &b_jetphi3);
   fChain->SetBranchAddress("jetphi4", &jetphi4, &b_jetphi4);
   fChain->SetBranchAddress("jetphi5", &jetphi5, &b_jetphi5);
   fChain->SetBranchAddress("jetphi6", &jetphi6, &b_jetphi6);
   fChain->SetBranchAddress("jetphi7", &jetphi7, &b_jetphi7);
   fChain->SetBranchAddress("jetpt1", &jetpt1, &b_jetpt1);
   fChain->SetBranchAddress("jetpt2", &jetpt2, &b_jetpt2);
   fChain->SetBranchAddress("jetpt3", &jetpt3, &b_jetpt3);
   fChain->SetBranchAddress("jetpt4", &jetpt4, &b_jetpt4);
   fChain->SetBranchAddress("jetpt5", &jetpt5, &b_jetpt5);
   fChain->SetBranchAddress("jetpt6", &jetpt6, &b_jetpt6);
   fChain->SetBranchAddress("jetpt7", &jetpt7, &b_jetpt7);
   fChain->SetBranchAddress("jetptd1", &jetptd1, &b_jetptd1);
   fChain->SetBranchAddress("jetptd2", &jetptd2, &b_jetptd2);
   fChain->SetBranchAddress("jetptd3", &jetptd3, &b_jetptd3);
   fChain->SetBranchAddress("jetptd4", &jetptd4, &b_jetptd4);
   fChain->SetBranchAddress("jetptd5", &jetptd5, &b_jetptd5);
   fChain->SetBranchAddress("jetptd6", &jetptd6, &b_jetptd6);
   fChain->SetBranchAddress("jetptd7", &jetptd7, &b_jetptd7);
   fChain->SetBranchAddress("jettche1", &jettche1, &b_jettche1);
   fChain->SetBranchAddress("jettche2", &jettche2, &b_jettche2);
   fChain->SetBranchAddress("jettche3", &jettche3, &b_jettche3);
   fChain->SetBranchAddress("jettche4", &jettche4, &b_jettche4);
   fChain->SetBranchAddress("jettche5", &jettche5, &b_jettche5);
   fChain->SetBranchAddress("jettche6", &jettche6, &b_jettche6);
   fChain->SetBranchAddress("jettche7", &jettche7, &b_jettche7);
   fChain->SetBranchAddress("jettchp1", &jettchp1, &b_jettchp1);
   fChain->SetBranchAddress("jettchp2", &jettchp2, &b_jettchp2);
   fChain->SetBranchAddress("jettchp3", &jettchp3, &b_jettchp3);
   fChain->SetBranchAddress("jettchp4", &jettchp4, &b_jettchp4);
   fChain->SetBranchAddress("jettchp5", &jettchp5, &b_jettchp5);
   fChain->SetBranchAddress("jettchp6", &jettchp6, &b_jettchp6);
   fChain->SetBranchAddress("jettchp7", &jettchp7, &b_jettchp7);
   fChain->SetBranchAddress("kfW", &kfW, &b_kfW);
   fChain->SetBranchAddress("leptonGeneta1", &leptonGeneta1, &b_leptonGeneta1);
   fChain->SetBranchAddress("leptonGeneta2", &leptonGeneta2, &b_leptonGeneta2);
   fChain->SetBranchAddress("leptonGeneta3", &leptonGeneta3, &b_leptonGeneta3);
   fChain->SetBranchAddress("leptonGenphi1", &leptonGenphi1, &b_leptonGenphi1);
   fChain->SetBranchAddress("leptonGenphi2", &leptonGenphi2, &b_leptonGenphi2);
   fChain->SetBranchAddress("leptonGenphi3", &leptonGenphi3, &b_leptonGenphi3);
   fChain->SetBranchAddress("leptonGenpid1", &leptonGenpid1, &b_leptonGenpid1);
   fChain->SetBranchAddress("leptonGenpid2", &leptonGenpid2, &b_leptonGenpid2);
   fChain->SetBranchAddress("leptonGenpid3", &leptonGenpid3, &b_leptonGenpid3);
   fChain->SetBranchAddress("leptonGenpt1", &leptonGenpt1, &b_leptonGenpt1);
   fChain->SetBranchAddress("leptonGenpt2", &leptonGenpt2, &b_leptonGenpt2);
   fChain->SetBranchAddress("leptonGenpt3", &leptonGenpt3, &b_leptonGenpt3);
   fChain->SetBranchAddress("lh1", &lh1, &b_lh1);
   fChain->SetBranchAddress("lh2", &lh2, &b_lh2);
   fChain->SetBranchAddress("lh3", &lh3, &b_lh3);
   fChain->SetBranchAddress("lh4", &lh4, &b_lh4);
   fChain->SetBranchAddress("mctruth", &mctruth, &b_mctruth);
   fChain->SetBranchAddress("metGeneta", &metGeneta, &b_metGeneta);
   fChain->SetBranchAddress("metGenphi", &metGenphi, &b_metGenphi);
   fChain->SetBranchAddress("metGenpt", &metGenpt, &b_metGenpt);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("mll", &mll, &b_mll);
   fChain->SetBranchAddress("mpmet", &mpmet, &b_mpmet);
   fChain->SetBranchAddress("mth", &mth, &b_mth);
   fChain->SetBranchAddress("mtw1", &mtw1, &b_mtw1);
   fChain->SetBranchAddress("mtw2", &mtw2, &b_mtw2);
   fChain->SetBranchAddress("nbjet", &nbjet, &b_nbjet);
   fChain->SetBranchAddress("nbjetcsvloose", &nbjetcsvloose, &b_nbjetcsvloose);
   fChain->SetBranchAddress("nbjetcsvmedium", &nbjetcsvmedium, &b_nbjetcsvmedium);
   fChain->SetBranchAddress("nbjetcsvtight", &nbjetcsvtight, &b_nbjetcsvtight);
   fChain->SetBranchAddress("nbjettche", &nbjettche, &b_nbjettche);
   fChain->SetBranchAddress("nbrem1", &nbrem1, &b_nbrem1);
   fChain->SetBranchAddress("nbrem2", &nbrem2, &b_nbrem2);
   fChain->SetBranchAddress("nbrem3", &nbrem3, &b_nbrem3);
   fChain->SetBranchAddress("nbrem4", &nbrem4, &b_nbrem4);
   fChain->SetBranchAddress("neutrinoGeneta1", &neutrinoGeneta1, &b_neutrinoGeneta1);
   fChain->SetBranchAddress("neutrinoGeneta2", &neutrinoGeneta2, &b_neutrinoGeneta2);
   fChain->SetBranchAddress("neutrinoGeneta3", &neutrinoGeneta3, &b_neutrinoGeneta3);
   fChain->SetBranchAddress("neutrinoGenphi1", &neutrinoGenphi1, &b_neutrinoGenphi1);
   fChain->SetBranchAddress("neutrinoGenphi2", &neutrinoGenphi2, &b_neutrinoGenphi2);
   fChain->SetBranchAddress("neutrinoGenphi3", &neutrinoGenphi3, &b_neutrinoGenphi3);
   fChain->SetBranchAddress("neutrinoGenpid1", &neutrinoGenpid1, &b_neutrinoGenpid1);
   fChain->SetBranchAddress("neutrinoGenpid2", &neutrinoGenpid2, &b_neutrinoGenpid2);
   fChain->SetBranchAddress("neutrinoGenpid3", &neutrinoGenpid3, &b_neutrinoGenpid3);
   fChain->SetBranchAddress("neutrinoGenpt1", &neutrinoGenpt1, &b_neutrinoGenpt1);
   fChain->SetBranchAddress("neutrinoGenpt2", &neutrinoGenpt2, &b_neutrinoGenpt2);
   fChain->SetBranchAddress("neutrinoGenpt3", &neutrinoGenpt3, &b_neutrinoGenpt3);
   fChain->SetBranchAddress("nextra", &nextra, &b_nextra);
   fChain->SetBranchAddress("njet", &njet, &b_njet);
   fChain->SetBranchAddress("njetid", &njetid, &b_njetid);
   fChain->SetBranchAddress("njetvbf", &njetvbf, &b_njetvbf);
   fChain->SetBranchAddress("pchmet", &pchmet, &b_pchmet);
   fChain->SetBranchAddress("pdfid1", &pdfid1, &b_pdfid1);
   fChain->SetBranchAddress("pdfid2", &pdfid2, &b_pdfid2);
   fChain->SetBranchAddress("pdfscalePDF", &pdfscalePDF, &b_pdfscalePDF);
   fChain->SetBranchAddress("pdfx1", &pdfx1, &b_pdfx1);
   fChain->SetBranchAddress("pdfx1PDF", &pdfx1PDF, &b_pdfx1PDF);
   fChain->SetBranchAddress("pdfx2", &pdfx2, &b_pdfx2);
   fChain->SetBranchAddress("pdfx2PDF", &pdfx2PDF, &b_pdfx2PDF);
   fChain->SetBranchAddress("peaking", &peaking, &b_peaking);
   fChain->SetBranchAddress("pfSumEt", &pfSumEt, &b_pfSumEt);
   fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
   fChain->SetBranchAddress("pfmetMEtSig", &pfmetMEtSig, &b_pfmetMEtSig);
   fChain->SetBranchAddress("pfmetSignificance", &pfmetSignificance, &b_pfmetSignificance);
   fChain->SetBranchAddress("pfmetTypeI", &pfmetTypeI, &b_pfmetTypeI);
   fChain->SetBranchAddress("pfmetTypeIphi", &pfmetTypeIphi, &b_pfmetTypeIphi);
   fChain->SetBranchAddress("pfmetphi", &pfmetphi, &b_pfmetphi);
   fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
   fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
   fChain->SetBranchAddress("phi3", &phi3, &b_phi3);
   fChain->SetBranchAddress("phi4", &phi4, &b_phi4);
   fChain->SetBranchAddress("ppfmet", &ppfmet, &b_ppfmet);
   fChain->SetBranchAddress("predmet", &predmet, &b_predmet);
   fChain->SetBranchAddress("pt1", &pt1, &b_pt1);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("pt3", &pt3, &b_pt3);
   fChain->SetBranchAddress("pt4", &pt4, &b_pt4);
   fChain->SetBranchAddress("pt5", &pt5, &b_pt5);
   fChain->SetBranchAddress("pt6", &pt6, &b_pt6);
   fChain->SetBranchAddress("ptll", &ptll, &b_ptll);
   fChain->SetBranchAddress("redmet", &redmet, &b_redmet);
   fChain->SetBranchAddress("sceta1", &sceta1, &b_sceta1);
   fChain->SetBranchAddress("sceta2", &sceta2, &b_sceta2);
   fChain->SetBranchAddress("sceta3", &sceta3, &b_sceta3);
   fChain->SetBranchAddress("sceta4", &sceta4, &b_sceta4);
   fChain->SetBranchAddress("softbjpb", &softbjpb, &b_softbjpb);
   fChain->SetBranchAddress("softtche", &softtche, &b_softtche);
   fChain->SetBranchAddress("tightmu", &tightmu, &b_tightmu);
   fChain->SetBranchAddress("triggAW", &triggAW, &b_triggAW);
   fChain->SetBranchAddress("triggBW", &triggBW, &b_triggBW);
   fChain->SetBranchAddress("trigger", &trigger, &b_trigger);
   fChain->SetBranchAddress("worstJetLepPt", &worstJetLepPt, &b_worstJetLepPt);
   fChain->SetBranchAddress("yll", &yll, &b_yll);
   fChain->SetBranchAddress("itpu", &itpu, &b_itpu);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("ootpum1", &ootpum1, &b_ootpum1);
   fChain->SetBranchAddress("ootpup1", &ootpup1, &b_ootpup1);
   fChain->SetBranchAddress("puAW", &puAW, &b_puAW);
   fChain->SetBranchAddress("puBW", &puBW, &b_puBW);
   fChain->SetBranchAddress("trpu", &trpu, &b_trpu);
   fChain->SetBranchAddress("bveto", &bveto, &b_bveto);
   fChain->SetBranchAddress("bveto_ip", &bveto_ip, &b_bveto_ip);
   fChain->SetBranchAddress("bveto_mu", &bveto_mu, &b_bveto_mu);
   fChain->SetBranchAddress("bveto_munj", &bveto_munj, &b_bveto_munj);
   fChain->SetBranchAddress("bveto_munj05", &bveto_munj05, &b_bveto_munj05);
   fChain->SetBranchAddress("bveto_munj30", &bveto_munj30, &b_bveto_munj30);
   fChain->SetBranchAddress("bveto_munj3005", &bveto_munj3005, &b_bveto_munj3005);
   fChain->SetBranchAddress("bveto_nj", &bveto_nj, &b_bveto_nj);
   fChain->SetBranchAddress("bveto_nj05", &bveto_nj05, &b_bveto_nj05);
   fChain->SetBranchAddress("bveto_nj30", &bveto_nj30, &b_bveto_nj30);
   fChain->SetBranchAddress("bveto_nj3005", &bveto_nj3005, &b_bveto_nj3005);
   fChain->SetBranchAddress("dphiveto", &dphiveto, &b_dphiveto);
   fChain->SetBranchAddress("pass2012ICHEP1", &pass2012ICHEP1, &b_pass2012ICHEP1);
   fChain->SetBranchAddress("pass2012ICHEP2", &pass2012ICHEP2, &b_pass2012ICHEP2);
   fChain->SetBranchAddress("pass2012ICHEP3", &pass2012ICHEP3, &b_pass2012ICHEP3);
   fChain->SetBranchAddress("pass2012ICHEP4", &pass2012ICHEP4, &b_pass2012ICHEP4);
   fChain->SetBranchAddress("pass2012ICHEP5", &pass2012ICHEP5, &b_pass2012ICHEP5);
   fChain->SetBranchAddress("pass2012ICHEP6", &pass2012ICHEP6, &b_pass2012ICHEP6);
   fChain->SetBranchAddress("pass2012ICHEPLOOSE1", &pass2012ICHEPLOOSE1, &b_pass2012ICHEPLOOSE1);
   fChain->SetBranchAddress("pass2012ICHEPLOOSE2", &pass2012ICHEPLOOSE2, &b_pass2012ICHEPLOOSE2);
   fChain->SetBranchAddress("pass2012ICHEPLOOSE3", &pass2012ICHEPLOOSE3, &b_pass2012ICHEPLOOSE3);
   fChain->SetBranchAddress("pass2012ICHEPLOOSE4", &pass2012ICHEPLOOSE4, &b_pass2012ICHEPLOOSE4);
   fChain->SetBranchAddress("pass2012ICHEPLOOSE5", &pass2012ICHEPLOOSE5, &b_pass2012ICHEPLOOSE5);
   fChain->SetBranchAddress("pass2012ICHEPLOOSE6", &pass2012ICHEPLOOSE6, &b_pass2012ICHEPLOOSE6);
   fChain->SetBranchAddress("sameflav", &sameflav, &b_sameflav);
   fChain->SetBranchAddress("zveto", &zveto, &b_zveto);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("triggW", &triggW, &b_triggW);
   fChain->SetBranchAddress("effW", &effW, &b_effW);
   fChain->SetBranchAddress("effWUp", &effWUp, &b_effWUp);
   fChain->SetBranchAddress("effWDown", &effWDown, &b_effWDown);
   fChain->SetBranchAddress("puW", &puW, &b_puW);
   Notify();
}

Bool_t TZJets::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TZJets::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TZJets::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TZJets_cxx
