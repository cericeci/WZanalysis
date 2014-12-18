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
  
  if (gotVarName) {
    variable = variableName;
  } else {
    variable = "LeadJetPt";
  }


  bool printBLUEmatrix(true);


  TFile * finput= TFile::Open("/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/unfoldingFinalResults/systematics.root");
  std::ostringstream outfilename;
  outfilename<<"/users/ltikvica/CMSSW_4_2_9_HLT1/src/latinosAnalysis/rezultati/rootFiles/unfaoldingFinalResults/combination_"<<variable<<".root";
  TFile * fout= new TFile(outfilename.str().c_str(), "RECREATE");  

  if (gotHistoBinning) {

    UnfoldingHistogramFactory * histoFac = UnfoldingHistogramFactory::GetInstance();
    histoFac->SetBinning(binningFileName);

  }

  const int nChannels(4);
  TH1D * h_qcdScale[nChannels];
  TH1D * h_PDFsys[nChannels];
  TH1D * h_leptTrgEff[nChannels];
  TH1D * h_Etsys[nChannels];
  TH1D * h_muMomScale[nChannels];
  TH1D * h_elEnScale[nChannels];
  TH1D * h_pileupSys[nChannels];
  TH1D * h_ZZxs[nChannels];
  TH1D * h_Zgammaxs[nChannels];
  TH1D * h_dataDrivensys[nChannels];
  TH1D * h_bckgSys[nChannels];
  TH1D * h_JESsys[nChannels];
  TH1D * h_JERsys[nChannels];
  TH1D * h_crossSection[nChannels];

  TH1D * h_crossSection_combination;
  TH1D * h_combStat;
  TH1D * h_combSyst;

  //ovo treba deklarirat!!!!!!!!!
  if (variable=="Njets"){
    h_crossSection_combination= UnfoldingHistogramFactory::createNjetsHistogram("h_xs_comb", "h_xs_comb");
    h_combStat = UnfoldingHistogramFactory::createNjetsHistogram("h_combStat", "h_combStat"); 
    h_combSyst = UnfoldingHistogramFactory::createNjetsHistogram("h_combSyst", "h_combSyst");
  }
  else if (variable=="LeadingJetPt"){
    h_crossSection_combination= UnfoldingHistogramFactory::createLeadingJetHistogram("h_xs_comb", "h_xs_comb");
    h_combStat = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combStat", "h_combStat"); 
    h_combSyst = UnfoldingHistogramFactory::createLeadingJetHistogram("h_combSyst", "h_combSyst");
  }
  else if (variable=="Zpt"){
    h_crossSection_combination= UnfoldingHistogramFactory::createZPtHistogram("h_xs_comb", "h_xs_comb");
    h_combStat = UnfoldingHistogramFactory::createZPtHistogram("h_combStat", "h_combStat"); 
    h_combSyst = UnfoldingHistogramFactory::createZPtHistogram("h_combSyst", "h_combSyst");
  }
  else
    std::cout<<"UNKNOWN VARIABLE!!!"<<std::endl;


  for (int hist=0; hist<nChannels; hist++){
    std::ostringstream qcdScaleName,PDFsysName,leptTrgEffName,EtsysName, 
      muMomScaleName, elEnScaleName, pileupSysName, ZZxsName, ZgammaxsName, 
      dataDrivenName, bckgSysName, xsName, JESsysName, JERsysName;
    qcdScaleName<<"h_qcdScale_"<<hist;
    PDFsysName<<"h_PDFsys_"<<hist;
    leptTrgEffName<<"h_leptTrgEff_"<<hist;
    EtsysName<<"h_Etsys_"<<hist;
    muMomScaleName<<"h_muMomScale_"<<hist;
    elEnScaleName<<"h_elEnScale_"<<hist;
    pileupSysName<<"h_pileupSys_"<<hist;
    ZZxsName<<"h_ZZxs_"<<hist;
    ZgammaxsName<<"h_Zgammaxs_"<<hist;
    dataDrivenName<<"h_dataDrivensys_"<<hist;
    bckgSysName<<"h_bckgSys_"<<hist;
    xsName<<"h_crossSection_"<<hist;
    JESsysName<<"h_JESsys_"<<hist;
    JERsysName<<"h_JERsys_"<<hist;

    h_qcdScale[hist] = (TH1D*) (finput ->Get(qcdScaleName.str().c_str())->Clone(qcdScaleName.str().c_str()));
    h_PDFsys[hist] = (TH1D*) (finput ->Get(PDFsysName.str().c_str())->Clone(PDFsysName.str().c_str()));
    h_leptTrgEff[hist] = (TH1D*) (finput ->Get(leptTrgEffName.str().c_str())->Clone(leptTrgEffName.str().c_str()));
    h_Etsys[hist] = (TH1D*) (finput ->Get(EtsysName.str().c_str())->Clone(EtsysName.str().c_str()));
    h_muMomScale[hist] = (TH1D*) (finput ->Get(muMomScaleName.str().c_str())->Clone(muMomScaleName.str().c_str()));
    h_elEnScale[hist] = (TH1D*) (finput ->Get(elEnScaleName.str().c_str())->Clone(elEnScaleName.str().c_str()));
    h_pileupSys[hist] = (TH1D*) (finput ->Get(pileupSysName.str().c_str())->Clone(pileupSysName.str().c_str()));
    h_ZZxs[hist] = (TH1D*) (finput ->Get(ZZxsName.str().c_str())->Clone(ZZxsName.str().c_str()));
    h_Zgammaxs[hist] = (TH1D*) (finput ->Get(ZgammaxsName.str().c_str())->Clone(ZgammaxsName.str().c_str()));
    h_dataDrivensys[hist] = (TH1D*) (finput ->Get(dataDrivenName.str().c_str())->Clone(dataDrivenName.str().c_str()));
    h_bckgSys[hist] = (TH1D*) (finput ->Get(bckgSysName.str().c_str())->Clone(bckgSysName.str().c_str()));
    h_crossSection[hist] = (TH1D*) (finput ->Get(xsName.str().c_str())->Clone(xsName.str().c_str()));
    
    if (variable!="Zpt"){
      std::cout<<variable<<"  JET SYSTEMATICS :D"<<std::endl;
      h_JESsys[hist] = (TH1D*) (finput->Get(JESsysName.str().c_str())->Clone(JESsysName.str().c_str()));
      h_JERsys[hist] = (TH1D*) (finput->Get(JERsysName.str().c_str())->Clone(JERsysName.str().c_str()));
    }
    
  }

  //looop over each bin 
  
  for (int bin=0; bin< (h_crossSection[0]->GetNbinsX()+1); bin++){
    double elements[16];
    for (int el=0; el<16; el++) elements[el]=0;

    double statisticError[nChannels];
    double systematicError2[nChannels];
    double systematicError[nChannels];
    //statistic and systematic errors:
    for (int nCh=0; nCh<nChannels; nCh++){
      statisticError[nCh]=h_crossSection[nCh]->GetBinError(bin);
      systematicError2[nCh]=(pow(h_qcdScale[nCh]->GetBinContent(bin),2)+
			     pow(h_PDFsys[nCh]->GetBinContent(bin), 2)+
			     pow(h_leptTrgEff[nCh]->GetBinContent(bin),2)+
			     pow(h_Etsys[nCh]->GetBinContent(bin), 2)+
			     pow(h_muMomScale[nCh]->GetBinContent(bin),2)+
			     pow(h_elEnScale[nCh]->GetBinContent(bin),2)+
			     pow(h_pileupSys[nCh]->GetBinContent(bin),2)+
			     pow(h_ZZxs[nCh]->GetBinContent(bin),2)+
			     pow(h_Zgammaxs[nCh]->GetBinContent(bin),2)+
			     pow(h_dataDrivensys[nCh]->GetBinContent(bin),2)+
			     pow(h_bckgSys[nCh]->GetBinContent(bin),2));
      if (variable!="Zpt"){
	systematicError2[nCh]+= pow (h_JESsys[nCh]->GetBinContent(bin),2)+
	  pow (h_JERsys[nCh]->GetBinContent(bin),2);
      }
      systematicError[nCh]= sqrt(systematicError2[nCh]);
    }
    //common elements
    double commonSys[4][4];
    
    for (int cha=0; cha<4; cha++){
      for (int chb=0; chb<4; chb++){
	commonSys[cha][chb]=(h_Etsys[cha]->GetBinContent(bin))*(h_Etsys[chb]->GetBinContent(bin))
	  + (h_pileupSys[cha]->GetBinContent(bin))*(h_pileupSys[chb]->GetBinContent(bin))
	  + (h_PDFsys[cha]->GetBinContent(bin))*(h_PDFsys[chb]->GetBinContent(bin))
	  + (h_qcdScale[cha]->GetBinContent(bin))*(h_qcdScale[chb]->GetBinContent(bin))
	  + (h_bckgSys[cha]->GetBinContent(bin))*(h_bckgSys[chb]->GetBinContent(bin));
	
	if (variable!="Zpt"){
	  commonSys[cha][chb]+=(h_JESsys[cha]->GetBinContent(bin))*(h_JESsys[chb]->GetBinContent(bin));
	}
      }
    }
    
    //diagonal elements
    
    elements[0]=pow(systematicError[0],2) + pow(statisticError[0],2);
    elements[5]=pow(systematicError[1],2) + pow(statisticError[1],2);
    elements[10]=pow(systematicError[2],2) + pow(statisticError[2],2);
    elements[15]=pow(systematicError[3],2) + pow(statisticError[3],2);
    
    //matrix is symetric
    //channels 0 and 1
    elements[4]= elements[1]= (h_crossSection[0]->GetBinContent(bin))*(h_crossSection[1]->GetBinContent(bin))*
      (commonSys[0][1] + (h_elEnScale[0]->GetBinContent(bin))*(h_elEnScale[1]->GetBinContent(bin))+
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[1]->GetBinContent(bin))+
       (h_leptTrgEff[0]->GetBinContent(bin))*(sqrt(2/3)*(h_leptTrgEff[1]->GetBinContent(bin))));
    
    
    //channels 0 and 2
    elements[8]=elements[2]=(h_crossSection[0]->GetBinContent(bin))*(h_crossSection[2]->GetBinContent(bin))*
      (commonSys[0][2] + (h_elEnScale[0]->GetBinContent(bin))*(h_elEnScale[2]->GetBinContent(bin)) +
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[2]->GetBinContent(bin))+
       (h_leptTrgEff[0]->GetBinContent(bin))*(sqrt(1/3)*(h_leptTrgEff[2]->GetBinContent(2))));
    
    //channels 0 and 3
    elements[12]=elements[3]=(h_crossSection[0]->GetBinContent(bin))*(h_crossSection[3]->GetBinContent(bin))*
      (commonSys[0][3]+ (h_elEnScale[0]->GetBinContent(bin))*(h_elEnScale[3]->GetBinContent(bin))+
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[3]->GetBinContent(bin)));
    
    //channels 1 and 2
    elements[9]=elements[6]= (h_crossSection[1]->GetBinContent(bin))*(h_crossSection[2]->GetBinContent(bin))*
      (commonSys[1][2]+ (h_elEnScale[1]->GetBinContent(bin))*(h_elEnScale[2]->GetBinContent(bin))
       + (h_muMomScale[1]->GetBinContent(bin))*(h_muMomScale[2]->GetBinContent(bin)) 
       + 2*sqrt(1/3)*(h_leptTrgEff[1]->GetBinContent(bin))*sqrt(2/3)*(h_leptTrgEff[2]->GetBinContent(bin)));
    
    //channels 1 and 3
    elements[13]=elements[7]= (h_crossSection[1]->GetBinContent(bin))*(h_crossSection[3]->GetBinContent(bin))*
      (commonSys[1][3]+ (h_elEnScale[1]->GetBinContent(bin))*(h_elEnScale[3]->GetBinContent(bin)) 
       + (h_muMomScale[1]->GetBinContent(bin))*(h_muMomScale[3]->GetBinContent(bin))
       + (h_leptTrgEff[1]->GetBinContent(bin))* (sqrt(1/3)* (h_leptTrgEff[3]->GetBinContent(bin))));
    
    //channels 2 and 3
    elements[14]=elements[11]= (h_crossSection[2]->GetBinContent(bin))*(h_crossSection[3]->GetBinContent(bin))*
      (commonSys[2][3]+ (h_elEnScale[2]->GetBinContent(bin))*(h_elEnScale[3]->GetBinContent(bin)) 
       + (h_muMomScale[2]->GetBinContent(bin))*(h_muMomScale[3]->GetBinContent(bin))
       + (h_leptTrgEff[2]->GetBinContent(bin))* (sqrt(2/3)* (h_leptTrgEff[3]->GetBinContent(bin))));
    
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
    double final_Xsec = alphaCH[0]*(h_crossSection[0]->GetBinContent(bin)) + alphaCH[1]*(h_crossSection[1]->GetBinContent(bin)) 
      + alphaCH[2]*(h_crossSection[2]->GetBinContent(bin)) + alphaCH[3]*(h_crossSection[3]->GetBinContent(bin));
    
    Double_t combined_error=0;
    Double_t *copyRef = errMatCopy.GetMatrixArray();
    for (int ier=0;ier<4;ier++)
      for (int jer=0;jer<4;jer++) {
	combined_error+=alphaCH[ier]*alphaCH[jer]*copyRef[ier*4+jer];
      }
    Double_t stat_err_tot2 =  pow(alphaCH[0]*statisticError[0],2) + pow(alphaCH[1]*statisticError[1],2)
      +pow(alphaCH[2]*statisticError[2],2) + pow(alphaCH[3]*statisticError[3],2);
    
    Double_t stat_err_tot = sqrt(stat_err_tot2);
    Double_t syst_err_tot = sqrt (combined_error - stat_err_tot2);
    
    combined_error=sqrt(combined_error);
    
    if(printBLUEmatrix){
      std::cout << "combined sigma(WZ) = " << final_Xsec << " +- " << stat_err_tot << "(stat.) +- " << syst_err_tot << "(syst) +- "
		<< 0.026*final_Xsec << " (lumi) pb " << endl;
      std::cout << endl;
    }
    double total_error=sqrt(stat_err_tot*stat_err_tot + syst_err_tot*syst_err_tot);
    h_crossSection_combination->SetBinContent(bin, final_Xsec);
    h_crossSection_combination->SetBinError(bin, final_Xsec);
    h_combStat->SetBinContent(bin, stat_err_tot);
    h_combSyst->SetBinContent(bin, syst_err_tot);
  }
  
  fout->cd();
  h_crossSection_combination->Write();
  h_combStat->Write();
  h_combSyst->Write();
 
}
