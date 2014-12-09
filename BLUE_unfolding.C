#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

void BLUE_unfolding()
{
  //read histograms form .root files
  //definirati file-ove iz kojih citam i u koje pisem
  const int nChannels(4);
  TH1D * h_qcdScale[nChannels];
  TH1D * h_PDFsys[nChannels];
  TH1D * h_leptTrgEff[nChannels];
  TH1D * h_Etsys[nChannels];
  TH1D * h_muMoScale[nChannels];
  TH1D * h_elEnScale[nChannels];
  TH1D * h_pileupSys[nChannels];
  TH1D * h_ZZxs[nChannels];
  TH1D * h_Zgammaxs[nChannels];
  TH1D * h_dataDrivensys[nChannels];
  TH1D * h_bckgSys[nChannels];
  
  TH1D * h_crossSection[nChannels];
  
  //results:
  //treba ih kreirati da imaju isti binning kao i originalni (staviti neke opcije npr. pa napraviti taj binning)
  TH1D * h_crossSection_combination;
  TH1D * h_combStat;
  TH1D * h_combSyst;

  //looop over each bin 
  
  for (int bin=0; bin< (h_crossSection[0]->GetNbinsX()+1); bin++){
    double elements[16];
    for (int el=0; el<16; el++) elements[el]=0;

    //statistic and systematic errors:
    for (int nCh=0; nCh<4; nCh++){
      statisticError[nCh]=h_crossSection[nCh]->GetBinError(bin);
      systematicError[nCh]=sqrt(pow(h_qcdScale[nCh]->GetBinContent(bin),2)+
				pow(h_PDFsys[nCh]->GetBinContent(bin), 2)+
				pow(h_leptTrgEff[nCh]->GetBinContent(bin),2)+
				pow(h_Etsys[nCh]->GetBinContent(bin), 2)+
				pow(h_muMoScale[nCh]->GetBinContent(bin),2)+
				pow(h_elEnScale[nCh]->GetbinContent(bin),2)+
				pow(h_pileupSys[nCh]->GetBinContent(bin),2)+
				pow(h_ZZxs[nCh]->GetBinContent(bin),2)+
				pow(h_Zgammaxs[nCh]->GetBinContent(bin),2)+
				pow(h_dataDrivensys[nCh]->GetBinContent(bin),2)+
				pow(h_bckgSys[nCh]->GetBinContent(bin),2));
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
      (commonSys[0][2] + (h_elEnScale[0]->GetBinContent(bin))*(h+elEnScalse[2]->GetBinContent(bin)) +
       (h_muMomScale[0]->GetBinContent(bin))*(h_muMomScale[2]->GetBinContent(bin))+
       (h_leptTrgEff[0]->GetBincontent(bin))*(sqrt(1/3)*(h_leptTrgEff[2]->GetBinContent(2))));
    
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
