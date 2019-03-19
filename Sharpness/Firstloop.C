/*
 * In this code, you correct the difference between two energy bands applying one factor called tau (One simulation)
 */

#include<iostream>
#include<fstream>
#include<istream>
#include<TGraph.h>
#include<TGraph2D.h>
#include<TF1.h>
#include<TFile.h>
#include<TH1D.h>
#include<TStyle.h>
#include<TMath.h>
#include<TRandom.h>
#include<TCanvas.h>
#include<TLegend.h>
#include <TVector.h>
#include <TVectorT.h>
#include<TMath.h>
#include<TKDE.h>
#include<TTree.h>

using namespace std;


const Int_t Np = 12520;   // Number of events (hardcoded)
Double_t t[Np];
Double_t E[Np];

Double_t tmin=0;
Double_t tmax=13000.;

//QG inputs
Double_t z=0.031;
Double_t H0=70.;
Double_t PlanckScale=1.2*TMath::Power(10,19);
Double_t PcTom=3.085*TMath::Power(10,19);
Double_t QGfixedPart=z*PcTom/H0; //(ETA)

void Firstloop(char *fname)
{
	ifstream in;
	in.open(fname);
	Double_t v1,v2;
	Int_t Npoints = 0;
	while(1)
	{
		if (in.eof()) break;
		in >> v1 >> v2;
		t[Npoints] = v1;
		if(v1>tmax)tmax=v1;
		E[Npoints] = v2;
		Npoints++;
	}

  std::cout << "Npoints: " << Npoints-1 << std::endl;
  std::cout << "tmax: " << tmax << std::endl;

  
  //Divide in time for low and high energy (x time, y energy)
  Double_t LE_t[Np], HE_t[Np], LE_E[Np], HE_E[Np];
  Int_t iLE = 0 ;
  Int_t  iHE = 0;
  Double_t cutInE =250.; //Choose the energy for the cut (GeV)
  for (Int_t i=0; i<Np; i++)
    {
      if(E[i]<=cutInE)
      {
    	  LE_t[iLE] = t[i];
    	  LE_E[iLE] = E[i];
    	  iLE++;
      }
      if(E[i]>cutInE)
      {
    	  HE_t[iHE] = t[i];
    	  HE_E[iHE] = E[i];
    	  iHE++;
      }
    }
  std::cout << " LE photons = " << iLE-1 << std::endl;
  std::cout << " HE photons = " << iHE-1 << std::endl;

  
  //Mean energie and timess
  Double_t Sum_E_LE, Sum_E_HE, Sum_t_LE, Sum_t_HE;
  Sum_E_LE =0;
  Sum_E_HE =0;
  Sum_t_LE = 0;
  Sum_t_HE = 0;

  for( Int_t a=0; a<iLE; a++)
    {
	  Sum_E_LE+= LE_E[a];
	  Sum_t_LE+= LE_t[a];
    }
  
  for( Int_t b=0; b<iHE; b++)
    {
	  Sum_E_HE+= HE_E[b];
	  Sum_t_HE+= HE_t[b];
    }

  std::cout << "Mean Energy in LE (TeV) = " << Sum_E_LE/(iLE-1) << std::endl;
  std::cout << "Mean Energy in HE (TeV) = " << Sum_E_HE/(iHE-1) << std::endl;
  std::cout << "Mean time in LE (s) = " << Sum_t_LE/(iLE-1) << std::endl;
  std::cout << "Mean time in HE (s) = " << Sum_t_HE/(iHE-1) << std::endl;


  //Create and fill histograms
  const Int_t bins = 40;
  TH1D* nLE = new TH1D ("nLE", "nLE", bins, tmin, tmax);
  TH1D* nHE = new TH1D ("nHE", "nHE", bins, tmin, tmax);

  for ( Int_t a = 0; a < iLE; a++)
    {
      nLE->Fill(LE_t[a]);
    }

  for ( Int_t b = 0; b < iHE; b++)
    {
      nHE->Fill(HE_t[b]);
    }
  

  //Build CDFs

  TH1D* nLE_CDF = new TH1D ("nLE_CDF", "nLE_CDF", bins, tmin, tmax);
  TH1D* nHE_CDF = new TH1D ("nHE_CDF", "nHE_CDF", bins,tmin, tmax);

  Double_t add_LE = 0;
  Double_t add_HE = 0;
  for(Int_t i = 1; i<=bins; i++)
    {
      Double_t value_LE = nLE->GetBinContent(i);
      Double_t value_HE = nHE->GetBinContent(i);
      nLE_CDF->SetBinContent(i,value_LE + add_LE);
      nHE_CDF->SetBinContent(i,value_HE + add_HE);
      add_LE+= value_LE;
      add_HE+= value_HE;
    }

  /* TCanvas* cosa = new TCanvas("cosa", "Histrograms", 800, 800);
    cosa->Divide(2,3);
    cosa->cd(1);
    nLE->SetXTitle("time(s)");
    nLE->Draw();
    cosa->cd(2);
    nHE->SetXTitle("time(s)");
    nHE->Draw();
    cosa->cd(3);
    nLE_CDF->SetXTitle("time(s)");
    nLE_CDF->Draw();
    cosa->cd(4);
    nHE_CDF->SetXTitle("time(s)");
    nHE_CDF->Draw("same");
    cosa->cd(5);
    //nLE_CDF->SetXTitle("time(s)");
    nLE->Draw();
    nHE->SetLineColor(2);
    nHE->Draw("same");
    cosa->cd(6);
    //nHE_CDF->SetXTitle("time(s)");
    nLE_CDF->Draw();
    nHE_CDF->SetLineColor(2);
    nHE_CDF->Draw("same");*/

    //Comparison of CDFs
    
    Double_t Comp[bins];
    for (Int_t c = 0; c<bins; c++)
      {
    	Comp[c]= TMath::Abs((nLE_CDF->GetBinContent(c+1))-(nHE_CDF->GetBinContent(c+1)));
      }
    Double_t  max = 0.;
    Int_t ind = 0;
    for (Int_t c = 0; c<bins; c++)
      {
    	if (max<=Comp[c])
    	{
    		max = Comp[c];
    		ind = c;
    	}
      }
      std::cout << "The maximum or Dk is: " << max << " located at the bin " << ind << std::endl;
      

      //First correction (For a given value of alpha)
      Double_t alpha = 60.;
      Double_t QGDelay= QGfixedPart*alpha/PlanckScale;
      Double_t tcorr[Np];
      for (Int_t i = 0; i<Np; i++)
      {
    	  tcorr[i]= t[i]-QGDelay*E[i];
      }
      

      //Redivide in energy and time
      Double_t LE_t2[Np], HE_t2[Np], LE_E2[Np], HE_E2[Np];
      Int_t iLE2 = 0 ;
      Int_t  iHE2 = 0;
      for ( Int_t i=0; i<Np; i++)
      {
    	  if(E[i]<=cutInE)
    	  {
    		  LE_t2[iLE2] = tcorr[i];
    		  LE_E2[iLE2] = E[i];
    		  iLE2++;
    	  }
    	  if(E[i]>cutInE)
    	  {
    		  HE_t2[iHE2] = tcorr[i];
    		  HE_E2[iHE2] = E[i];
    		  iHE2++;
    	  }
      }
      std::cout << " Corrected LE photons = " << iLE2-1 << std::endl;
      std::cout << " Corrected HE photons = " << iHE2-1 << std::endl;

      
      //Recompute mean energie and times
      Double_t Sum_E_LE2, Sum_E_HE2, Sum_t_LE2, Sum_t_HE2;
      Sum_E_LE2 =0;
      Sum_E_HE2 =0;
      Sum_t_LE2 = 0;
      Sum_t_HE2 = 0;

      for( Int_t a=0; a<iLE2; a++)
      {
    	  Sum_E_LE2+= LE_E2[a];
    	  Sum_t_LE2+= LE_t2[a];
      }
  
      for( Int_t b=0; b<iHE2; b++)
      {
    	  Sum_E_HE2+= HE_E2[b];
    	  Sum_t_HE2+= HE_t2[b];
      }

      std::cout << "Corrected Mean Energy in LE (TeV) = " << Sum_E_LE2/(iLE2-1) << std::endl;
      std::cout << "Corrected Mean Energy in HE (TeV) = " << Sum_E_HE2/(iHE2-1) << std::endl;
      std::cout << "Corrected Mean time in LE (s) = " << Sum_t_LE2/(iLE2-1) << std::endl;
      std::cout << "Corrected Mean time in HE (s) = " << Sum_t_HE2/(iHE2-1) << std::endl;


    //Recreate and fill histograms
    TH1D* nLE2 = new TH1D ("nLE2", "nLE2", bins, tmin, tmax);
    TH1D* nHE2 = new TH1D ("nHE2", "nHE2", bins, tmin, tmax);
    for (Int_t a = 0; a < iLE2; a++)
    {
      nLE2->Fill(LE_t2[a]);
    }

    for (Int_t b = 0; b < iHE; b++)
    {
      nHE2->Fill(HE_t2[b]);
    }
  


  //Rebuild CDFs
  TH1D* nLE_CDF2 = new TH1D ("nLE_CDF2", "nLE_CDF2", bins, tmin, tmax);
  TH1D* nHE_CDF2 = new TH1D ("nHE_CDF2", "nHE_CDF2", bins, tmin, tmax);
  Double_t add_LE2 = 0;
  Double_t add_HE2 = 0;
  for(Int_t i = 1; i<=bins; i++)
    {
      Double_t value_LE2 = nLE2->GetBinContent(i);
      Double_t value_HE2 = nHE2->GetBinContent(i);
      nLE_CDF2->SetBinContent(i,value_LE2 + add_LE2);
      nHE_CDF2->SetBinContent(i,value_HE2 + add_HE2);
      add_LE2+= value_LE2;
      add_HE2+= value_HE2;
    }

  /*TCanvas* cosa2 = new TCanvas("cosa2", "Histrograms2", 800, 800);
    cosa2->Divide(2,3);
    cosa2->cd(1);
    nLE2->SetXTitle("time(s)");
    nLE2->Draw();
    cosa2->cd(2);
    nHE2->SetXTitle("time(s)");
    nHE2->Draw();
    cosa2->cd(3);
    nLE_CDF2->SetXTitle("time(s)");
    nLE_CDF2->Draw();
    cosa2->cd(4);
    nHE_CDF2->SetXTitle("time(s)");
    nHE_CDF2->Draw("same");
    cosa2->cd(5);
    //nLE_CDF->SetXTitle("time(s)");
    nLE2->Draw();
    nHE2->SetLineColor(2);
    nHE2->Draw("same");
    cosa2->cd(6);
    //nHE_CDF->SetXTitle("time(s)");
    nLE_CDF2->Draw();
    nHE_CDF2->SetLineColor(2);
    nHE_CDF2->Draw("same");*/

    //Comparison of CDFs
    
    Double_t Comp2[bins];
    for ( Int_t c = 0; c<bins; c++)
      {
	Comp2[c]= TMath::Abs((nLE_CDF2->GetBinContent(c+1))-(nHE_CDF2->GetBinContent(c+1)));
      }
    Double_t  max2 = 0.;
    Int_t ind2 = 0;
    for (Int_t c = 0; c<bins; c++)
      {
    	if (max2<=Comp2[c])
    	{
    		max2 = Comp2[c];
    		ind2 = c;
    	}
      }
      std::cout << "The corrected  maximum or Dk is: " << max2 << " located at the corrected bin " << ind2 << std::endl;
	
      TCanvas* cosa3 = new TCanvas("cosa3", "Histrograms3", 800, 800);                                                                      
      cosa3->Divide(2,2);                                                                                                                    
      cosa3->cd(1);
      nLE->SetTitle("Lightcurve before correction");
      nLE->SetXTitle("time(s)");                                                                                                            
      nLE->Draw();
      nHE->SetLineColor(2);
      nHE->SetXTitle("time(s)");                                                                                                             
      nHE->Draw("same");                                                                                                                     
      cosa3->cd(2);
      nLE_CDF->SetTitle("CDF before correction");                                                                                          
      nLE_CDF->SetXTitle("time(s)");                                                                                                         
      nLE_CDF->Draw();
      nHE_CDF->SetLineColor(2);
      nHE_CDF->SetXTitle("time(s)");                                                                                                         
      nHE_CDF->Draw("same");                                                                                                                
      cosa3->cd(3);
      nLE2->SetTitle("Lightcurve after correction");
      nLE2->SetXTitle("time(s)");                                                                                                            
      nLE2->Draw();
      nHE2->SetLineColor(2);
      nHE2->SetXTitle("time(s)");                                                                                                             
      nHE2->Draw("same");
      cosa3->cd(4);
      nLE_CDF2->SetTitle("CDF after correction");                                                                                          
      nLE_CDF2->SetXTitle("time(s)");                                                                                                         
      nLE_CDF2->Draw();
      nHE_CDF2->SetLineColor(2);
      nHE_CDF2->SetXTitle("time(s)");                                                                                                         
      nHE_CDF2->Draw("same"); 

}

