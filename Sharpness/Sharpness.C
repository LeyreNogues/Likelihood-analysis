/*
 * In this code, you just check the difference between two energy bands (One simulation)
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
#include<TStyle.h>

using namespace std;

const Int_t Np = 12520;   // Number of events (hardcoded)
Double_t t[Np];
Double_t E[Np];
Int_t ind = 0;

Double_t tmin=0;
Double_t tmax=13000.;

void Sharpness(char *fname)
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


  //Divide in time for low and high energy
  Double_t LE_t[Np], HE_t[Np], LE_E[Np], HE_E[Np];
  Int_t iLE = 0 ;
  Int_t  iHE = 0;
  Double_t cutInE =250.; //Choose the energy for the cut (GeV)


  for ( Int_t i=0; i<Np; i++)
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


  //Mean energies and times
  Double_t Sum_E_LE, Sum_E_HE, Sum_t_LE, Sum_t_HE;
  Sum_E_LE =0;
  Sum_E_HE =0;
  Sum_t_LE = 0;
  Sum_t_HE = 0;
  for(Int_t a=0; a<iLE; a++)
    {
	  Sum_E_LE+= LE_E[a];
	  Sum_t_LE+= LE_t[a];
    }
  for(Int_t b=0; b<iHE; b++)
    {
	  Sum_E_HE+= HE_E[b];
	  Sum_t_HE+= HE_t[b];
    }
  std::cout << "Mean Energy in LE (TeV) = " << Sum_E_LE/(iLE-1) << std::endl;
  std::cout << "Mean Energy in HE (TeV) = " << Sum_E_HE/(iHE-1) << std::endl;
  std::cout << "Difference in energy (TeV) = " << (Sum_E_HE/(iHE-1))-(Sum_E_LE/(iLE-1)) << std::endl;
   std::cout << "Mean time in LE (s) = " << Sum_t_LE/(iLE-1) << std::endl;
  std::cout << "Mean time in HE (s) = " << Sum_t_HE/(iHE-1) << std::endl;


  //Create and fill histograms
  const Int_t bins = 40;
  TH1D* nLE = new TH1D ("LE events", "LE events", bins, tmin, tmax);
  TH1D* nHE = new TH1D ("HE events", "HE events", bins, tmin, tmax);


  for ( Int_t a = 0; a < iLE; a++)
    {
      nLE->Fill(LE_t[a]);
    }

  for ( Int_t b = 0; b < iHE; b++)
    {
      nHE->Fill(HE_t[b]);
    }
  
  nLE->Scale(1./nLE->Integral());
  nHE->Scale(1./nHE->Integral());

  TH1D* nHE_copy = (TH1D*)nHE->Clone("name");


  //Build CDFs
  TH1D* nLE_CDF = new TH1D ("LE events CDF", "LE events CDF", bins, tmin, tmax);
  TH1D* nHE_CDF = new TH1D ("HE events CDF", "HE events CDF", bins, tmin, tmax);


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

  TH1D* nHE_CDF_copy = (TH1D*)nHE_CDF->Clone("name2");

  TCanvas* cosa = new TCanvas("cosa", "Histograms", 800, 800);
  cosa->Divide(2,3);
  cosa->cd(1);
  nLE->SetStats(0);
  nLE->SetXTitle("time(s)");
  nLE->Draw();
  cosa->cd(2);
  nHE->SetXTitle("time(s)");
  nHE->SetLineColor(2);
  nHE->SetStats(0);
  nHE->Draw();
  cosa->cd(3);
  nLE_CDF->SetXTitle("time(s)");
  nLE_CDF->SetStats(0);
  nLE_CDF->Draw();
  cosa->cd(4);
  nHE_CDF->SetXTitle("time(s)");
  nHE_CDF->SetLineColor(2);
  nHE_CDF->SetStats(0);
  nHE_CDF->Draw("same");
  cosa->cd(5);
  nHE_copy->SetXTitle("time(s)");
  nHE_copy->Draw("");
  nHE_copy->SetStats(0);
  nHE_copy->SetTitle("");
  nHE_copy->SetLineColor(2);
  nLE->Draw("same");
  cosa->cd(6);
  nHE_CDF_copy->SetXTitle("time(s)");
  nHE_CDF_copy->SetTitle("");
  nHE_CDF_copy->SetStats(0);
  nHE_CDF_copy->Draw("");
  nHE_CDF_copy->SetLineColor(2);
  nLE_CDF->Draw("same");


  //Comparison of CDFs
  Double_t Comp[bins];
  for(Int_t c = 0; c<bins; c++)
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
}

