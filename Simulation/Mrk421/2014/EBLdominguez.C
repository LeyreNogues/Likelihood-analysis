#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TStyle.h>

using namespace std;

TGraph *transvse;
TGraph *tauvse;
Double_t transfunc(Double_t *x, Double_t*)  // x = log10(E/GeV) //Computes the opacity
{
  return exp(-tauvse->Eval(x[0]));
}

Double_t transfunc2(Double_t *x, Double_t *par) // x = E/GeV
{
  //  return exp(-tauvse->Eval(log10(x[0]), 0, "S"));
  return par[0]*exp(-tauvse->Eval(log10(x[0])));
}

///////////

void EBLdominguez(Double_t redshift = 0.5)
{
  if (redshift > 1.99)
    {
      cout << "Sorry, this works only for z < 2" << endl;
      return;
    }
        
  Float_t z[39] = {0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895, 0.11684211, 0.13210526, 0.14736842, 0.16263158, 0.17789474, 0.19315789, 0.20842105, 0.22368421, 0.23894737, 0.25421053, 0.26947368, 0.28473684, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.2, 1.4, 1.6, 1.8, 2.};

  ifstream f("tau_dominguez11.dat");
  if (!f.is_open())
    {
      cout << "Cannot open input file" << endl;
      return;
    }

  Char_t ch[1000];
  Double_t data[50];

  TGraph2D* tau_vs_z_vs_e = new TGraph2D;
  Int_t ipoint = 0;


  // Set optical depth = 0 at all redshifts for E=1GeV
  tau_vs_z_vs_e->SetPoint(ipoint++, 0., 0., 0.);
  for (Int_t iz = 0; iz < 39; iz++)
    tau_vs_z_vs_e->SetPoint(ipoint++, 0., z[iz], 0.);

  Int_t i_energy = 1; // Just to count number of E values

  for (;;)  // Loop over energies
    {
      f.getline(ch, 1000);
	  
      if (f.eof())
	break;

      if (ch[0] == '#')
	continue;

      TString str = ch;
      //      cout << str.Data() << endl;

      TObjArray* arr;
      arr = str.Tokenize(" ");

      for (Int_t i = 0; i < arr->GetEntries(); i++)  // Get all optical depths for this energy (Reading one line)
	{
	  data[i] = ((TObjString*)arr->At(i))->GetString().Atof();
	}

      // Fill the optical depth (0) for redshift 0
      tau_vs_z_vs_e->SetPoint(ipoint++, log10(data[0]*1.e3), 0., 0.); // tau=0 for z=0
      
      for (Int_t i = 1; i < arr->GetEntries(); i++)
	{
	  tau_vs_z_vs_e->SetPoint(ipoint++, log10(data[0]*1.e3), z[i-1], data[i]);
	  // E in GeV!!
	  //	  cout << log10(data[0]*1.e3) << " = log10(E/GeV),  z = " << z[i-1] << ",  tau = " << data[i] << endl;
	}
      i_energy++;
    }

  tau_vs_z_vs_e->SetNpx(i_energy);
  tau_vs_z_vs_e->SetNpy(40);
  
  f.close();

  gStyle->SetPalette(1);

  if (! gROOT->FindObject("map") )
    {
      TCanvas* map = new TCanvas;
      map->SetName("map");
      tau_vs_z_vs_e->Draw("zcol");
    }

  transvse = new TGraph();
  tauvse = new TGraph();
  ipoint=0;


  for(Double_t log10e = 0.; log10e < log10(3.e4); log10e += 0.05)
    {
      transvse->SetPoint(ipoint, log10e, exp(-1.*tau_vs_z_vs_e->Interpolate(log10e, redshift)));
      tauvse->SetPoint(ipoint, log10e, tau_vs_z_vs_e->Interpolate(log10e, redshift));
      ipoint++;
    }

  
  if (!gROOT->FindObject("tau"))
    {
      TCanvas* tau = new TCanvas;
      tau->SetName("tau");
      tau->SetGridx();
      tau->SetGridy();
      tau->SetLogy();
      tau->DrawFrame(1., 1.e-4, 4.4, 2.);
      transvse->Draw("l");
    }
  else
    transvse->Draw("l");

  // Transform into a function:
  TF1 *f1 = new TF1("f1", transfunc, 1., 4.5, 0);
  f1->SetLineColor(2);
  f1->Draw("same");

  //  TF1 *f2 = new TF1("f2", transfunc2, 10., 3.e4, 0);
  TF1 *f2 = new TF1("f2", transfunc2, 10., 1.e5, 1);  // NOTE!! : NOT reliable beyond 30 TeV!! (just extrapolated!) Anyway tau too high...
  f2->SetParameter(0, 1.); // Not really needed here, but used as a test for flute.cc
  f2->SetLineColor(2);
  new TCanvas;
  f2->Draw();

  // IMPORTANT!! THIS FUNCTION WILL BE STORED IN NUMERICAL FORM!
  f2->SetNpx(100000); 
  // Attempt below to achieve logarithmic interpolation does not work:
  // TH1D* h2 = (TH1D*) f2->GetHistogram();
  // Double_t xedges[101];
  // xedges[0] = h2->GetXaxis()->GetXmin();
  // xedges[100] = h2->GetXaxis()->GetXmax();
  // for (Int_t i=1; i<100; i++)
  //   xedges[i] = pow(10., log10(xedges[0])+i*(log10(xedges[100])-log10(xedges[0]))/100.);
  // h2->SetBins(100, xedges);
  

  // TEST
//  TFile fout("out.root", "recreate");
//  f1->Write();
//  f2->Write();
//  fout.Close();
  
  return;
}
