/*
 * In this code, you correct the difference between two energy bands applying differect factors called d (One simulation)
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
Int_t QGorder = 1;
Double_t z=0.031;
Double_t H0=70.;
Double_t PlanckScale=1.2*TMath::Power(10,19);
Double_t PcTom=3.085*TMath::Power(10,19);
Double_t QGfixedPart=z*PcTom/H0; //(ETA)

void Loop(char *fname)
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

  std::cout << "t[0]: " << t[0] << std::endl;
  std::cout << "t[1]: " << t[1] << std::endl;
  std::cout << "t[2]: " << t[2] << std::endl;
  std::cout << "t[Npoints-4]: " << t[Npoints-4] << std::endl;
  std::cout << "t[Npoints-3]: " << t[Npoints-3] << std::endl;
  std::cout << "t[Npoints-2]: " << t[Npoints-2] << std::endl;
  std::cout << "t[Npoints-1]: " << t[Npoints-1] << std::endl;
  std::cout << "E[0]: " << E[0] << std::endl;
  std::cout << "E[1]: " << E[1] << std::endl;
  std::cout << "E[2]: " << E[2] << std::endl;
  std::cout << "E[Npoints-4]: " << E[Npoints-4] << std::endl;
  std::cout << "E[Npoints-3]: " << E[Npoints-3] << std::endl;
  std::cout << "E[Npoints-2]: " << E[Npoints-2] << std::endl;
  std::cout << "E[Npoints-1]: " << E[Npoints-1] << std::endl;

      
   //++++++++++++++++++++++++++++++++++++++++++  LOOP FOR ALPHA  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
  //Variables used in all the iterations
  Double_t alpha_min = 0;
  Double_t alpha_max = 0;
  Double_t alpha_step = 0;
  if(QGorder ==1){
	   alpha_min = -2000;
	   alpha_max = 2000;
	   alpha_step = 1.;
  }
  if(QGorder ==2){
	   alpha_min = -4;
	   alpha_max = 4;
	   alpha_step = 0.02;
  }
  Double_t alpha = alpha_min;
  Double_t Dk[6000] = {0};//Big number
  Double_t alphas[6000] = {0};
  Double_t tcorr[Np] = {0};
  Double_t LE_t[Np] = {0}, HE_t[Np]= {0}, LE_E[Np]= {0}, HE_E[Np]= {0};
  Double_t Sum_E_LE = 0, Sum_E_HE = 0, Sum_t_LE = 0, Sum_t_HE = 0;
  const Int_t bins = 200;
  Double_t Comp[bins] = {0};
  Int_t loop = 0;


  while(alpha<=alpha_max)
    {
//      std::cout << "alpha: " << alpha << " in loop number " << loop << std::endl;
      // if(alpha == 0.0)
      //{
      //The histograms are created and deleted in every iteration
      TH1D* nLE = new TH1D ("LE events", "nLE", bins, tmin, tmax);
      TH1D* nHE = new TH1D ("HE events", "nHE", bins, tmin, tmax);
      TH1D* nLE_CDF = new TH1D ("nLE_CDF", "nLE_CDF", bins, tmin, tmax);
      TH1D* nHE_CDF = new TH1D ("nHE_CDF", "nHE_CDF", bins, tmin, tmax);
	

      //Apply correction
      Double_t QGDelay= QGfixedPart*alpha/PlanckScale;
      Double_t tcorr[Np];
      for (Int_t i = 0; i<Np; i++)
      {
    	  if(QGorder == 1)tcorr[i]= t[i]-QGDelay*E[i];
    	  if(QGorder == 2)tcorr[i]= t[i]-QGDelay*E[i]*E[i];
      }
      
      //Divide in two energy bands
      Int_t iLE = 0 ;
      Int_t  iHE = 0;
      Double_t cutInE = 250.; //Choose the energy for the cut (GeV)
      for ( Int_t i=0; i<Np; i++)
	  {
	    if(E[i]<=cutInE)
	    {
	    	LE_t[iLE] = tcorr[i];
	    	LE_E[iLE] = E[i];
	    	iLE++;
	    }
	    if(E[i]>cutInE)
	    {
	    	HE_t[iHE] = tcorr[i];
	    	HE_E[iHE] = E[i];
	    	iHE++;
	    }
	  }
	
	
//      std::cout << " LE photons = " << iLE-1 << std::endl;
//      std::cout << " HE photons = " << iHE-1 << std::endl;
	

      //Compute mean energie and times
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
//      std::cout << "Mean Energy in LE (TeV) = " << Sum_E_LE/(iLE-1) << std::endl;
//      std::cout << "Mean Energy in HE (TeV) = " << Sum_E_HE/(iHE-1) << std::endl;
//      std::cout << "Difference in energy (TeV) = " << (Sum_E_HE/(iHE-1))-(Sum_E_LE/(iLE-1)) << std::endl;
//       std::cout << "Mean time in LE (s) = " << Sum_t_LE/(iLE-1) << std::endl;
//      std::cout << "Mean time in HE (s) = " << Sum_t_HE/(iHE-1) << std::endl;
	


      //Fill histograms
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


      //Build CDFs
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

	  
//      TCanvas* cosa_loop = new TCanvas("cosa_loop", "Histrograms_loop", 800, 800);
//      cosa_loop->Divide(2,3);
//      cosa_loop->cd(1);
//	  nLE->SetXTitle("time(s)");
//	  nLE->Draw();
//	  cosa_loop->cd(2);
//	  nHE->SetXTitle("time(s)");
//	  nHE->Draw();
//	  cosa_loop->cd(3);
//	  nLE_CDF->SetXTitle("time(s)");
//	  nLE_CDF->Draw();
//	  cosa_loop->cd(4);
//	  nHE_CDF->SetXTitle("time(s)");
//	  nHE_CDF->Draw("same");
//	  cosa_loop->cd(5);
//	  //nLE_CDF->SetXTitle("time(s)");
//	  nLE->Draw();
//	  nHE->SetLineColor(2);
//	  nHE->Draw("same");
//	  cosa_loop->cd(6);
//	  //nHE_CDF->SetXTitle("time(s)");
//	  nLE_CDF->Draw();
//	  nHE_CDF->SetLineColor(2);
//	  nHE_CDF->Draw("same");
	  


      //Comparison of CDFs
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
//	  std::cout << "The corrected  maximum or Dk is: " << max << " located at the corrected bin " << ind << std::endl;
	  Dk[loop]= max;
	  alphas[loop]= alpha;

//	  std::cout << "Dk: " << Dk[loop] << std::endl;

	  delete nLE;
	  delete nHE;
	  delete nLE_CDF;
	  delete nHE_CDF;

	  // }//This is put for tests

	  alpha = alpha+alpha_step;
	  loop++;
    }


  //++++++++++++++++++++++++++++++++++++++++++  END OF LOOP  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //Look for the value of alpha giving a minimum Dk
  Int_t index=0;
  Double_t bestalpha = 0, bestmax = 0;
  bestmax = 1000;
  for (Int_t y=0; y<loop; y++)
  {
	  cout << "For alpha: " << alphas[y] << " Dk is: " << Dk[y] << endl;
	  if (Dk[y] <= bestmax)
      {
		  bestmax = Dk[y];
    	  bestalpha = alphas[y];
      }
  }

  TGraph* Dkvsalpha = new TGraph(loop, alphas, Dk);
//  TCanvas* result = new TCanvas();
  Dkvsalpha->SetTitle("");
  Dkvsalpha->SetMarkerStyle(33);
  Dkvsalpha->GetXaxis()->SetTitle("alpha");
  Dkvsalpha->GetYaxis()->SetTitle("Kolmogorov distance");
//  Dkvsalpha->Draw("");


  std::cout << "FINAL RESULT: Best alpha is " << bestalpha << " with a Dk: " << bestmax << std::endl;
  TTree tree("t1","Tree with simple variables");
  Double_t final_Dk = 0;
  tree.Branch("alpha",&alpha,"px/D");
  tree.Branch("Dk",&final_Dk,"px/D");
  cout << "Tree created" << endl;
  alpha = bestalpha;
  final_Dk = bestmax;

  tree.Fill();
  //Create de TFile
  TFile* Likelihood_outputfile = new TFile("Sharpness_KazumaRandom_Linear_DATA.root","recreate");
  tree.Write();
  Dkvsalpha->Write();

  cout << "Root file created"  << endl;


}

