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

const Int_t Np = 1491;   // Number of events (hardcoded)
Double_t t[Np];
Double_t E[Np];

Double_t tmin=0;
Double_t tmax=2750.;

//QG inputs
Double_t z=0.034;
Double_t H0=70.;
Double_t PlanckScale=1.2*TMath::Power(10,19);
Double_t PcTom=3.085*TMath::Power(10,19);
Double_t QGfixedPart=z*PcTom/H0; //(ETA)

void Flareloop()
{
  //++++++++++++++++++++++++++++++++++++++++++++++  LOOP FOR FLARES  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  const Int_t numOfFlares = 100;
  Double_t results[numOfFlares];
  Int_t Npoints = 0; 
  Double_t alpha;
  Double_t Dk[numOfFlares];
  Double_t alphas[numOfFlares];
  Double_t tcorr[Np];
  Double_t LE_t[Np], HE_t[Np], LE_E[Np], HE_E[Np];
  Double_t Sum_E_LE, Sum_E_HE, Sum_t_LE, Sum_t_HE;
  const Int_t bins = 40;
  Double_t Comp[bins];
  Int_t loop;


  for(Int_t l=0; l<numOfFlares; l++)
  {
      char *fname;
      ifstream in;
      Double_t v1,v2;
//      if(l==0){                                                                                                                                             //Read file and fill the arrays t and E
      std::cout << "antes " << t[0] << " " << E[0] << std::endl;
      memset(t, 0, Np);
      memset(E, 0, Np);
      std::cout << "despues " << t[0] << " " << E[0] << std::endl;
      Npoints = 0;
      fname = Form("/home/lnogues/workspace_cpp/LIVanalysis/likelihood_code/Simulation/Mrk501/Manel_flares/alpha/with_Resol/Manel_flare_alpha20_%i.txt", l); //Poner dónde están los flares
      in.open(fname);
      std::cout << "Flare: " << fname << std::endl;
      
      while(1)
      {
    	  in >> v1 >> v2;
    	  if (!in.good()) break;
    	  t[Npoints] = v1;
    	  E[Npoints] = v2;
    	  Npoints++;
      }
      
      //++++++++++++++++++++++++++++++++++++++++++  LOOP FOR ALPHA  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      
      //Initialization of variables used of all the iterations
      Double_t alpha_min = -2.;
	  Double_t alpha_max = 45.;
	  Double_t alpha = alpha_min;
	  Double_t alpha_step = 1.;
      memset(Dk, 0, numOfFlares);
      memset(alphas, 0, numOfFlares);
      memset(tcorr, 0, Np);
      memset(Comp, 0, bins);
      loop = 0;
      
      while(alpha <=alpha_max)
      {
//    	  std::cout << "alpha: " << alpha << " in loop number " << loop << std::endl;
//     	  if(alpha == 0.0){

		  //The histograms are created and deleted in every iteration
		  TH1D* nLE = new TH1D ("nLE", "nLE", bins, tmin, tmax);
		  TH1D* nHE = new TH1D ("nHE", "nHE", bins, tmin, tmax);
		  TH1D* nLE_CDF = new TH1D ("nLE_CDF", "nLE_CDF", bins, tmin, tmax);
		  TH1D* nHE_CDF = new TH1D ("nHE_CDF", "nHE_CDF", bins, tmin, tmax);


    	  //Apply correction
    	  Double_t QGDelay= QGfixedPart*alpha/PlanckScale;
    	  Double_t tcorr[Np];
    	  for (Int_t i = 0; i<Np; i++)
    	  {
    		  tcorr[i]= t[i]-QGDelay*E[i];
    	  }
      

    	  //Divide in energy and time
    	  Int_t iLE = 0 ;
    	  Int_t  iHE = 0;
    	  Double_t cutInE = 230; //Choose the energy for the cut (GeV)
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

    	//std::cout << " Corrected LE photons = " << iLE-1 << std::endl;
    	//std::cout << " Corrected HE photons = " << iHE-1 << std::endl;
	  

	  
	  //Compute mean energie and times
	  Sum_E_LE = 0;
	  Sum_E_HE = 0;
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
	  //std::cout << "Corrected Mean Energy in LE (TeV) = " << Sum_E_LE/(iLE-1) << std::endl;
	  //std::cout << "Corrected Mean Energy in HE (TeV) = " << Sum_E_HE/(iHE-1) << std::endl;
	  //std::cout << "Corrected Mean time in LE (s) = " << Sum_t_LE/(iLE-1) << std::endl;
	  //std::cout << "Corrected Mean time in HE (s) = " << Sum_t_HE/(iHE-1) << std::endl;
	  

	  
	  //Fill histograms
	  for ( Int_t a = 0; a < iLE; a++)
	  {
	      nLE->Fill(LE_t[a]);
	  }
	  for ( Int_t b = 0; b < iHE; b++)
	  {
	      nHE->Fill(HE_t[b]);
	  }
	  

	  
	  //Build CDFs
	  Double_t add_LE = 0;
	  Double_t add_HE = 0;
	  for(Int_t i = 1; i<=bins; i++)
	  {
	      Double_t value_LE = nLE->GetBinContent(i);
	      Double_t value_HE = nHE->GetBinContent(i);
	      nLE_CDF->SetBinContent(i,value_LE+add_LE);
	      nHE_CDF->SetBinContent(i,value_HE+add_HE);
	      add_LE+= value_LE;
	      add_HE+= value_HE;
	  }
	  nLE_CDF->Scale(1./nLE_CDF->Integral());
	  nHE_CDF->Scale(1./nHE_CDF->Integral());
	  
//	  TCanvas* cosa = new TCanvas("cosa", "Histrograms", 800, 800);
//	  cosa->Divide(2,3);
//	  cosa->cd(1);
//	  nLE->SetXTitle("time(s)");
//	  nLE->Draw();
//	  cosa->cd(2);
//	  nHE->SetXTitle("time(s)");
//	  nHE->Draw();
//	  cosa->cd(3);
//	  nLE_CDF->SetXTitle("time(s)");
//	  nLE_CDF->Draw();
//	  cosa->cd(4);
//	  nHE_CDF->SetXTitle("time(s)");
//	  nHE_CDF->Draw("same");
//	  cosa->cd(5);
//	  //nLE_CDF->SetXTitle("time(s)");
//	  nLE->Draw();
//	  nHE->SetLineColor(2);
//	  nHE->Draw("same");
//	  cosa->cd(6);
//	  //nHE_CDF->SetXTitle("time(s)");
//	  nLE_CDF->Draw();
//	  nHE_CDF->SetLineColor(2);
//	  nHE_CDF->Draw("same");

	  
	  //Comparison of CDFs
	  for(Int_t c = 0; c<bins; c++)
	  {
	      Comp[c]= TMath::Abs((nLE_CDF->GetBinContent(c+1))-(nHE_CDF->GetBinContent(c+1)));
	  }
	  Double_t max =0.;
	  Int_t ind;
	  for (Int_t c=0; c<bins; c++)
	  {
	      if(max<=Comp[c])
	      {
	    	  max = Comp[c];
	    	  ind = c;
	      }
	  }
//	   std::cout << "The corrected  maximum or Dk is: " << max << " located at the corrected bin " << ind << std::endl;
	  Dk[loop]= max;
	  alphas[loop]= alpha;
//	   std::cout << "Dk: " << Dk[loop] << std::endl;
	  
	  delete nLE;
	  delete nHE;
	  delete nLE_CDF;
	  delete nHE_CDF;
	  
//	   }//This is put for tests of alpha
	  
	  alpha = alpha + alpha_step;
	  loop++;
	}

      //++++++++++++++++++++++++++++++++++++++++++  END LOOP FOR ALPHA  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    //Look for the value of alpha giving a minimum Dk
    Int_t index=0;
    Double_t bestalpha, bestmax;
    bestmax = 1000;
    for (Int_t y=0; y<loop; y++)
    {
    	if (Dk[y]<=bestmax)
    	{
    		bestmax = Dk[y];
    		bestalpha = alphas[y];
    	}
    }
    std::cout << "FINAL RESULT: Best alpha is " << bestalpha << " with a Dk: " << bestmax <<  " in flare " << l << std::endl;
    results[l] = bestalpha;
    std::cout << "Result alpha loop: " << results[l] << std::endl;
//     }//This is for tests for l
  }

  //++++++++++++++++++++++++++++++++++++++++++  END LOOP FOR FLARES  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  //Computation of the final results, the errors and the histogram
  TH1D* final = new TH1D("All flares results", "All flares results", bins,-4., 50.);
  Double_t mean=0;
  Double_t sum=0;
  Double_t dev = 0;
  Double_t error = 0.;
  for(int i=0; i<numOfFlares; i++){
	  sum+= results[i];
	  final->Fill(results[i]);
  }
  std::cout << "sum: " << sum << std::endl;
  mean = sum/numOfFlares;
  for(int i=0; i<numOfFlares; i++){
	  dev+= pow((results[i]-mean),2);
   }
   error = sqrt(dev/numOfFlares);
   std::cout << "Mean: " << mean << " Error: " << error << std::endl;
   final->Draw();






}


