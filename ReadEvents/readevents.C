/*
 * readevents.C
 *
 *  Created on: Feb 17, 2016
 *      Author: lnogues
 */

#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TStyle.h>
#include<iostream>
#include<fstream>
#include<istream>
using namespace std;

const Int_t Np = 100000;   // Number of events (hardcoded)
float row[Np];
float alpha[Np];
float Ep[Np];
float tp[Np];
float MJD[Np];

Double_t powerLawFit(Double_t* x, Double_t* par)
{
  Double_t func2 = par[0]*TMath::Power(x[0],-par[1]);
  return func2;
}

void readevents(char *fname)
{

	/*//settings AJ
	gStyle->SetOptStat(0);
	gROOT->SetStyle("Plain");
	gStyle->SetPalette(1);
	gStyle->SetTitleAlign(13);

	// Canvas
	gStyle->SetCanvasColor(10);

	// Frame
	 gStyle->SetFrameBorderMode(0);
	 gStyle->SetFrameFillColor(0);

	// Pad
	 gStyle->SetPadBorderMode(0);
	 gStyle->SetPadColor(0);
	 gStyle->SetPadTopMargin(0.07);
	 gStyle->SetPadLeftMargin(0.13);
	 gStyle->SetPadRightMargin(0.11);
	 gStyle->SetPadBottomMargin(0.1);
	 gStyle->SetPadTickX(1); //make ticks be on all 4 sides.
	 gStyle->SetPadTickY(1);

	// histogram
	 gStyle->SetHistFillStyle(0);
	 gStyle->SetOptTitle(0);

	// histogram title
	 gStyle->SetTitleSize(0.22);
	 gStyle->SetTitleFontSize(2);
	 gStyle->SetTitleFont(42);
	 gStyle->SetTitleFont(62,"xyz");
	 gStyle->SetTitleYOffset(1.0);
	 gStyle->SetTitleXOffset(1.0);
	 gStyle->SetTitleXSize(0.04);
	 gStyle->SetTitleYSize(0.04);
	 gStyle->SetTitleX(.15);
	 gStyle->SetTitleY(.98);
	 gStyle->SetTitleW(.70);
	 gStyle->SetTitleH(.05);

	// statistics box
	 gStyle->SetStatFont(42);
	 gStyle->SetStatX(.91);
	 gStyle->SetStatY(.90);
	 gStyle->SetStatW(.15);
	 gStyle->SetStatH(.15);

	// axis labels
	 gStyle->SetLabelFont(42,"xyz");
	 gStyle->SetLabelSize(0.035,"xyz");
	 gStyle->SetGridColor(16);
	 gStyle->SetLegendBorderSize(0);
     gStyle->SetOptStat(111111);
     gStyle->SetOptFit(111111);*/

	gStyle->SetOptStat(0);
    gStyle->SetTitleXSize(0.05);
	 gStyle->SetTitleYSize(0.05);
	 gStyle->SetTitleSize(0.30);


//---------------------------------------Read the events--------------------------------------

  ifstream in;
  in.open(fname);
  float v1,v2,v3,v4,v5;
  int Npoints = 0;
  while(1)
    {
      in >> v1 >> v2 >> v3 >> v4 >> v5;
      if (!in.good()) break;
      row[Npoints] = v1;
      alpha[Npoints] = v2;
      Ep[Npoints] = v3;
      tp[Npoints] = v4;
      MJD[Npoints] = v5;
      Npoints++;
    }
  cout << "Npoints " << Npoints << endl;
  /*cout << row[0] << " " << alpha[0] << " " << Ep[0] << " " << tp[0] << " " << MJD[0] << endl;
  cout << row[1] << " " << alpha[1] << " " << Ep[1] << " " << tp[1] << " " << MJD[1] << endl;
  cout << row[2] << " " << alpha[2] << " " << Ep[2] << " " << tp[2] << " " << MJD[2] << endl;*/

  //-----------------------------Maximum and minimum time and energy---------------------------

  //TIME
  float maxt, mint;
  maxt = mint = tp[0];
  for ( int a = 1; a < Npoints; a++)
     {
       if(tp[a]<mint)
 	{
           mint=tp[a];
 	}
       if(tp[a]>maxt)
         {
           maxt=tp[a];
         }
     }
   cout << "The minimum time is " << mint << " and the maximum time is " << maxt << "s"<< endl;

   //ENERGY
     float maxE, minE;
     maxE = minE = Ep[0];
     for ( int a = 1; a < Npoints; a++)
        {
          if(Ep[a]<minE)
    	{
              minE=Ep[a];
    	}
          if(Ep[a]>maxE)
            {
              maxE=Ep[a];
            }
        }

     cout << "The minimum energy is " << minE << " and the maximum energy is " << maxE << "GeV"<< endl;

     //ALPHA
          float maxalpha, minalpha;
          maxalpha = minalpha = alpha[0];
          for ( int a = 1; a < Npoints; a++)
             {
               if(alpha[a]<minalpha)
         	{
            	   minalpha=alpha[a];
         	}
               if(alpha[a]>maxalpha)
                 {
            	   maxalpha=alpha[a];
                 }
             }

      cout << "The minimum alpha is " << minalpha << " and the maximum alpha is " << maxalpha << " degrees"<< endl;


      //-----------------------------Histograms for time and energy---------------------------


      //TF1* plfit = new TF1("plfit", powerLawFit, 0.15, 7., 2);
     TH1F* time = new TH1F("time", "Signal (150GeV<E<12000GeV)", 12, mint, ceil(maxt));
      time->SetLineWidth(2);
//      TH1F* timebkg = new TH1F("timebkg", "Background (150GeV<E<12000GeV)", 12, mint, ceil(maxt));
//      timebkg->SetLineWidth(2);
//      TH1F* energy = new TH1F("Baseline", "Baseline", 50, floor(minE/1000), ceil(maxE/1000));

//      FILE *fout;
//      fout = fopen("Mrk501_gammadata_factor10.txt","wb");
      int numberofevents= 0;
      for ( int a = 0; a < Npoints; a++)
      {
    		  if(Ep[a]>=150 && alpha[a]<10)
    		  {
//    			  fprintf(fout,"%0.8lf    %0.8lf   \n", tp[a], Ep[a]*1.1);
    			  numberofevents++;
    			  time->Fill(tp[a]);
    		  }
//    		  if(Ep[a]>=150 && Ep[a]<=12000 && alpha[a]>10)
//			  {
//				  timebkg->Fill(tp[a]);
//			  }

    	  }

      cout<<numberofevents<<endl;


      TCanvas* t = new TCanvas();
//      t->SetTitle("Time distributions");
//      t->Divide(1,2);
//      t->cd(1);
      time->GetXaxis()->SetTitle("Time(s)");
      time->Draw();
//      t->cd(2);
//      timebkg->GetXaxis()->SetTitle("Time(s)");
//      timebkg->Draw();

//      TCanvas* E = new TCanvas();
      //E->SetTitle("Baseline energy fit");
      //E->SetLogx();
      //E->SetLogy();
//      energy->GetXaxis()->SetTitle("Energy(TeV)");
      //energy->Fit(plfit, "QR");
//      energy->Draw();

      //-----------------------------Energy bands---------------------------

     /*TH1F* energy1 = new TH1F("energy1", "energy1", 12, mint, ceil(maxt));
      energy1->SetLineWidth(2);
      TH1F* energy2 = new TH1F("energy2", "energy2", 12, mint, ceil(maxt));
      energy2->SetLineWidth(2);
      TH1F* energy3 = new TH1F("energy3", "energy3", 12, mint, ceil(maxt));
      energy3->SetLineWidth(2);
      TH1F* energy4 = new TH1F("energy3", "energy4", 12, mint, ceil(maxt));
      energy4->SetLineWidth(2);


      for ( int a = 0; a < Npoints; a++)
         {
       		  if(alpha[a]<10) //Only signal
       		  {
       			  if(Ep[a]>150&&Ep[a]<250){
       				energy1->Fill(tp[a]);
       			  }
       			if(Ep[a]>250&&Ep[a]<600){
					energy2->Fill(tp[a]);
				  }
       			if(Ep[a]>600&&Ep[a]<1200){
					energy3->Fill(tp[a]);
				  }
       			if(Ep[a]>1200){
					energy4->Fill(tp[a]);
				  }

   			  }

       	  }*/



      	/*TCanvas* bands = new TCanvas();
      	bands->SetTitle("Time distributions");
      	bands->Divide(1,4);
      	bands->cd(1);
		energy1->GetXaxis()->SetTitle("Time(s)");
		energy1->SetTitle("150-250 GeV");
		energy1->SetLineColor(6);
		energy1->Draw();
		bands->cd(2);
		energy2->GetXaxis()->SetTitle("Time(s)");
		energy2->SetTitle("250-600 GeV");
		energy2->SetLineColor(4);
		energy2->Draw();
		bands->cd(3);
		energy3->GetXaxis()->SetTitle("Time(s)");
		energy3->SetTitle("600-1200 GeV");
		energy3->SetLineColor(3);
		energy3->Draw();
		bands->cd(4);
		energy4->GetXaxis()->SetTitle("Time(s)");
		energy4->SetTitle("1200-12000 GeV");
		energy4->SetLineColor(2);
		energy4->Draw();*/



    //LC functions

      /*TF1 *baseline = new TF1("baseline", "pol0", mint, maxt);
      baseline->SetParameter(0,38.93);
      TF1* flare = new TF1("flare", "gaus", mint, maxt);
      flare->SetParameters(76.96, 2009., 689.6);

      TCanvas* comLC = new TCanvas();
      baseline->SetTitle("Time Distribution");
      baseline->GetXaxis()->SetTitle("Time(s)");
      baseline->SetLineColor(3);
      baseline->Draw();
      flare->SetLineColor(2);
      flare->Draw("same");

      cout << "Area ratio: " << baseline->Integral(mint, maxt)/flare->Integral(mint, maxt) << endl;*/



      //-----------------------------Signal to Background ratio---------------------------

      /*TH1F* final = new TH1F("final", "final", 8, TMath::Log10(300), TMath::Log10(8000));
      TH1F* alph = new TH1F("alpha", "alpha", 20, 0, 90);

      Double_t Nsig[10], Nbkg[10], SBratio[10];
      Int_t times = 1;
      Double_t suma = 0;

      for(Int_t i=1; i<=8; i++)
      {
    	  Int_t bin = i;
    	  Nsig[bin-1]= 0;
    	  Nbkg[bin-1]= 0;
    	  SBratio[bin-1]= 0;

    	  for ( int a = 0; a < Npoints; a++)
    	  {
    		  if(1200.<=tp[a]&&tp[a]<=maxt)
    		  {
    			  if(final->GetBinLowEdge(bin)<=TMath::Log10(Ep[a]) && TMath::Log10(Ep[a])<=final->GetBinLowEdge(bin+1))
    			  {
    				  //alph->Fill(alpha[a]);

    				  if(alpha[a]<=10.)
    				  {
    					  Nsig[bin-1]++;
    				  }
    				  if(alpha[a]>10.)
    				  {
    					  Nbkg[bin-1]++;
    				  }

    			  }
    		  }
    	  }


    	  SBratio[bin-1]= Nsig[bin-1]/Nbkg[bin-1];
    	  cout << "Bin: " << bin << " Nsig: " << Nsig[bin-1]<< " Nbkg: " << Nbkg[bin-1] << endl;
    	  cout << "SBratio: " << SBratio[bin-1] << endl;
    	  if(SBratio[bin-1] <=1000)
    	  {
    	  final->SetBinContent(bin, SBratio[bin-1]);
    	  suma+= SBratio[bin-1];
    	  times++;
    	  }
      }


      cout<< "Mean SB ratio: " << suma/times << endl;

      TGraph* funene = new TGraph();

      for(Int_t i=1; i<=9; i++)
      {
    	  funene->SetPoint(i-1, TMath::Power(10,final->GetBinCenter(i)), final->GetBinContent(i));
      }


      TCanvas* alphcanvas = new TCanvas();
      alphcanvas->SetTitle("SBratio");
      final->SetTitle("S/B ratio vs Energy");
      final->SetName("SBfit");
      final->GetXaxis()->SetTitle("Log10(E)(GeV)");
      final->GetYaxis()->SetTitle("S/B ratio");
      final->Fit("pol1");
      final->Draw("EP");

      TF1* line = final->GetFunction("pol1");
      line->SetLineColor(2);
      line->Draw("same");

      /*Double_t ini = TMath::Log10(300.);
      Double_t fin = TMath::Log10(5000.);
      cout << ini << fin << endl;
      Double_t fitMean = line->Integral(ini, fin)/(fin-ini);
      cout << "FitMean: " << fitMean << endl;

      TString rootoutfile("/home/lnogues/workspace_cpp/LIVanalysis/Likelihood/Mrk501/SBratio.root");
      TFile* outputfile = new TFile(rootoutfile, "recreate");
      final->Write();*/

      /*TCanvas* funenecanvas = new TCanvas();
      funenecanvas->SetTitle("SBratio vs Energy");
      funenecanvas->SetLogx();
      //funenecanvas->SetLogy();
      funene->SetLineColor(2);
      funene->SetLineWidth(2);
      funene->SetTitle("S/B vs Energy");
      funene->GetXaxis()->SetTitle("E(GeV)");
      funene->GetYaxis()->SetTitle("S/B ratio");
      funene->Draw("AL");*/


}
