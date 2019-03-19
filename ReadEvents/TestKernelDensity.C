/*
 * KernelDensity.C
 *
 *  Created on: Apr 24, 2018
 *      Author: lnogues
 */

#include "TH1.h"
#include "TF1.h"
#include "TKDE.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
//#include "TStopwatch.h"
#include "TRandom.h"
#ifndef __CINT__
#include "Math/DistFunc.h"
#endif
#include "TLegend.h"
#include<iostream>
#include<fstream>
#include<istream>

using namespace std;


void exampleTKDE(int n = 12520) { //Kazuma events

 // generate some gaussian points


   int nbin = 200;
   Double_t extra_time = 0.;
   double xmin = 0-extra_time;
   double xmax = 13100+extra_time;

   TH1D* h1 = new TH1D("h1","h1",nbin,xmin,xmax);

   // generate some points with bi-gaussian distribution

//   std::vector<double> data;
//   gRandom->SetSeed(0);
//   for(int i = 0; i < n; ++i){
//	 data.push_back(gRandom->Uniform(xmin,xmax));
////	 h1->Fill(data[i]);
//   }




	TCanvas* TKDE_canvas = new TCanvas();
   // scale histogram
//   h1->Scale(1./h1->Integral(),"width" );
   h1->SetTitle("Histo");
   h1->SetMaximum(0.1e-3);
   h1->SetMinimum(0.05e-3);
   h1->Draw("");




   ///////////////////////////////////////////////////////////////////////////

   TGraphErrors* Qgraph =  new TGraphErrors();
   Int_t graph_index = 0;
   std::vector<double> data;
   const Int_t numberOfSamples = 100;
   double rho;
   Double_t Q_factor[numberOfSamples];
   Double_t meanQfactor = 0.;
   Double_t RMSQfactor = 0.;

   gRandom->SetSeed(0);

   for(Double_t i = 0.05; i<=1.1 ; i=i+0.05){

	   meanQfactor=0.;
	   RMSQfactor=0.;

	   for(Int_t k = 0; k<numberOfSamples; k++){

		   Q_factor[k] = 0.;

		   cout << "HERE STARTS ONE SAMPLE " << endl;
		   cout << "K =  " << k << endl;


		   for(int p = 0; p < n; ++p){
			 data.push_back(gRandom->Uniform(xmin,xmax));
		   }

//		   cout << "First component: " << data.at(0) << endl;

		   // create TKDE class
		   rho = i; //default value
		   cout << "rho: " << rho << endl;
		   TKDE * kde = new TKDE(n, &data[0], xmin,xmax, "", rho);


		//   //alternative way instead of using kde->Draw()
			TF1 * const hk = kde->GetFunction(5000);
			hk->SetLineColor(kRed);
		//    hk->SetNpx(2000);
			hk->Draw("same");

			//----Getting Maximum and Minimum peaks and Mean Value
			Double_t up_bound = xmax-3000.;
			Double_t down_bound = xmin+3000;
			Double_t max_Peak = 0.;
			Double_t min_Peak = 10000.;
			Double_t mean_value = 0.;
			Int_t numberOfPoints = 100;
			Double_t step = (up_bound-down_bound)/numberOfPoints;
//			cout << "Bounds: " << down_bound << " " << up_bound << endl;
//			cout << "Number Of Points: " << numberOfPoints << endl;
//			cout << "Step: " << step << endl;

			for(Int_t j=0; j<numberOfPoints; j++){
				Double_t point = down_bound+j*step;
				Double_t value = hk->Eval(point);
				if(max_Peak<value) max_Peak = value;
				if(min_Peak>value) min_Peak =  value;
				mean_value+=value;
			}

			mean_value = mean_value/numberOfPoints;
			Q_factor[k]=(max_Peak-min_Peak)/mean_value;
			meanQfactor+=Q_factor[k];
			RMSQfactor+=(Q_factor[k])*(Q_factor[k]);

			cout << "min_Peak: " << min_Peak << endl;
			cout << "max_Peak: " << max_Peak << endl;
			cout << "mean_value: " << mean_value << endl;
			cout << "Q_factor: " << Q_factor[k] << endl;
			cout << "meanQfactor: " << meanQfactor << endl;
			cout << "RMSQfactor: " << RMSQfactor << endl;


			kde->Delete();
			data.clear();
	   }

	   	meanQfactor = (double)meanQfactor/numberOfSamples;
	   	Double_t Stand_dev = TMath::Sqrt((RMSQfactor-(meanQfactor*meanQfactor*numberOfSamples))/(numberOfSamples-1.));
	   	cout << "Final point: " << meanQfactor << endl;
	   	cout << "Final error: " << Stand_dev << endl;
		Qgraph->SetPoint(graph_index, rho, meanQfactor);
		Qgraph->SetPointError(graph_index, 0., Stand_dev);
		graph_index++;

   }

   gStyle->SetOptFit(1101);
   TCanvas* Qcanvas =  new TCanvas();
   Qgraph->SetMarkerStyle(20);
   Qgraph->GetXaxis()->SetTitle("Smoothing parameter (rho)");
   Qgraph->GetYaxis()->SetTitle("Noise Estimator (%)");
   Qgraph->Draw("ap");
   Qgraph->Fit("pol0", "","", 0.45, 1.1);



}




