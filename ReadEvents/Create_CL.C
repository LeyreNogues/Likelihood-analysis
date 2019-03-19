/*
 * Create_CL.C
 *
 *  Created on: Jul 28, 2018
 *      Author: lnogues
 *
 *      This codes get a chain made with the distribution of results, normalize it
 *      and gives us the CL for a given measurement
 */

#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include<iostream>
#include<fstream>
#include<istream>

using namespace std;


void Create_CL(){

	//Read the histogram

	TFile* HistoFile = TFile::Open("./Coverage/Canvas_Quadratic_Paper_wide.root", "READ");
	TCanvas *c1 = (TCanvas*)HistoFile->Get("c1");
	TH1F* original_histo = dynamic_cast<TH1F*>(c1->FindObject("alpha"));

	TCanvas* canvas = new TCanvas();
//	original_histo->SetStats(0);
//	original_histo->Rebin(4);
//	original_histo->SetBinContent(14,197);
//	original_histo->Draw();


	//Normalize the histogram
	original_histo->Scale(1./original_histo->Integral());
	original_histo->SetTitle("Error PDF");
	original_histo->GetXaxis()->SetTitle("#alpha");
	original_histo->GetYaxis()->SetTitle("p(1)");
//	original_histo->GetIntegral();
	original_histo->Draw();

	Int_t bin1 = original_histo->FindBin(-1.26);
	Int_t bin2 = original_histo->FindBin(1.70);
	Double_t integral_before = original_histo->Integral(bin1,bin2);
	cout << "Integral before is: " << integral_before << endl;



	//Compute CL
   const Int_t nq = 2;
   Double_t xq[nq];  // position where to compute the quantiles in [0,1]
   Double_t yq[nq];  // array to contain the quantiles
   xq[0]=0.05;
//   xq[1]=0.1;
//   xq[2]=0.9;
   xq[1]=0.95;
   original_histo->GetQuantiles(nq,yq,xq);
   //show the original histogram in the top pad
   TCanvas *c2 = new TCanvas("c2","demo quantiles",10,10,700,900);
   c2->Divide(1,2);
   c2->cd(1);
   original_histo->Draw();
   // show the quantiles in the bottom pad
   c2->cd(2);
   gPad->SetGrid();
   TGraph *gr = new TGraph(nq,xq,yq);
   gr->SetTitle("Quantiles");
   gr->GetYaxis()->SetTitle("#alpha");
   gr->SetMarkerStyle(21);
   gr->Draw("alp");

   for(Int_t i=0; i<nq; i++){
	   std::cout << "Quantil: " << xq[i] << " tau value: " << yq[i] << std::endl;
   }






}


