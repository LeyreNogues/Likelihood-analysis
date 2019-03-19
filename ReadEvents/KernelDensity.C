/*
 * KernelDensity.C
 *
 *  Created on: Apr 24, 2018
 *      Author: lnogues
 */

#include <Rtypes.h>
#include <TAttLine.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TKDE.h>
#include <vector>

#ifndef __CINT__
#include "Math/DistFunc.h"
#endif
#include "TLegend.h"
#include<iostream>
#include<fstream>
#include<istream>

using namespace std;


void exampleTKDE(int n = 10) {

 // generate some gaussian points


   int nbin = 200;
   Double_t extra_time = 8000.;
   double xmin = 0-extra_time;
   double xmax = 13100+extra_time;
//
   TH1D* h1 = new TH1D("h1","h1",nbin,xmin,xmax);

//   // generate some points with bi-gaussian distribution
//
   std::vector<double> data;
//   gRandom->SetSeed(0);
//   for(int i = 0; i < n; ++i){
//	 data.push_back(gRandom->Uniform(xmin,xmax));
//	 h1->Fill(data[i]);
//   }



   //Read real data
//	TFile* HistoFile = TFile::Open("../ReadEvents/FinalEventLists/All_Events/AllEvents_ONhistograms_Mrk421_2014_Mireia.root", "READ");
//	if (!HistoFile) {
//		cout << "ERROR: file " << HistoFile << " not found..."	<< endl;
//	}
//	TH1D* lightcurve = dynamic_cast<TH1D*>(HistoFile->Get("LC_insecs"));
//	if (!lightcurve) {
//		cout << "LC_insecs not found" << endl;
//	}
//
//	Int_t numberOfBins = lightcurve->GetNbinsX();
//	cout << "Number of bins in initial histo: " << numberOfBins << endl;
//	const Double_t* arrayBins = lightcurve->GetXaxis()->GetXbins()->GetArray();


	ifstream in;
	in.open("./CompleteList_Kazuma_120_OFF.txt");
	Double_t v1,v2;
//	TH1D* h1 = new TH1D("Time distribution", "Time distribution", numberOfBins, arrayBins);
	int Npoints = 0;
	while(1)
	{
		in >> v1 >> v2;
		if (!in.good()) break;
		data.push_back(v1);
		h1->Fill(v1);
		Npoints++;
	}

	cout << "Npoints(Number of events): " << Npoints << endl;
	n=Npoints;

//	TCanvas* canvas = new TCanvas();
//	h1->Draw();

	TCanvas* TKDE_canvas = new TCanvas();
   // scale histogram
   h1->Scale(1./h1->Integral(),"width" );
   h1->SetStats(false);
   h1->GetXaxis()->SetTitle("Time(s)");
//   h1->SetTitle("Histo");
////   h1->GetYaxis()->SetRange(0, 80e-6);
////   h1->SetMaximum(80e-6);
////   h1->SetMinimum(70e-6);
   h1->Draw("");

//	  TF1 * f1 = new TF1("f1","TMath::Landau(x,[0],[1],0)+TMath::Landau(x,[2],[3],0)+TMath::Landau(x,[4],[5],0)",xmin,xmax);
//	  f1->SetParameters(0.2,1., 3., 0.7, -3, 1.2);
//	  f1->SetLineColor(kGreen+2);
//	  f1->Draw("");



//	   for (int i = 0; i < n; ++i) {
//	         data[i] = f1->GetRandom();
//	         h1->Fill(data[i]);
//	   }


   // drawn true normalized density
//   TF1 * f1 = new TF1("f1","0.4*ROOT::Math::normal_pdf(x,1,2)+0.6*ROOT::Math::normal_pdf(x,1.5,7)",xmin,xmax);
//   f1->SetLineColor(kGreen+2);
//   f1->Draw("SAME");


   ///////////////////////////////////////////////////////////////////////////

//    create TKDE class
   double rho = 0.5; //default value
//   double rho2 = 0.1; //default value
   TKDE * kde = new TKDE(n, &data[0], xmin,xmax, "", rho);
//   TKDE * kde2 = new TKDE(n, &data[0], xmin,xmax, "", rho2);
   //kde->Draw("ConfidenceInterval@0.95 Same");
//
//
//   // KDE options
////    int nbinKDE = 1000;
//   // TKDE::EBinning bin;
    TKDE::EIteration iter = TKDE::kFixed;

//////    TKDE::EMirror mirror = TKDE::kMirrorBoth; //Problem
////    TKDE::EMirror mirror = TKDE::kNoMirror;
//////    TKDE::EMirror mirror = TKDE::kMirrorAsymBoth;
//////      TKDE::EMirror mirror = TKDE::kMirrorAsymLeft;
//////      TKDE::EMirror mirror = TKDE::kMirrorAsymLeftRight; //Problem
//////      TKDE::EMirror mirror = TKDE::kMirrorAsymRight;
//////      TKDE::EMirror mirror = TKDE::kMirrorLeft;
//////      TKDE::EMirror mirror = TKDE::kMirrorRight;
//////      TKDE::EMirror mirror = TKDE::kMirrorLeftAsymRight;
////
////
////
////
////
//////    if (nbinKDE> 0) kde->SetNBins(nbinKDE);
    if (iter != TKDE::kAdaptive) kde->SetIteration(iter);
//    if (iter != TKDE::kAdaptive) kde2->SetIteration(iter);

//////    if (mirror != TKDE::kNoMirror) kde->SetMirror(mirror);
////    kde->SetMirror(mirror);
//    kde->SetIteration(iter);
//
//    kde->Draw("SAME");
//
//
//   //alternative way instead of using kde->Draw()
    TF1 * const hk = kde->GetFunction(5000);
//    TF1 * const hk2 = kde2->GetFunction(5000);

    hk->SetLineColor(kRed);
//    hk2->SetLineColor(kBlue);

//    hk->SetNpx(2000);
    hk->Draw("same");
//    hk2->Draw("same");
//
////    TF1 * fl = kde->GetLowerFunction(0.684);
////    TF1 * fu = kde->GetUpperFunction(0.684);
////    fl->SetNpx(2000);
////    fu->SetNpx(2000);
////    fl->Draw("SAME");
////    fl->SetLineColor(kBlue-5);
////    fu->SetLineColor(kBlue-5);
////    fu->Draw("SAME");
//
//
//   TLegend * legend = new TLegend(0.6,0.7,0.9,0.95);
//   legend->AddEntry(h1,"Data histogram");
//   legend->AddEntry(hk,"KDE, #rho = 1.0");
//   legend->AddEntry(hk2,"KDE, #rho = 0.1");
////   legend->AddEntry(fl,"TKDE - #sigma");
////   legend->AddEntry(fu,"TKDE + #sigma");
//   legend->Draw();
//   return;
//
//   gPad->Update();


    //Save TF1 for model
    TFile* TKDE_outputfile = new TFile("TKDE_Kazuma_120_0.5_OFF.root","recreate");
    hk->SetName("TKDE_model");
    hk->Write();
    TKDE_canvas->SetName("TKDE_canvas");
    TKDE_canvas->Write();
    TKDE_outputfile->Close();


}




