/*
 * DrawResults.C
 *
 *  Created on: Jul 12, 2018
 *      Author: lnogues
 */

#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TStyle.h>
#include<iostream>
#include<fstream>
#include<istream>
#include<TVectorD.h>
#include<TSystem.h>

using namespace std;

void DrawResults(){

	 gSystem->Load("libMathMore");

	Double_t sigma[10], width[10];
	Double_t Diff[10];
	Double_t Se[10];
	Double_t Se_RMS[10];
	Double_t sigma_RMS[10] = {0};
	Double_t ratio = 1.0279;

	sigma[0] = 3;
	sigma[1] = 4;
	sigma[2] = 5;
	sigma[3] = 6;
	sigma[4] = 7;
	sigma[5] = 8;
	sigma[6] = 9;
	sigma[7] = 10;
	sigma[8] = 11;
	sigma[9] = 12;

	for (Int_t i=0; i<10;i++){
		width[i]=sigma[i]*sigma[i]/ratio;
	}


	Diff[0] = 74.3048;
	Diff[1] = 27.7645;
	Diff[2] = 11.7618;
	Diff[3] = 6.0738;
	Diff[4] = 4.6647;
	Diff[5] = 2.1892;
	Diff[6] = 1.5813;
	Diff[7] = 1.1456;

	Se[0] = 0.095022;
	Se[1] = 0.07742;
	Se[2] = 0.066256;
	Se[3] = 0.058748;
	Se[4] = 0.052673;
	Se[5] = 0.045649;
	Se[6] = 0.041486;
	Se[7] = 0.037302;
	Se[8] = 0.03503;
	Se[9] = 0.029986;

	Se_RMS[0] = 0.00572*2;
	Se_RMS[1] = 0.00403*2;
	Se_RMS[2] = 0.00337*2;
	Se_RMS[3] = 0.00313*2;
	Se_RMS[4] = 0.00334*2;
	Se_RMS[5] = 0.00306*2;
	Se_RMS[6] = 0.00294*2;
	Se_RMS[7] = 0.00311*2;
	Se_RMS[8] = 0.00326*2;
	Se_RMS[9] = 0.00308*2;

//	TGraph* diffGraph =  new TGraph(8, sigma, Diff);
//	diffGraph->GetXaxis()->SetTitle("sigma");
//	diffGraph->GetYaxis()->SetTitle("#sum#Delta c");
//	diffGraph->SetMarkerStyle(20);
//	diffGraph->SetTitle("");
//	TGraphErrors* SeGraph =  new TGraphErrors(10, sigma, Se, sigma_RMS, Se_RMS);
//	SeGraph->SetTitle("Se - fixed-time-rate histogram");
//	SeGraph->GetXaxis()->SetTitle("sigma");
//	SeGraph->GetYaxis()->SetTitle("Se");
//	SeGraph->SetMarkerStyle(20);
//	SeGraph->SetTitle("");
//
//	TCanvas* patata = new TCanvas();
////	patata->Divide(2,1);
////	patata->cd(1);
////	diffGraph->Draw("AP");
//	patata->cd(1);
//	SeGraph->Draw("AP");

    TF1* bessel_0= new TF1("J_0", "ROOT::Math::cyl_bessel_i(0,2*x)*TMath::Exp(-2*x)", 0, 110);
    TF1* bessel_1= new TF1("J_1", "ROOT::Math::cyl_bessel_i(1,2*x)*TMath::Exp(-2*x)", 0, 110);
    TF1* sumbessel =  new TF1("SUM" , "ROOT::Math::cyl_bessel_i(0,2*x)*TMath::Exp(-2*x)+ROOT::Math::cyl_bessel_i(1,2*x)*TMath::Exp(-2*x)", 0 ,110);
    TCanvas* besselcanvas =  new TCanvas();
    bessel_0->Draw();
    bessel_1->SetLineColor(1);
    bessel_1->Draw("same");
    sumbessel->SetLineColor(4);
    sumbessel->Draw("same");

    cout << "Bessel_0 in 100: " << bessel_0->Eval(100.)<< endl;
    cout << "Bessel_1 in 100: " << bessel_1->Eval(100.)<< endl;
    cout << "Bessel_suma in 100: " << sumbessel->Eval(100.)<< endl;

    cout << "Bessel_0 in 200: " << bessel_0->Eval(110.)<< endl;
    cout << "Bessel_1 in 200: " << bessel_1->Eval(110.)<< endl;
    cout << "Bessel_suma in 200: " << sumbessel->Eval(110.)<< endl;

	TCanvas* rosa = new TCanvas();
	TGraph* graph = new TGraph(9);
	for(Int_t i = 0; i< 9; i++){
		Double_t mu=sigma[i]*sigma[i];
		cout <<"mu: " << mu << endl;
		Double_t y = (2*mu/width[i])*sumbessel->Eval(mu);
		graph->SetPoint(i, sigma[i], y);
	}
	graph->SetLineColor(2);
	graph->Draw("AL");

//	TGraph* graph_ON = new TGraph(10);
//	for(Int_t i = 0; i< 10; i++){
//		Double_t mu=sigma[i]*sigma[i];
//		Double_t y = (2*mu/width[i])*TMath::Exp(-2*mu)*(bessel_0->Eval(2*mu)+bessel_1->Eval(2*mu));
//		graph_ON->SetPoint(i, sigma[i], y);
//	}
//	TCanvas* canvas2 = new TCanvas();
//	graph_ON->SetLineColor(3);
//
	Int_t N=100;
	TGraph* graph_ON2 = new TGraph(N);
	for(Int_t i = 0; i< N; i++){
		Double_t mu=i*i;
		Double_t y = (2*ratio)*TMath::Sqrt(1./TMath::Pi()/mu);
//		Double_t y = (4/wi)*TMath::Exp(-mu)*TMath::Sqrt(mu/2/TMath::Pi());
		graph_ON2->SetPoint(i, i, y);
	}
////	TCanvas* canvas22 = new TCanvas();
	graph_ON2->SetLineColor(9);
//    graph_ON2->Draw("AL");
	graph_ON2->Draw("L");

//	TCanvas* canvas2 = new TCanvas();
//	graph->Draw("ALP");





}

