/*
 * Fit_histos.C
 *
 *  Created on: Apr 12, 2018
 *      Author: lnogues
 */

#include<TFile.h>
#include<TCanvas.h>
#include<TH1.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TF1.h>


using namespace std;

Double_t limit_QG(Int_t type, Double_t lag, Double_t stat, Double_t sys, Double_t nsig = 2.){

	Double_t z=0.031;
	Double_t H0=70.;
	Double_t PlanckScale=1.2*TMath::Power(10,19);
	Double_t QGpar=0.;
	Double_t PcTom=3.085*TMath::Power(10,19);
	Double_t QGfixedPart=z*PcTom/H0;

	Double_t lag_nsig = lag + nsig * sqrt(stat * stat + sys * sys); //Final alpha

	Double_t res = 0.; //resolution
	  if(type==1){ //Linear
		res= PlanckScale/lag_nsig;
	  }
	  if(type==2){ //Quadratic
		res= TMath::Sqrt((3/2.)*(PlanckScale/lag_nsig));
	  }
	  return res;
}

 Double_t Double_Gaus(Double_t* x, Double_t*  par){
	Double_t simgaus = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	Double_t simgaus2 = par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5],2));
	return simgaus + simgaus2;
}

 Double_t One_Gaus(Double_t* x, Double_t*  par){
	Double_t simgaus = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	return simgaus;
}

void Fit_histos(){

TFile* RootFile = TFile::Open("histo_alpha-10_Qua.root", "READ");
if (!RootFile) {
	cout << "ERROR: file " << RootFile << " not found..."	<< endl;
}
TCanvas *c1 = (TCanvas*)RootFile->Get("c1");
TH1D* histo = (TH1D*)c1->FindObject("htemp");
histo->SetName("alpha_histo");
TCanvas* canvas = new TCanvas();
histo->Draw();

//Fit with a double Gaussian
TF1* function = new TF1("function", Double_Gaus, -30, 30, 6);
function->SetNpx(1000);
function->SetParameters(45, -10, 2, 10, 0, 2);
function->SetParName(0, "C1");
function->SetParName(1, "Mean1");
function->SetParName(2, "Sigma1");
function->SetParName(3, "C2");
function->SetParName(4, "Mean2");
function->SetParName(5, "Sigma2");
histo->Fit(function);
//TCanvas* canvas2 = new TCanvas();
//function->Draw("");




}



