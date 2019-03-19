/*
 * 2Dsim.C
 *
 *  Created on: Mar 7, 2016
 *      Author: lnogues
 *
 *  Simulation in 2D for LIV studies. All the events are used in the FIT, all used to fit the LC?
 */

#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TRandom.h>
using namespace std;

TGraph* CollAr;

Double_t tmin = 0;
Double_t tmax = 2731;
Double_t Emin = 0.15;
Double_t Emax = 11.210;

const Double_t numev = 2000;

Double_t signal2D(Double_t *x, Double_t *par)
{
	Double_t signal = TMath::Power(x[1],-2.7)+par[0]*TMath::Power(x[1],-2.4)*exp(-0.5*TMath::Power(((x[0]-par[1]-par[3]*x[1])/par[2]),2));
	Double_t area = CollAr->Eval(x[1]);
	return signal*area;
}

Double_t powerLawFit(Double_t* x, Double_t* par)
{
  Double_t func2 = par[0]*TMath::Power(x[0],-par[1]);
  return func2;
}

void twoDsim()
{
	//ROOT OPTION
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(1111);
		gStyle->SetOptTitle(1);
		gStyle->SetCanvasColor(0);
		gStyle->SetLabelOffset(0.003,"y");
		gStyle->SetLabelOffset(0.01,"x");
		gStyle->SetTitleOffset(1.0,"x");
		gStyle->SetTitleOffset(1.0,"y");
		gStyle->SetTitleXSize(0.05);
		gStyle->SetTitleYSize(0.05);
		gStyle->SetPadLeftMargin(0.15);
		gStyle->SetPadBottomMargin(0.15);
		gStyle->SetCanvasBorderMode(0);
		gStyle->SetTickLength(0.01,"y");
		gStyle->SetStatX(0.99);
		gStyle->SetStatY(0.99);
		gStyle->SetStatH(0.1);
		gStyle->SetStatW(0.1);
		gStyle->SetStatH(0.2);
		gStyle->SetStatW(0.2);
		gStyle->SetStatFontSize(0.03);

	//-------------------------------Reading of CollAr----------------------------------------

		TString myPath = "/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Mrk501/";
		const TString sNewFile = myPath+"MonoCollAr.root";
		TFile* fNewFile = new TFile(sNewFile);
		CollAr = (TGraph*)fNewFile->Get("Mono");
		/*TCanvas* area = new TCanvas();
		area->SetLogx();
		area->SetLogy();
		CollAr->Draw();*/


	//-------------------------------Events creation (Without CollAr)----------------------------------------


	double XX_step = 0.2;
	int n = ceil((log10(Emax)-log10(Emin))/XX_step);
	double *XX = new double[n+1];
	for(int i=0;i<n+1;i++)
	{
		XX[i]= pow(10,(log10(Emin)+i*XX_step));
	}

	double XX2_step = 0.2;
	int n2 = ceil((log10(Emax)-(log10(Emin)-0.6))/XX2_step);
	double *XX2 = new double[n2+1];
	for(int i=0;i<n2+1;i++)
	{
		XX2[i]= pow(10,((log10(Emin)-0.6)+i*XX2_step));
	}



	Double_t E[numev], t[numev];

	TH1D* energyflare = new TH1D("energyflare", "energyflare", n, XX);
	TH1D* energybaseline = new TH1D("energybaseline", "energybaseline", n, XX);

	TH1D* time = new TH1D("time", "time", 30, tmin, tmax);


	Double_t norm, mean, sigma, LIV;
	norm=1.98;
	mean=2009.;
	sigma=686.9;
	LIV=0.;

	/*TF2* simpler = new TF2("simpler","TMath::Power(y,-2.7)+[0]*TMath::Power(y,-2.4)*exp(-0.5*TMath::Power(((x-[1]-[3]*y)/[2]),2))"
		, tmin, tmax, Emin, Emax);
	simpler->SetParameters(norm, mean, sigma, LIV);
	simpler->SetNpx(2000);
	simpler->SetNpy(2000);*/


	TF1* spectrumfit = new TF1("spectrumfit", powerLawFit, Emin, 10., 2);

	/*TCanvas* can = new TCanvas();
	simpler->Draw("surf1");*/

	/*gRandom->SetSeed(0);
	for(Int_t i=0; i<numev; i++)
	{
		simpler->GetRandom2(t[i], E[i]);
		if(t[i]<1200)energybaseline->Fill(E[i]);
		else energyflare->Fill(E[i], 1/E[i]);
		time->Fill(t[i]);
	}*/

	/*TCanvas* canvas1 = new TCanvas();
	canvas1->SetLogx();
	canvas1->SetLogy();
	energyflare->SetTitle("Energy distribution (Flare)");
	energyflare->GetXaxis()->SetTitle("Energy (TeV)");
	energyflare->Fit(spectrumfit, "QR");
	energyflare->Draw();
	TCanvas* canvas2 = new TCanvas();
	time->SetTitle("Time distribution");
	time->GetXaxis()->SetTitle("Time (s)");
	time->Fit("gaus");
	time->Draw();*/


	//-------------------------------Events creation (With CollAr)----------------------------------------


	TF2* complex = new TF2("complex", signal2D, tmin, tmax, Emin, Emax, 4);
	complex->SetParameters(norm, mean, sigma, LIV);
	complex->SetNpx(2000);
	complex->SetNpy(2000);

	/*TCanvas* can = new TCanvas();
	complex->Draw("surf1");*/

	gRandom->SetSeed(0);
	for(Int_t i=0; i<numev; i++)
	{
		complex->GetRandom2(t[i], E[i]);
		//if(t[i]<1200)energybaseline->Fill(E[i], 1/E[i]);
		energyflare->Fill(E[i], 1/E[i]);//else
		time->Fill(t[i]);
	}


	/*TCanvas* canvas1 = new TCanvas();
	canvas1->SetLogx();
	canvas1->SetLogy();*/
	//energyflare->SetTitle("Energy distribution (Flare+Baseline)");
	//energyflare->GetXaxis()->SetTitle("Energy (TeV)");
	//energybaseline->Fit(spectrumfit, "QR");
	//energyflare->Draw("EP");
	/*TCanvas* canvas2 = new TCanvas();
	time->SetTitle("Time distribution");
	time->GetXaxis()->SetTitle("Time (s)");
	time->Fit("gaus");
	time->Draw();*/

	//-------------------------------Smearing of the events----------------------------------------


	TH1D* energyrec = new TH1D("energyrec", "energyrec", n2, XX2);

	Double_t res = 0.215;
	cout << "Starting smearing of events..." << endl;
	cout << "Resolution " << res << endl;

	Double_t Esmeared[numev];

	TF1 *gausian = new TF1("gausian", "gaus", Emin-2., Emax+2.);

	for(Int_t a=0; a<numev; a++)
	{

		gausian->SetParameters(1, E[a], res*E[a]);
		Esmeared[a] = gausian->GetRandom();
		energyrec->Fill(Esmeared[a], 1/Esmeared[a]);
		//cout << "True energy (TeV) " << E[a] << endl;
		//cout << "Reconstructed energy(TeV) " << Esmeared[a] << endl;

	}

	/*energyrec->SetLineColor(2);
	energyrec->SetTitle("Energy distribution (Flare+Baseline)");
	energyrec->GetXaxis()->SetTitle("Energy (TeV)");
	energyrec->Draw("EP");
	energyflare->Draw("sameEP");*/

	/*TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(energyflare,"Simulated events","le");
	leg1->AddEntry(energyrec,"Smeared events","le");
	leg1->Draw();*/


	//-------------------------------Save the events in a file----------------------------------------

	//Write a txt file with the evetns
	Int_t num = 4;
	FILE *fout;
	fout = fopen(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Mrk501/2D/tau0/flaresme%i.txt", num),"wb");

	Int_t total=0;
	for(int i = 0; i<numev; i++)
	{
		if(Esmeared[i]>0.15)
		{
			fprintf(fout,"%0.3lf    %0.8lf   \n", t[i], Esmeared[i]);
			total++;

		}
	}
	cout << "Final number of events: " << total << endl;
}
