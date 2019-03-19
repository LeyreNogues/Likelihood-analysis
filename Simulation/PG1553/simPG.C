/*
 * shortsim.C
 *
 *  Created on: Mar 25, 2016
 *      Author: lnogues
 *
 *  In this simulation, the effect of the acceptance is applied computing weights for
 *  every event. This weights would be used in the Likelihood
 */

#include <TF1.h>
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
#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>
#include <iostream>

using namespace std;

const Int_t simnum = 100;
Double_t tinit = 0.;
Double_t tfin = 8000.;

Double_t tmean[simnum], ind[simnum], chipl[simnum], Nlc[simnum], Nsp[simnum], Emean[simnum], TimeMean[simnum], SBratio[simnum];

Double_t powerLaw(Double_t* x, Double_t* par)
{
	Double_t func = TMath::Power(x[0],-par[0]);
	 return func;
}

Double_t doblegaus(Double_t* x, Double_t* par)
{

	Double_t bkg, bkg2, tmax, sigma,tmax2, sigma2, tinit, tfin;
	bkg = 18.32;
	tmax = 2540.;
	sigma = 798;
	bkg2 = 15.62;
	tmax2 = 5898.;
	sigma2 = 1070.;

	TF1 *gaus1 = new TF1("gaus1", gaus, tinit, tfin);
	gaus1->SetParameters(bkg, tmax, sigma);
	TF1 *gaus2 = new TF1("gaus2", gaus, tinit, tfin);
	gaus2->SetParameters(bkg2, tmax2, sigma2);

	Double_t returnValue = gaus1->Eval(x[0])+gaus2->Eval(x[0]);
	delete gaus1;
	delete gaus2;
	return returnValue;
}

Double_t gaus(Double_t* x, Double_t* par) //with background
{
	Double_t funci = par[0]*exp(-(pow(x[0]-par[1],2)/(2*pow(par[2],2))));
	return funci;
}

void simulation(Int_t num)
{
	//settings AJ
	//gStyle->SetOptStat(0);
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
	gStyle->SetTitleSize(0.40);
	gStyle->SetTextSize(0.9);
	gStyle->SetTitleFontSize(2);
	gStyle->SetTitleFont(42);
	gStyle->SetTitleFont(62,"xyz");
	gStyle->SetTitleYOffset(1.0);
	gStyle->SetTitleXOffset(1.0);
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetTitleX(.15);
	gStyle->SetTitleY(.98);
	gStyle->SetTitleW(.70);
	gStyle->SetTitleH(.05);
	// statistics box
	gStyle->SetStatFont(42);
	gStyle->SetStatX(.91);
	gStyle->SetStatY(.90);
	gStyle->SetStatW(.23);
	gStyle->SetStatH(.18);
	// axis labels
	gStyle->SetLabelFont(42,"xyz");
	gStyle->SetLabelSize(0.056,"xyz");
	gStyle->SetGridColor(16);
	gStyle->SetLegendBorderSize(0);
	gStyle->SetOptStat(111111);
	gStyle->SetOptFit(111111);

	const int numev = 950;

	//------------------ --------------------Time generation------------------------------------------

	/*Double_t bkg, bkg2, tmax, sigma,tmax2, sigma2, tinit, tfin;
	bkg = 18.32;
	tmax = 2540.;
	sigma = 798;
	bkg2 = 15.62;
	tmax2 = 5898.;
	sigma2 = 1070.;*/

	TF1 *lightcurve = new TF1("lightcurve",doblegaus, tinit, tfin, 0); //FOr the signal
	/*TF1 *lightcurve2 = new TF1("lightcurve2","[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[4])/[5])**2)", tinit, tfin, 6);
	lightcurve2->SetParameters(bkg, tmax, sigma, bkg2, tmax2, sigma2);*/

	TH1D* signaltimehisto = new TH1D("time events", "time events", 50, tinit, tfin);
	TH1D* bkgtimehisto = new TH1D("bkgtimehisto", "bkgtimehisto", 50, tinit, tfin);
	TH1D* timehisto = new TH1D("timehisto", "timehisto", 50, tinit, tfin);


	TRandom* bkg = new TRandom();
	gRandom->SetSeed(0);
	Double_t t[numev] = 0;
	for(Int_t i = 0; i<numev; i++)
	{
		if(i < numev-220) //Only signal
		{
			t[i]=lightcurve->GetRandom(tinit, tfin);
			signaltimehisto->Fill(t[i]);
		}
		else //Only bkg
		{
			t[i]=bkg->Uniform(tinit, tfin);
			bkgtimehisto->Fill(t[i]);
		}
		timehisto->Fill(t[i]); //all the events (signal+bkg)
	}

	/*TCanvas* canvas = new TCanvas();
	lightcurve->SetTitle("Time Distribution");
	lightcurve->GetXaxis()->SetTitle("time(s)");*/
	//signaltimehisto->SetLineColor(2);
	//signaltimehisto->Draw("EP");
	//bkgtimehisto->SetLineColor(2);
	//bkgtimehisto->Draw("EP");
	//timehisto->SetLineColor(2);
	//timehisto->Draw("EP");


	//----------------------------------Energy generation (Bkg and signal)------------------------------------------

	Double_t Emin,Emax;
	Emin = 0.2;
	Emax = 1.2;

	Double_t indexs, indexb;
	indexs = 4.8;
	indexb = 2.5;

	TF1* signalspectrum = new TF1("signalspectrum", powerLaw, Emin, Emax, 1);
	signalspectrum->SetParameter(0, indexs);
	TF1* bkgspectrum = new TF1("bkgspectrum", powerLaw, Emin, Emax, 1);
	bkgspectrum->SetParameter(0, indexb);
	TF1* fitspectrum = new TF1("fitspectrum", powerLaw, 0.4, 1., 1);
	fitspectrum->SetName("fitspectrum");

	double XX_step = 0.1;
	int n = ceil((log10(Emax)-log10(Emin))/XX_step);
	double *XX = new double[n+1];
	for(int i=0;i<n+1;i++)
	{
		XX[i]= pow(10,(log10(Emin)+i*XX_step));
	}

	TH1D* energyhisto = new TH1D("energy events", "energy events", n, XX);
	TH1D* signalenergyhisto = new TH1D("signalenergyhisto", "signalenergyhisto", n, XX);
	TH1D* bkgenergyhisto = new TH1D("bkgenergyhisto", "bkgenergyhisto", n, XX);
	energyhisto->Sumw2();

	Double_t E[numev] = 0;
	for(Int_t i = 0; i<numev; i++)
	{
		if(i<numev-220)
		{
			E[i]=signalspectrum->GetRandom(Emin, Emax);
			signalenergyhisto->Fill(E[i], 1./E[i]);
		}
		else
		{
			E[i]=bkgspectrum->GetRandom(Emin, Emax);
			bkgenergyhisto->Fill(E[i], 1./E[i]);
		}
		energyhisto->Fill(E[i], 1./E[i]);
	}

	/*TCanvas* canvas2 = new TCanvas();
	canvas2->SetLogx();
	canvas2->SetLogy();
	signalenergyhisto->SetLineColor(2);
	signalenergyhisto->Draw("EP");
	bkgenergyhisto->Draw("sameEP");*/
	/*signalspectrum->SetLineColor(4);
	signalspectrum->GetXaxis()->SetTitle("Energy(TeV)");
	signalspectrum->Draw();
	bkgspectrum->SetLineColor(2);
	bkgspectrum->Draw("same");
	energyhisto->GetXaxis()->SetTitle("Energy(TeV)");*/
	//energyhisto->Draw();

	/*TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(signalspectrum,"Signal","l");
	leg1->AddEntry(bkgspectrum,"Background","l");
	leg1->Draw();*/

	//-------------------------------------Flare generation-----------------------------------


	//Lag injection (Only in the signal)
	Double_t tau =-2000;
	Double_t tlag[numev] = 0;

	TH1D* laghisto = new TH1D("laghisto", "laghisto", 50, tinit, tfin);
	TH1D* signallaghisto = new TH1D("signallaghisto", "signallaghisto", 50, tinit, tfin);
	TH1D* bkglaghisto = new TH1D("bkglaghisto", "bkglaghisto", 50, tinit, tfin);
	for(int i = 0; i<numev; i++)
	{
		if(i<numev-220)
			{
			tlag[i]=t[i]+tau*E[i];
			signallaghisto->Fill(tlag[i]);
			}
		else
			{
			tlag[i]=t[i];
			bkglaghisto->Fill(tlag[i]);
			}
		laghisto->Fill(tlag[i]);
	}

	//signallaghisto->SetLineColor(4);
	//signallaghisto->Draw("sameEP");
	//bkglaghisto->SetLineColor(4);
	//bkglaghisto->Draw("sameEP");
	//laghisto->SetLineColor(4);
	//laghisto->Draw("same EP");

	//Smear of the events
	Double_t Esmeared[numev] = 0;

	double XX2_step = 0.1;
	int n2 = ceil((log10(Emax)-(log10(Emin)-0.1))/XX2_step);
	double *XX2 = new double[n2+1];
	for(int i=0;i<n2+1;i++)
	{
		XX2[i]= pow(10,((log10(Emin)-0.1)+i*XX2_step));
	}

	TF1 *gausian = new TF1("gausian", "gaus", Emin-2., Emax+2.);
	TH1D* smearhisto = new TH1D("smearhisto", "smearhisto", n2, XX2);
	for(Int_t a=0; a<numev; a++)
	{
		gausian->SetParameters(1, E[a], 0.1*E[a]);
		Esmeared[a] = gausian->GetRandom();
		smearhisto->Fill(Esmeared[a], 1./Esmeared[a]);
	}

	/*TCanvas* canvas3 = new TCanvas();
	canvas3->SetLogx();
	canvas3->SetLogy();
	smearhisto->GetXaxis()->SetTitle("Energy(TeV)");
	smearhisto->SetLineColor(2);
	smearhisto->Draw("EP");
	energyhisto->Draw("sameEP");
	TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(energyhisto,"True","le");
	leg1->AddEntry(smearhisto,"Smeared","le");
	leg1->Draw();*/

	//Number of signal and bkg events at the end (Should be 317+144=461)
	Int_t hadrons, gammas;
	hadrons=0;
	gammas=0;
	Double_t ratio;

	for(Int_t a=0; a<numev; a++)
		{
			if(Esmeared[a]>=0.25 && Esmeared[a]<=1.)
			{
				if(a<numev-220) gammas++;
				else hadrons++;
			}
		}
	ratio=gammas/(double)hadrons;
	//cout << "Over 250 there are: "  << gammas <<" gammas and " << hadrons << " hadrons." << endl;
	//cout << "S/B ratio: " << ratio << endl;
	SBratio[num]= ratio;



	//PL template ( Only signal in the fit energy region)
	double XX3_step = 0.08;
	int n3 = ceil((log10(1.)-log10(0.4))/XX3_step);
	double *XX3 = new double[n3+1];
	for(int i=0;i<n3+1;i++)
	{
	XX3[i]= pow(10,(log10(0.4)+i*XX3_step));
	}
	TH1D* plhisto = new TH1D("plhisto", "plhisto", n3, XX3);
	for(Int_t a=0; a<numev-220; a++) //Only signal
	{
		if(Esmeared[a]>=0.4 && Esmeared[a]<=1.)
		{
			plhisto->Fill(Esmeared[a], 1./Esmeared[a]);
		}
	}

	/*TCanvas* plcanvas = new TCanvas();
	plcanvas->SetLogx();
	plcanvas->SetLogy();
	plhisto->GetXaxis()->SetTitle("Energy(TeV)");*/
	plhisto->Fit(fitspectrum, "0QR");
	//plhisto->Draw("EP");
	TF1* fitpl = plhisto->GetFunction("fitspectrum");
	ind[num] = fitpl->GetParameter(0);
	chipl[num] = (fitpl->GetChisquare())/(fitpl->GetNDF());



	//Division in LC and Fit
	Int_t LCphotons = 0;
	Int_t Fitphotons = 0;
	TH1D* LCtemplater = new TH1D("LCtemplater", "LCtemplater", 50, tinit, tfin);
	TH1D* energyevents = new TH1D("energyevents", "energyevents", n, XX);
	TH1D* timeevents = new TH1D("timeevents", "timeevents", 50, tinit, tfin);
	Double_t Esum = 0;
	Double_t Esimr[numev] = 0;
	Double_t tsimr[numev] = 0;
	Int_t i3 = 0;

	for(int i = 0; i<numev; i++)
	{
		if(Esmeared[i] >= 0.3 && Esmeared[i] < 0.4)
		{
			LCtemplater->Fill(tlag[i]);
			LCphotons++;
		}
		if(Esmeared[i] >= 0.4 && Esmeared[i] <= 1.)
		{
			Esimr[i3] = Esmeared[i];
			tsimr[i3] = tlag[i];
			Esum+=Esimr[i3];
			timeevents->Fill(tsimr[i3]);
			energyevents->Fill(Esimr[i3]);
			Fitphotons++;
			i3++;
		}
	}

	//cout << "LC photons: " << LCphotons << endl;
	//cout << "Fit photons: " << Fitphotons << endl;

	Emean[num] = Esum/Fitphotons;
	tmean[num] = timehisto->GetMean();
	TimeMean[num] = timeevents->GetMean();

	Nlc[num] = LCphotons;
	Nsp[num] = Fitphotons;

	//Write a txt file with the evetns
	FILE *fout;
	fout = fopen(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau-2000/flarepg%i.txt", num),"wb");
	fprintf(fout,"%0.3lf \n", LCphotons);
	fprintf(fout,"%0.3lf \n", Fitphotons);
	fprintf(fout,"%0.3lf \n", chipl[num]);
	for(int i = 0; i<i3; i++)
	{
		fprintf(fout,"%0.3lf    %0.8lf\n", tsimr[i], Esimr[i]);
	}

	//Write a root file with the important plots
	TString rootoutfile(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau-2000/Simulationpg%i.root", num));
	TFile* outputfile = new TFile(rootoutfile, "recreate");
	signaltimehisto->SetTitle(" Signal Time generated events");
	signaltimehisto->Write();
	energyhisto->SetTitle("Raw energy generated events");
	energyhisto->Write();
	laghisto->SetTitle(Form("Events with lag %1.2f", tau));
	laghisto->Write();
	smearhisto->SetTitle("Smeared events");
	smearhisto->Write();
	plhisto->SetTitle("PL template: 0.4-1 TeV");
	plhisto->SetName("PLtemplate");
	plhisto->Write();
	timeevents->SetTitle("Flare time events: 0.3-8 TeV");
	timeevents->Write();
	energyevents->SetTitle("Flare energy events: 0.3-8 TeV");
	energyevents->Write();
	outputfile->Close();


	delete signaltimehisto;
	delete bkgtimehisto;
	delete energyhisto;
	delete signalenergyhisto;
	delete bkgenergyhisto;
	delete laghisto;
	delete signallaghisto;
	delete bkglaghisto;
	delete smearhisto;
	delete plhisto;
	delete LCtemplater;
	delete energyevents;
	delete timeevents;
	delete timehisto;
	delete signaltimehisto;
	delete bkgtimehisto;

}


void simPG()
{
	Int_t iter;

	for(Int_t i = 0; i < simnum; i++)
	{
		iter = i;
		cout << "----------------------------- SIMULATION " << iter << " --------------------------" << endl;
		simulation(iter);
		/*cout << "tmean " <<tmean[iter] << endl;
		cout << "ind " << ind[iter] << endl;
		cout << "chipl " << chipl[iter] << endl;
		cout << "Nlc " << Nlc[iter] << endl;
		cout << "Nsp " << Nsp[iter] << endl;
		cout << "tevents " << TimeMean[iter] << endl;
		cout << "Eevents " << Emean[iter] << endl;*/
	}

	//Mean computation
	Double_t sumtmean, sumtime, sumE, sumpl, sumratio, sumLC, sumPL;
	sumtmean = 0;
	sumtime = 0;
	sumE = 0;
	sumpl = 0;
	sumratio = 0;
	sumLC = 0;
	sumPL = 0;
	for(Int_t i=0; i<simnum; i++)
	{
		sumtmean+=tmean[i];
		sumtime+=TimeMean[i];
		sumE+=Emean[i];
		sumpl+=ind[i];
		sumratio+=SBratio[i];
		sumLC+=Nlc[i];
		sumPL+=Nsp[i];
	}

	Double_t meantmean, meantime, meanE, meanpl, meanratio, meanLC, meanPL;
	meantmean=sumtmean/simnum;
	meantime=sumtime/simnum;
	meanE=sumE/simnum;
	meanpl=sumpl/simnum;
	meanratio=sumratio/simnum;
	meanLC=sumLC/simnum;
	meanPL=sumPL/simnum;


	TH1D* tmeanhisto = new TH1D("tmeanhisto", "tmeanhisto", 50, meantmean-200, meantmean+200);
	TH1D* indhisto = new TH1D("indhisto", "indhisto", 50, meanpl-0.6, meanpl+0.8);
	TH1D* chiplhisto = new TH1D("chiplhisto", "chiplhisto", 50, 0, 5);
	TH1D* Nlchisto = new TH1D("Nlchisto", "Nlchisto", 50, meanLC-100, meanLC+100);
	TH1D* Nsphisto = new TH1D("Nsphisto", "Nsphisto", 50, meanPL-100, meanPL+100);
	TH1D* TimeMeanhisto = new TH1D("TimeMeanhisto", "TimeMeanhisto", 50, meantime-300, meantime+300);
	TH1D* MeanEhisto = new TH1D("MeanEhisto", "MeanEhisto", 50, meanE-0.3, meanE+0.3);
	TH1D* ratiohisto = new TH1D("ratiohisto", "ratiohisto", 50, meanratio-1, meanratio+1);


	for(Int_t i = 0; i<=iter; i++)
	{
		tmeanhisto->Fill(tmean[i]);
		indhisto->Fill(ind[i]);
		chiplhisto->Fill(chipl[i]);
		Nlchisto->Fill(Nlc[i]);
		Nsphisto->Fill(Nsp[i]);
		TimeMeanhisto->Fill(TimeMean[i]);
		MeanEhisto->Fill(Emean[i]);
		ratiohisto->Fill(SBratio[i]);
	}

	gStyle->SetOptStat(1111);

	TCanvas* big = new TCanvas("tau-2000simpg", "tau-2000simpg", 1000, 1000);
	big->Divide(2,4);
	big->cd(1);
	gPad->SetTitle("LC Tmax");
	tmeanhisto->SetName("LC Tmax");
	tmeanhisto->GetXaxis()->SetTitle("s");
	tmeanhisto->Draw("EP");
	big->cd(2);
	gPad->SetTitle("S/B ratio");
	ratiohisto->SetName("S/B ratio");
	ratiohisto->GetXaxis()->SetTitle("ratio(0.25-1TeV)");
	ratiohisto->Draw("EP");
	big->cd(3);
	gPad->SetTitle("PL index");
	indhisto->SetName("PL index");
	indhisto->GetXaxis()->SetTitle("index");
	indhisto->Fit("gaus", "Q");
	indhisto->Draw("EP");
	big->cd(4);
	gPad->SetTitle("Chi2/NDF(PL)");
	chiplhisto->SetName("#Chi^{2}/NDF(PL)");
	chiplhisto->GetXaxis()->SetTitle("#Chi^{2}/NDF");
	chiplhisto->Draw("EP");
	big->cd(5);
	Nlchisto->SetTitle("Events LC");
	Nlchisto->SetName("Events LC");
	Nlchisto->GetXaxis()->SetTitle("events");
	Nlchisto->Draw("EP");
	big->cd(6);
	gPad->SetTitle("Events Fit");
	Nsphisto->SetName("Events Fit");
	Nsphisto->GetXaxis()->SetTitle("events");
	Nsphisto->Draw("EP");
	big->cd(7);
	gPad->SetTitle("Mean time Events");
	TimeMeanhisto->SetName("Mean time Events");
	TimeMeanhisto->GetXaxis()->SetTitle("s");
	TimeMeanhisto->Draw("EP");
	big->cd(8);
	gPad->SetTitle("Mean E Events");
	MeanEhisto->SetName("Mean E Events");
	MeanEhisto->GetXaxis()->SetTitle("TeV");
	MeanEhisto->Draw("EP");
	big->SaveAs("./tau-2000/Canvastau-2000pg.png");


}

