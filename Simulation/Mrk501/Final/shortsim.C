/*
 * shortsim.C
 *
 *  Created on: Mar 25, 2016
 *      Author: lnogues
 *
 *      In this simulation, the effect of the acceptance is applied created a new
 *      template for the energy (the weighted power law), from where the energy events
 *      are sampled, so the weights are not necessary)
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


const Int_t simnum = 300;
Double_t tmean[simnum], chigaus[simnum], ind[simnum], chipl[simnum], Nlc[simnum], Nsp[simnum], Emean[simnum], TimeSigma[simnum], TimeMean[simnum];

Double_t powerLaw(Double_t* x, Double_t* par)
{
  Double_t func = par[0]*TMath::Power(x[0]/par[1],-par[2]);
  return func;
}

Double_t powerLawFit(Double_t* x, Double_t* par)
{
  Double_t func2 = par[0]*TMath::Power(x[0],-par[1]);
  return func2;
}

Double_t gaus(Double_t* x, Double_t* par)
{
	Double_t funci = 1.*exp(-(pow(x[0]-par[0],2)/(2*pow(par[1],2))));
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
	 gStyle->SetTitleSize(0.45);
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

	const int numev = 1490;

	//------------------ --------------------Time generation------------------------------------------

	Double_t tmax, sigma, tinit, tfin;
	tmax = 805;
	sigma = 219;
	tinit = 0;
	tfin = 1531;

	TF1 *lightcurve = new TF1("lightcurve", gaus, tinit, tfin, 2);
	lightcurve->SetParameters(tmax, sigma);
	TH1D* timehisto = new TH1D("timehisto", "timehisto", 50, tinit, tfin);

	gRandom->SetSeed(0);
	Double_t t[numev] = 0;
	for(Int_t i = 0; i<numev; i++)
	{
		t[i]=lightcurve->GetRandom(tinit, tfin);
		timehisto->Fill(t[i]);
	}

	//--------------------------------------Energy generation------------------------------------------

	Double_t Emin,Emax;
	Emin = 0.15;
	Emax = 11.210;

	Double_t norm, index, E0;
	norm = 14.3e-6;
	index = 2.2;
	E0 = 0.3;

	TF1* spectrum = new TF1("spectrum", powerLaw, Emin, Emax, 3);
	spectrum->SetParameters(norm, E0,index);

	TF1* spectrumfit = new TF1("spectrumfit", powerLawFit, 0.3, 8., 2);
	spectrumfit->SetName("spectrumfit");
	spectrumfit->SetParName(0, "norm");
	spectrumfit->SetParName(1, "index");

	double XX_step = 0.1;
	int n = ceil((log10(Emax)-log10(Emin))/XX_step);
	double *XX = new double[n+1];
	for(int i=0;i<n+1;i++)
	{
		XX[i]= pow(10,(log10(Emin)+i*XX_step));
	}


	TH1D* energyhisto = new TH1D("energyhisto", "energyhisto", n, XX);
	energyhisto->Sumw2();

	Double_t E[numev] = 0;
	for(Int_t i = 0; i<numev; i++)
	{
		E[i]=spectrum->GetRandom(Emin, Emax);
		energyhisto->Fill(E[i], 1./E[i]);
	}


	//--------------------------------------Collection Area------------------------------------------

	TString myPath = "/home/lnogues/workspace_cpp/LIVanalysis/Simulation/";
	const TString sNewFile = myPath+"MonoCollAr.root";
	TFile* fNewFile = new TFile(sNewFile);
	collar = (TGraph*)fNewFile->Get("Mono");

	Double_t fmax = 96409.73;

	Double_t factor[1490] = 0;
	Double_t factor2[1490] = 0;
	Double_t total = 0;
	for(int i = 0; i<1490; i++)
	{
		factor[i] = collar->Eval(E[i])/fmax;
		total += factor[i];
	}

	Double_t meanfactor = total/1490;

	TH1D* weighthisto = new TH1D("weighthisto", "Energy template",  n, XX);
	for(Int_t i = 0; i<numev; i++)
	{
		factor2[i] = factor[i]/meanfactor;
		weighthisto->Fill(E[i], factor2[i]/E[i]);
	}


	//--------------------------------------Mono Migration Matrix-------------------------------------

	TFile* fileStatusFluxlc = TFile::Open("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Status_Output_fluxlc.root","READ");
	MStatusArray* StatusOutputFluxlc = fileStatusFluxlc->Get("MStatusDisplay");
	TCanvas * MigMatrix = (TCanvas *) StatusOutputFluxlc->FindCanvas("Migration Matrix");
	TH2D* mm = (TH2D*)MigMatrix->GetPrimitive("mimgplot");

	Double_t oldmean[17], oldSigma[17], oldEtrue[17], oldbias[17], oldsigmaE[17] = 0;
	Int_t oldbin;
	Int_t initbin = 5;
	TH1D* oldprojX;
	TF1* oldmyfunc;

	for(Int_t i = initbin; i<=18; i++)
	{
		oldbin = i;
		oldprojX = mm->ProjectionX("oldprojX", oldbin, oldbin, "");
		oldprojX->Fit("gaus", "0Q");

		oldmyfunc = oldprojX->GetFunction("gaus");
		Double_t oldpar1 = oldmyfunc->GetParameter(0);
		Double_t oldpar2 = oldmyfunc->GetParameter(1);
		Double_t oldpar3 = oldmyfunc->GetParameter(2);

		oldEtrue[oldbin-initbin] = mm->GetYaxis()->GetBinCenter(oldbin);
		oldmean[oldbin-initbin] = oldpar2;
		oldSigma[oldbin-initbin] = oldpar3;
		oldbias[oldbin-initbin] = oldmean[oldbin-initbin]/oldEtrue[oldbin-initbin];
		oldsigmaE[oldbin-initbin] = oldSigma[oldbin-initbin]/oldEtrue[oldbin-initbin];

	}

	Int_t count = 0;
	Double_t oldadd = 0;
	for(Int_t point = 0; point <=oldbin-initbin; point++)
	{
		if(oldEtrue[point]>300 &&oldEtrue[point]<8000)
			{
				oldadd+= oldsigmaE[point];
				count++;
			}
	}

	Double_t res = oldadd/count;

	Double_t biasLC, biasFit = 0;
	Int_t iLC, iFit = 0;
	for(Int_t point = 0; point <=oldbin-initbin; point++)
	{
		if(oldEtrue[point]< 300)
		{
			iLC++;
			biasLC += oldbias[point];
		}
		if(oldEtrue[point]>=300 && oldEtrue[point]<= Emax)
		{
			iFit++;
			biasFit += oldbias[point];
		}
	}

	Double_t finalbiasLC= biasLC/iLC;
	Double_t finalbiasFit= biasFit/iFit;



	//-------------------------------------Flare generation-----------------------------------

	//Energy events simulation (using the weighted power law)
	Double_t Ew[numev] = 0;
	for(int i = 0; i<numev; i++)
	{
		Ew[i] = weighthisto->GetRandom();
	}

	//Lag injection
	Double_t tau = -80;
	Double_t tlag[numev] = 0;

	TH1D* laghisto = new TH1D("laghisto", "laghisto", 50, tinit, tfin);
	for(int i = 0; i<numev; i++)
	{
		tlag[i]=t[i]+tau*Ew[i];
		laghisto->Fill(tlag[i]);
	}

	//Smear of the events
	Double_t Esmeared[numev] = 0;

	double XX2_step = 0.1;
	int n2 = ceil((log10(Emax)-(log10(Emin)-0.6))/XX2_step);
	double *XX2 = new double[n2+1];
	for(int i=0;i<n2+1;i++)
	{
		XX2[i]= pow(10,((log10(Emin)-0.6)+i*XX2_step));
	}

	TF1 *gausian = new TF1("gausian", "gaus", Emin-2., Emax+2.);
	for(Int_t a=0; a<numev; a++)
	{
		gausian->SetParameters(1, Ew[a], res*Ew[a]);
		Esmeared[a] = gausian->GetRandom();
	}

	TH1D* smearhisto = new TH1D("smearhisto", "smearhisto", n2, XX2);
	for(int i = 0; i<numev; i++)
	{
		smearhisto->Fill(Esmeared[i], 1/Esmeared[i]);
	}

	smearhisto->Fit(spectrumfit, "QR");
	TF1* fitpl = smearhisto->GetFunction("spectrumfit");
	ind[num] = fitpl->GetParameter(1);
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
		if(Esmeared[i] >= 0.15 && Esmeared[i] < 0.3)
		{
			LCtemplater->Fill(tlag[i]);
			LCphotons++;
		}
		if(Esmeared[i] >= 0.3 && Esmeared[i] <= 8.)
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

	Emean[num] = Esum/Fitphotons;

	LCtemplater->Sumw2();
	LCtemplater->Scale(1./LCtemplater->Integral());
	LCtemplater->Fit("gaus", "Q");
	TF1* fitlc = LCtemplater->GetFunction("gaus");
	tmean[num] = fitlc->GetParameter(1);
	chigaus[num] = (fitlc->GetChisquare())/(fitlc->GetNDF());

	timeevents->Sumw2();
	timeevents->Scale(1/timeevents->Integral());
	timeevents->Fit("gaus", "Q");
	TF1* fitevt = timeevents->GetFunction("gaus");
	TimeMean[num] = fitevt->GetParameter(1);
	TimeSigma[num] = fitevt->GetParameter(2);


	Nlc[num] = LCphotons;
	Nsp[num] = Fitphotons;

	//Write a txt file with the evetns
	FILE *fout;
	fout = fopen(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Final/tau-80/flarer%i.txt", num),"wb");
	fprintf(fout,"%0.3lf \n", LCphotons);
	fprintf(fout,"%0.3lf \n", Fitphotons);
	fprintf(fout,"%0.3lf \n", chipl[num]);
	for(int i = 0; i<i3; i++)
	{
		fprintf(fout,"%0.3lf    %0.8lf   \n", tsimr[i], Esimr[i]);
	}

	//Write a root file with the important plots
	TString rootoutfile(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Final/tau-80/Simulationr%i.root", num));
	TFile* outputfile = new TFile(rootoutfile, "recreate");
	timehisto->SetTitle("Time generated events");
	timehisto->Write();
	energyhisto->SetTitle("Raw energy generated events");
	energyhisto->Write();
	weighthisto->SetTitle("Weigted events: Energy template");
	weighthisto->Write();
	laghisto->SetTitle(Form("Events with lag %f", tau));
	laghisto->Write();
	smearhisto->SetTitle("Smeared events");
	smearhisto->Write();
	LCtemplater->SetTitle("LC template: 0.15-0.3 TeV");
	LCtemplater->SetName("LCtemplate");
	LCtemplater->Write();
	timeevents->SetTitle("Flare time events: 0.3-8 TeV");
	timeevents->Write();
	energyevents->SetTitle("Flare energy events: 0.3-8 TeV");
	energyevents->Write();
	outputfile->Close();


	delete timehisto;
	delete energyhisto;
	delete weighthisto;
	delete laghisto;

}


void shortsim()
{
	Int_t iter;

	for(Int_t i = 0; i < simnum; i++)
	{
		iter = i;
		cout << "----------------------------- SIMULATION " << iter << " --------------------------" << endl;
		simulation(iter);
		/*cout << "tmean " <<tmean[iter] << endl;
		cout << "chigaus " << chigaus[iter] << endl;
		cout << "ind " << ind[iter] << endl;
		cout << "chipl " << chipl[iter] << endl;
		cout << "Nlc " << Nlc[iter] << endl;
		cout << "Nsp " << Nsp[iter] << endl;*/
	}

	//Mean computation
	Double_t sumtmean, sumtime, sumE, sumpl;
	sumtmean = 0;
	sumtime = 0;
	sumE = 0;
	sumpl = 0;
	for(Int_t i=0; i<simnum; i++)
	{
		sumtmean+=tmean[i];
		sumtime+=TimeMean[i];
		sumE+=Emean[i];
		sumpl+=ind[i];
	}

	Double_t meantmean, meantime, meanE, meanpl;
	meantmean=sumtmean/simnum;
	meantime=sumtime/simnum;
	meanE=sumE/simnum;
	meanpl=sumpl/simnum;


	TH1D* tmeanhisto = new TH1D("tmeanhisto", "tmeanhisto", 50, meantmean-200, meantmean+200);
	TH1D* chigaushisto = new TH1D("chigaushisto", "chigaushisto", 50, 0, 2);
	TH1D* indhisto = new TH1D("indhisto", "indhisto", 50, meanpl-0.5, meanpl+0.8);
	TH1D* chiplhisto = new TH1D("chiplhisto", "chiplhisto", 50, 0, 5);
	TH1D* Nlchisto = new TH1D("Nlchisto", "Nlchisto", 50, 500, 900);
	TH1D* Nsphisto = new TH1D("Nsphisto", "Nsphisto", 50, 300, 800);
	TH1D* TimeMeanhisto = new TH1D("TimeMeanhisto", "TimeMeanhisto", 50, meantime-200, meantime+200);
	TH1D* MeanEhisto = new TH1D("MeanEhisto", "MeanEhisto", 50, meanE-0.3, meanE+0.3);

	for(Int_t i = 0; i<=iter; i++)
	{
		tmeanhisto->Fill(tmean[i]);
		chigaushisto->Fill(chigaus[i]);
		indhisto->Fill(ind[i]);
		chiplhisto->Fill(chipl[i]);
		Nlchisto->Fill(Nlc[i]);
		Nsphisto->Fill(Nsp[i]);
		TimeMeanhisto->Fill(TimeMean[i]);
		MeanEhisto->Fill(Emean[i]);
	}

	gStyle->SetOptStat(111111);

	TCanvas* big = new TCanvas("tau-80simr", "tau-80simr", 1000, 1000);
	big->Divide(2,4);
	big->cd(1);
	gPad->SetTitle("LC Tmax");
	tmeanhisto->SetName("LC Tmax");
	tmeanhisto->GetXaxis()->SetTitle("s");
	tmeanhisto->Draw("EP");
	big->cd(2);
	gPad->SetTitle("Chi2/NDF(LC)");
	chigaushisto->SetName("#Chi^{2}/NDF(LC)");
	chigaushisto->GetXaxis()->SetTitle("#Chi^{2}/NDF");
	chigaushisto->Draw("EP");
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
	big->SaveAs("./tau-80/Canvastau-80r.png");

}

