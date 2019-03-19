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

Double_t tmean[simnum], inds[simnum], indb[simnum], chipl[simnum], Nlc[simnum], Nsp[simnum], Emean[simnum], TimeMean[simnum], SBratio[simnum], Ecut[simnum];

Double_t powerLaw(Double_t* x, Double_t* par)
{
	Double_t func = TMath::Power(x[0],-par[0]);
	 return func;
}

Double_t powerLawfit(Double_t* x, Double_t* par)
{
	Double_t func = par[0]*TMath::Power(x[0],-par[1]);
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

	const int numev = 900;
	const int numbkg = 100;

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


	TRandom* bkg = new TRandom();
	gRandom->SetSeed(0);
	Double_t tsig[numev] = 0;
	Double_t tbkg[numbkg] = 0;

	for(Int_t i = 0; i<numev; i++) //Signal time events
	{
		tsig[i]=lightcurve->GetRandom(tinit, tfin);
		signaltimehisto->Fill(tsig[i]);
	}
	for(Int_t i = 0; i<numbkg; i++) //Background time events
	{
		tbkg[i]=bkg->Uniform(tinit, tfin);
		bkgtimehisto->Fill(tbkg[i]);
	}

	/*TCanvas* canvas = new TCanvas();
	signaltimehisto->SetTitle("Time Distribution");
	signaltimehisto->SetName("Time of events");
	signaltimehisto->GetXaxis()->SetTitle("time(s)");
	signaltimehisto->SetLineColor(2);
	signaltimehisto->Draw("EP");*/
	//bkgtimehisto->Draw("sameEP");

	/*TLegend* leg0 = new TLegend(0.1,0.7,0.4,0.9);
	leg0->AddEntry(signaltimehisto,"Signal","le");
	leg0->AddEntry(bkgtimehisto,"Background","le");
	leg0->Draw();*/

	//----------------------------------Energy generation (Bkg and signal)------------------------------------------

	Double_t Emin,Emax;
	Emin = 0.2;
	Emax = 1.2;

	Double_t indexs, indexb;
	indexs = 4.8;
	indexb = 4.5; //2.5

	TF1* signalspectrum = new TF1("signalspectrum", powerLaw, Emin, Emax, 1);
	signalspectrum->SetParameter(0, indexs);
	TF1* bkgspectrum = new TF1("bkgspectrum", powerLaw, Emin, Emax, 1);
	bkgspectrum->SetParameter(0, indexb);
	//Always two parameters to fit, dear
	TF1* fitspectrum = new TF1("fitspectrum", powerLawfit, 0.4, 1., 2);
	fitspectrum->SetName("fitspectrum");
	fitspectrum->SetParName(0, "norm");
	fitspectrum->SetParName(1, "index");

	double XX_step = 0.08;
	int n = ceil((log10(Emax)-log10(Emin))/XX_step);
	double *XX = new double[n+1];
	for(int i=0;i<n+1;i++)
	{
		XX[i]= pow(10,(log10(Emin)+i*XX_step));
	}

	TH1D* signalenergyhisto = new TH1D("signalenergyhisto", "signalenergyhisto", n, XX);
	TH1D* bkgenergyhisto = new TH1D("bkgenergyhisto", "bkgenergyhisto", n, XX);

	Double_t Esig[numev] = 0;
	Double_t Ebkg[numbkg] = 0;
	for(Int_t i = 0; i<numev; i++)
	{
		Esig[i]=signalspectrum->GetRandom(Emin, Emax);
		signalenergyhisto->Fill(Esig[i], 1./Esig[i]);
	}
	for(Int_t i = 0; i<numbkg; i++)
	{
		Ebkg[i]=bkgspectrum->GetRandom(Emin, Emax);
		bkgenergyhisto->Fill(Ebkg[i], 1./Ebkg[i]);
	}

	/*TCanvas* canvas2 = new TCanvas();
	canvas2->SetLogx();
	canvas2->SetLogy();*/
	/*signalenergyhisto->SetName("Energy of events");
	signalenergyhisto->SetLineColor(2);
	signalenergyhisto->Draw("EP");
	bkgenergyhisto->Draw("sameEP");*/
	/*signalspectrum->SetLineColor(4);
	signalspectrum->GetXaxis()->SetTitle("Energy(TeV)");
	signalspectrum->Draw();
	bkgspectrum->SetLineColor(2);
	bkgspectrum->Draw("same");*/

	/*TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(signalenergyhisto,"Signal","le");
	leg1->AddEntry(bkgenergyhisto,"Background","le");
	leg1->Draw();*/

	//-------------------------------------FLARE GENERATION-----------------------------------


	//Lag injection (Only in the signal)
	Double_t tau = 0;
	Double_t tlag[numev] = 0;

	TH1D* laghisto = new TH1D("laghisto", "laghisto", 50, tinit, tfin);

	for(int i = 0; i<numev; i++)
	{
			tlag[i]=tsig[i]+tau*Esig[i];
			laghisto->Fill(tlag[i]);
	}

	//laghisto->SetLineColor(4);
	//laghisto->Draw("sameEP");

	//Smear of the events
	Double_t Esmearedsig[numev] = 0;
	Double_t Esmearedbkg[numbkg] = 0;

	double XX2_step = 0.1;
	int n2 = ceil((log10(Emax)-(log10(Emin)-0.1))/XX2_step);
	double *XX2 = new double[n2+1];
	for(int i=0;i<n2+1;i++)
	{
		XX2[i]= pow(10,((log10(Emin)-0.1)+i*XX2_step));
	}

	TF1 *gausian = new TF1("gausian", "gaus", Emin-2., Emax+2.);
	TH1D* smearhistosignal = new TH1D("smearhistosignal", "smearhistosignal", n2, XX2);
	TH1D* smearhistobkg = new TH1D("smearhistobkg", "smearhistobkg", n2, XX2);

	for(Int_t a=0; a<numev; a++) //Smear of signal
	{
		gausian->SetParameters(1, Esig[a], 0.1*Esig[a]);
		Esmearedsig[a] = gausian->GetRandom();
		smearhistosignal->Fill(Esmearedsig[a], 1./Esmearedsig[a]);
	}
	for(Int_t a=0; a<numbkg; a++) //Smear of background
	{
		gausian->SetParameters(1, Ebkg[a], 0.1*Ebkg[a]);
		Esmearedbkg[a] = gausian->GetRandom();
		smearhistobkg->Fill(Esmearedbkg[a], 1./Esmearedbkg[a]);
	}

	/*TCanvas* canvas3 = new TCanvas();
	canvas3->SetLogx();
	canvas3->SetLogy();
	smearhistosignal->GetXaxis()->SetTitle("Energy(TeV)");
	smearhistosignal->SetLineColor(2);
	smearhistosignal->Draw("EP");
	signalenergyhisto->Draw("sameEP");*/
	/*smearhistobkg->GetXaxis()->SetTitle("Energy(TeV)");
	smearhistobkg->SetLineColor(2);
	smearhistobkg->Draw("EP");
	bkgenergyhisto->Draw("sameEP");
	TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(bkgenergyhisto,"True","le");
	leg1->AddEntry(smearhistobkg,"Smeared","le");
	leg1->Draw();*/

	//Number of  measured signal and bkg events (Should be 461)
	Int_t hadrons, gammas;
	hadrons=0;
	gammas=0;
	Double_t ratio;

	for(Int_t a=0; a<numev; a++)
	{
		if(Esmearedsig[a]>=0.25 && Esmearedsig[a]<=1.)
		{
			gammas++;
		}
	}
	for(Int_t a=0; a<numbkg; a++)
	{
		if(Esmearedbkg[a]>=0.25 && Esmearedbkg[a]<=1.)
		{
			hadrons++;
		}
	}

	//cout << "Measured events over 250 GeV are: " << gammas+hadrons << endl;


	//-------------------------------------PL template----------------------------------------

	//SIGNAL ONLY
	double XX3_step = 0.08;
	int n3 = ceil((log10(1.)-log10(0.4))/XX3_step);
	double *XX3 = new double[n3+1];
	for(int i=0;i<n3+1;i++)
	{
	XX3[i]= pow(10,(log10(0.4)+i*XX3_step));
	}
	TH1D* plhisto = new TH1D("plhistosig", "plhistosig", n3, XX3);
	for(Int_t a=0; a<numev; a++) //Only signal
	{
		if(Esmearedsig[a]>=0.4 && Esmearedsig[a]<=1.)
		{
			plhisto->Fill(Esmearedsig[a], 1./Esmearedsig[a]);
		}
	}

	//TCanvas* plcanvas = new TCanvas();
	//plcanvas->SetLogx();
	//plcanvas->SetLogy();
	//plhisto->GetXaxis()->SetTitle("Energy(TeV)");
	plhisto->Fit(fitspectrum, "0QR");
	//plhisto->Draw("EP");
	TF1* fitpl = plhisto->GetFunction("fitspectrum");
	inds[num] = fitpl->GetParameter(1);
	chipl[num] = (fitpl->GetChisquare())/(fitpl->GetNDF());

	//BACKGROUND
	/*TH1D* plhistobkg = new TH1D("plhistobkg", "plhistobkg", n3, XX3);
	for(Int_t a=0; a<numbkg; a++) //Only signal
		{
			if(Esmearedbkg[a]>=0.4 && Esmearedbkg[a]<=1.)
			{
				//plhisto->Fill(Esmearedbkg[a], 1./Esmearedbkg[a]);
				plhistobkg->Fill(Esmearedbkg[a], 1./Esmearedbkg[a]);
			}
		}*/

	/*TCanvas* plcanvas2 = new TCanvas();
	plcanvas2->SetTitle("signal");
	plcanvas2->SetLogx();
	plcanvas2->SetLogy();*/
	//plhisto->GetXaxis()->SetTitle("Energy(TeV)");
	//plhisto->Fit(fitspectrum, "0QR");
	//plhisto->Draw("EP");
	//TCanvas* plcanvas3 = new TCanvas();
	//plcanvas3->SetTitle("background");
	//plcanvas3->SetLogx();
	//plcanvas3->SetLogy();
	//plhistobkg->Fit(fitspectrum, "0QR");
	//plhistobkg->Draw("EP");
	//TF1* fitpl = plhisto->GetFunction("fitspectrum");
	//TF1* fitpl2 = plhistobkg->GetFunction("fitspectrum");
	//inds[num] = fitpl->GetParameter(1);
	//indb[num] = fitpl2->GetParameter(1);
	//chipl[num] = (fitpl->GetChisquare())/(fitpl->GetNDF());

	/*TCanvas* cutcanvas = new TCanvas();
	cutcanvas->SetLogx();
	cutcanvas->SetLogy();
	fitpl->SetLineColor(2);
	fitpl->Draw();
	fitpl2->SetLineColor(4);
	fitpl2->Draw("same");

	TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(fitpl,"Signal fit","l");
	leg1->AddEntry(fitpl2,"Bkg fit","l");
	leg1->Draw();*/

	/*Double_t par0sig=fitpl->GetParameter(0);
	Double_t par1sig=fitpl->GetParameter(1);
	Double_t par0bkg=fitpl2->GetParameter(0);
	Double_t par1bkg=fitpl2->GetParameter(1);

	Double_t cutE = TMath::Exp((TMath::Log(par0sig/par0bkg)/(par1sig-par1bkg)));
	if(cutE<100)Ecut[num] = cutE;
	else Ecut[num] = 1.;*/

	//cout << "Cut: "  << cutE << endl;


	//---------------------------------Division in LC and Fit--------------------------------

	Int_t LCphotonssig = 0;
	Int_t Fitphotonssig = 0;
	Int_t LCphotonsbkg = 0;
	Int_t Fitphotonsbkg = 0;
	TH1D* LCtemplate = new TH1D("LCtemplate", "LCtemplate", 50, tinit, tfin);
	TH1D* energyeventssig = new TH1D("energyeventssig", "energyeventssig", n, XX);
	TH1D* timeeventssig = new TH1D("timeeventssig", "timeeventssig", 50, tinit, tfin);
	TH1D* energyeventsbkg = new TH1D("energyeventsbkg", "energyeventsbkg", n, XX);
	TH1D* timeeventsbkg = new TH1D("timeeventsbkg", "timeeventsbkg", 50, tinit, tfin);
	Double_t Esimsig[numev] = 0;
	Double_t tsimsig[numev] = 0;
	Double_t Esimbkg[numbkg] = 0;
	Double_t tsimbkg[numbkg] = 0;
	Int_t sigcount, bkgcount;
	sigcount = 0; //Number of gammas in the fit region
	bkgcount = 0; //Number of hadrons in the fit region

	for(int i = 0; i<numev; i++)
	{
		if(Esmearedsig[i] >= 0.3 && Esmearedsig[i] < 0.4)
		{
			LCtemplate->Fill(tlag[i]);
			LCphotonssig++;
		}
		if(Esmearedsig[i] >= 0.4 && Esmearedsig[i] <= 1.)
		{
			Esimsig[sigcount] = Esmearedsig[i];
			tsimsig[sigcount] = tlag[i];
			timeeventssig->Fill(tsimsig[sigcount]);
			energyeventssig->Fill(Esimsig[sigcount]);
			Fitphotonssig++;
			sigcount++;
		}
	}

	//cout << "LC gammas: " << LCphotonssig << endl;
	//cout << "Fit gammas: " << Fitphotonssig << endl;

	for(int i = 0; i<numbkg; i++)
	{
		if(Esmearedbkg[i] >= 0.3 && Esmearedbkg[i] < 0.4)
		{
			LCphotonsbkg++;
		}
		if(Esmearedbkg[i] >= 0.4 && Esmearedbkg[i] <= 1.)
		{
			Esimbkg[bkgcount] = Esmearedbkg[i];
			tsimbkg[bkgcount] = tbkg[i];
			timeeventsbkg->Fill(tsimbkg[bkgcount]);
			energyeventsbkg->Fill(Esimbkg[bkgcount]);
			Fitphotonsbkg++;
			bkgcount++;
		}
	}

	//Esum

	//cout << "LC hadrons: " << LCphotonsbkg << endl;
	//cout << "Fit hadrons: " << Fitphotonsbkg << endl;
	ratio=Fitphotonssig/(double)Fitphotonsbkg;
	//cout << " Fit S/B ratio: " << ratio << endl;
	SBratio[num]= ratio;


	//------------------------------Put gammas and hadrons together for photon list--------------------------------

	//Put signal and backgroun alternatively untill we run out of backgrond.

	const Int_t f = Fitphotonssig+Fitphotonsbkg;
	Double_t Esim[f], tsim[f];
	Int_t nb, ns;
	nb = 0;
	ns = 0;

	//Look for odd/even condition to fill the array
	for(Int_t i=0; i<f; i++)
	{
		if(nb<Fitphotonsbkg)
		{
			if(i%2 == 0) //Even numbers
			{
				Esim[i]=Esimsig[ns];
				tsim[i]=tsimsig[ns];
				ns++;
				//cout << "lleno par" << endl;
			}
			if(i%2 != 0) //Odd numbers
			{
				Esim[i]=Esimbkg[nb];
				tsim[i]=tsimbkg[nb];
				nb++;
				//cout << "lleno impar" << endl;
			}
		}
		else
		{
			Esim[i]=Esimsig[ns];
			tsim[i]=tsimsig[ns];
			ns++;
			//cout << "lleno resto de signal" << endl;
		}

		//cout << "ns: " << ns << endl;
		//cout << "nb: " << nb << endl;

	}

	//cout << "Signal: " << Esimsig[ns-3] << " " << Esimsig[ns-2] << " " << Esimsig[ns-1] << endl;
	//cout << "Bkg: " << Esimbkg[nb-3] << " " << Esimbkg[nb-2] << " " << Esimbkg[nb-1] << endl;
	//cout << "Together: " << Esim[f-6] << " " << Esim[f-5] << " " << Esim[f-4] << " " << Esim[f-3] << " " << Esim[f-2] << " " << Esim[f-1] << endl;

	Double_t Esum = 0;
	//Mean energy
	for(Int_t i=0; i<f; i++)
	{
		Esum+=Esim[i];
	}
	Emean[num] = Esum/(double)f;
	//cout << "Total fit events: " << f << endl;
	//cout << "Emean: " << Emean[num] << endl;


	//-------------------------------------Photon list and graphs-----------------------------------




	tmean[num] = signaltimehisto->GetMean();
	TimeMean[num] = timeeventssig->GetMean();

	Nlc[num] = LCphotonssig+LCphotonsbkg;
	Nsp[num] = Fitphotonssig+Fitphotonsbkg;

	//cout << "Nlc: " << Nlc[num] << endl;
	//cout << "Nsp: " << Nsp[num] << endl;

	//Write a txt file with the evetns
	FILE *fout;
	fout = fopen(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau0/flarepgplbkg45%i.txt", num),"wb");
	fprintf(fout,"%0.3lf \n", LCphotonssig+LCphotonsbkg);
	fprintf(fout,"%0.3lf \n", Fitphotonssig+Fitphotonsbkg);
	//fprintf(fout,"%0.3lf \n", Fitphotonssig);
	//fprintf(fout,"%0.3lf \n", Fitphotonsbkg);
	fprintf(fout,"%0.3lf \n", chipl[num]); //Indicate the fit region for the power laws
	for(int i = 0; i<f; i++)
	{
		fprintf(fout,"%0.3lf    %0.8lf\n", tsim[i], Esim[i]);
	}

	//Write signal and bkg separately
	/*for(int i = 0; i<Fitphotonssig; i++)
	{
		fprintf(fout,"%0.3lf    %0.8lf\n", tsimsig[i], Esimsig[i]);
	}
	for(int i = 0; i<Fitphotonsbkg; i++)
	{
		fprintf(fout,"%0.3lf    %0.8lf\n", tsimbkg[i], Esimbkg[i]);
	}*/

	//Write a root file with the important plots
	TString rootoutfile(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau0/Simulationplbkg45%i.root", num));
	TFile* outputfile = new TFile(rootoutfile, "recreate");
	signaltimehisto->SetTitle(" Signal Time generated events");
	signaltimehisto->Write();
	signalenergyhisto->SetTitle("Signal energy generated events");
	signalenergyhisto->Write();
	laghisto->SetTitle(Form("Events with lag %1.2f", tau));
	laghisto->Write();
	smearhistosignal->SetTitle("Signal Smeared events");
	smearhistosignal->Write();
	plhisto->SetTitle("Signal PL template: 0.4-1 TeV");
	//plhisto->SetName("SignalPLtemplate");
	plhisto->SetName("PLtemplate");
	plhisto->Write();
	/*plhistobkg->SetTitle("Bkg PL template: 0.4-1 TeV");
	plhistobkg->SetName("BkgPLtemplate");
	plhistobkg->Write();*/
	timeeventssig->SetTitle("Signal Flare time events: 0.4-1 TeV");
	timeeventssig->Write();
	energyeventssig->SetTitle("Signal Flare energy events: 0.4-1 TeV");
	energyeventssig->Write();
	outputfile->Close();


	delete signaltimehisto;
	delete bkgtimehisto;
	delete signalenergyhisto;
	delete bkgenergyhisto;
	delete laghisto;
	delete plhisto;
	//delete plhistobkg;
	delete LCtemplate;
	delete energyeventssig;
	delete timeeventssig;
	delete energyeventsbkg;
	delete timeeventsbkg;
	delete bkgtimehisto;
	delete smearhistosignal;
	delete smearhistobkg;

}


void sepsimPG()
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
	Double_t sumtmean, sumtime, sumE, sumpls, sumplb, sumratio, sumLC, sumPL, sumcut;
	sumtmean = 0;
	sumtime = 0;
	sumE = 0;
	sumpls = 0;
	sumplb = 0;
	sumratio = 0;
	sumLC = 0;
	sumPL = 0;
	sumcut = 0;

	for(Int_t i=0; i<simnum; i++)
	{
		sumtmean+=tmean[i];
		sumtime+=TimeMean[i];
		sumE+=Emean[i];
		sumpls+=inds[i];
		sumplb+=indb[i];
		sumratio+=SBratio[i];
		sumLC+=Nlc[i];
		sumPL+=Nsp[i];
		sumcut+=Ecut[i];
	}

	Double_t meantmean, meantime, meanE, meanpls, meanplb, meanratio, meanLC, meanPL, meancut;
	meantmean=sumtmean/simnum;
	meantime=sumtime/simnum;
	meanE=sumE/simnum;
	meanpls=sumpls/simnum;
	meanplb=sumplb/simnum;
	meanratio=sumratio/simnum;
	meanLC=sumLC/simnum;
	meanPL=sumPL/simnum;
	meancut=sumcut/simnum;

	TH1D* tmeanhisto = new TH1D("tmeanhisto", "tmeanhisto", 50, meantmean-300, meantmean+300);
	TH1D* indshisto = new TH1D("indshisto", "indshisto", 50, meanpls-4., meanpls+5.); //0.6
	//TH1D* indbhisto = new TH1D("indbhisto", "indbhisto", 50, meanplb-4., meanplb+5.); //0.6
	TH1D* chiplhisto = new TH1D("chiplhisto", "chiplhisto", 50, 0, 5);
	TH1D* Nlchisto = new TH1D("Nlchisto", "Nlchisto", 50, meanLC-100, meanLC+100);
	TH1D* Nsphisto = new TH1D("Nsphisto", "Nsphisto", 50, meanPL-100, meanPL+100);
	TH1D* TimeMeanhisto = new TH1D("TimeMeanhisto", "TimeMeanhisto", 50, meantime-750, meantime+750);
	TH1D* MeanEhisto = new TH1D("MeanEhisto", "MeanEhisto", 50, meanE-0.2, meanE+0.2);
	TH1D* ratiohisto = new TH1D("ratiohisto", "ratiohisto", 50, meanratio-5., meanratio+6.); //1
	//TH1D* cuthisto = new TH1D("cuthisto", "cuthisto", 50, meancut-5., meancut+7.); //1


	for(Int_t i = 0; i<=iter; i++)
	{
		tmeanhisto->Fill(tmean[i]);
		indshisto->Fill(inds[i]);
		//indbhisto->Fill(indb[i]);
		chiplhisto->Fill(chipl[i]);
		Nlchisto->Fill(Nlc[i]);
		Nsphisto->Fill(Nsp[i]);
		TimeMeanhisto->Fill(TimeMean[i]);
		MeanEhisto->Fill(Emean[i]);
		ratiohisto->Fill(SBratio[i]);
		//cuthisto->Fill(Ecut[i]);
	}

	gStyle->SetOptStat(111111);

	TCanvas* big = new TCanvas("tauplbkg45", "tauplbkg45", 1000, 1000);
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
	ratiohisto->Fit("gaus", "Q");
	ratiohisto->Draw("EP");
	big->cd(3);
	gPad->SetTitle("Signal PL index");
	indshisto->SetName("PL index (Signal)");
	indshisto->GetXaxis()->SetTitle("index");
	indshisto->Fit("gaus", "Q");
	indshisto->Draw("EP");
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
	big->SaveAs("./tau0/Canvastau0pgplbkg45.png");



}

