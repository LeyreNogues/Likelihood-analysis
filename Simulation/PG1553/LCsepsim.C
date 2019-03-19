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
#include <TRandom.h>
#include <iostream>

using namespace std;

const Int_t simnum = 100;
Double_t tinit = 0.;
Double_t tfin = 8000.;

Double_t tmean1[simnum], tmean2[simnum], ind[simnum], chipl[simnum], chilc[simnum], Nlc[simnum], Nsp[simnum], Emean[simnum], TimeMean[simnum], SBratio[simnum];

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

Double_t doblegaus(Double_t* x, Double_t* par) //with background
{
	return par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2))+par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5],2));

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

	Double_t bkg, bkg2, tmax, sigma,tmax2, sigma2;
	bkg = 18.32;
	tmax = 2540.;
	sigma = 798;
	bkg2 = 15.62;
	tmax2 = 5898.;
	sigma2 = 1070.;

	TF1 *lightcurve = new TF1("lightcurve",doblegaus, tinit, tfin, 6);
	lightcurve->SetParameters(bkg, tmax, sigma, bkg2, tmax2, sigma2);
	lightcurve->SetParName(0, "h_{1}");
	lightcurve->SetParName(1, "#nu_{1}");
	lightcurve->SetParName(2, "#sigma_{1}");
	lightcurve->SetParName(3, "h_{2}");
	lightcurve->SetParName(4, "#nu_{2}");
	lightcurve->SetParName(5, "#sigma_{2}");

	TH1D* signaltimehisto = new TH1D("time events", "time events", 50, tinit, tfin);
	TH1D* bkgtimehisto = new TH1D("bkgtimehisto", "bkgtimehisto", 50, tinit, tfin);


	TRandom* bkgr = new TRandom();
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
		tbkg[i]=bkgr->Uniform(tinit, tfin);
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
	indexb = 2.5;

	TF1* signalspectrum = new TF1("signalspectrum", powerLaw, Emin, Emax, 1);
	signalspectrum->SetParameter(0, indexs);
	TF1* bkgspectrum = new TF1("bkgspectrum", powerLaw, Emin, Emax, 1);
	bkgspectrum->SetParameter(0, indexb);
	TF1* fitspectrum = new TF1("fitspectrum", powerLawfit, 0.4, 1., 2);
	fitspectrum->SetName("fitspectrum");
	fitspectrum->SetParName(0, "norm");
	fitspectrum->SetParName(1, "index");

	double XX_step = 0.1;
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
	canvas2->SetLogy();
	signalenergyhisto->SetName("Energy of events");
	signalenergyhisto->SetLineColor(2);
	signalenergyhisto->Fit(fitspectrum);
	signalenergyhisto->Draw("EP");*/
	//bkgenergyhisto->Draw("sameEP");
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


	//-------------------------------------LC and PL template----------------------------------------

	//SIGNAL ONLY
	double XX3_step = 0.08;
	int n3 = ceil((log10(1.)-log10(0.4))/XX3_step);
	double *XX3 = new double[n3+1];
	for(int i=0;i<n3+1;i++)
	{
	XX3[i]= pow(10,(log10(0.4)+i*XX3_step));
	}
	TH1D* plhisto = new TH1D("plhistosig", "plhistosig", n3, XX3);
	TH1D* lchisto = new TH1D("lchisto", "lchisto", 20, tinit, tfin);

	for(Int_t a=0; a<numev; a++) //Only signal
	{
		if(Esmearedsig[a]>=0.3 && Esmearedsig[a]<=0.4)
		{
			lchisto->Fill(tlag[a]);
		}
		if(Esmearedsig[a]>=0.4 && Esmearedsig[a]<=1.)
		{
			plhisto->Fill(Esmearedsig[a], 1./Esmearedsig[a]);
		}
	}

	/*TCanvas* plcanvas = new TCanvas();
	plcanvas->SetLogx();
	plcanvas->SetLogy();*/
	//plhisto->GetXaxis()->SetTitle("Energy(TeV)");
	plhisto->Fit(fitspectrum, "0QR");
	//plhisto->Draw("EP");
	TF1* fitpl = plhisto->GetFunction("fitspectrum");
	ind[num] = fitpl->GetParameter(1);
	chipl[num] = (fitpl->GetChisquare())/(fitpl->GetNDF());

	//TCanvas* lccanvas = new TCanvas();
	//lchisto->GetXaxis()->SetTitle("Time(s)");
	lchisto->Fit("lightcurve", "0Q");
	//lchisto->Draw("EP");
	TF1* fitlc = lchisto->GetFunction("lightcurve");
	tmean1[num] = fitlc->GetParameter(1);
	tmean2[num] = fitlc->GetParameter(4);
	chilc[num] = (fitlc->GetChisquare())/(fitlc->GetNDF());


	//---------------------------------Division in LC and Fit--------------------------------

	Int_t LCphotonssig = 0;
	Int_t Fitphotonssig = 0;
	Int_t LCphotonsbkg = 0;
	Int_t Fitphotonsbkg = 0;
	TH1D* energyeventssig = new TH1D("energyeventssig", "energyeventssig", n, XX);
	TH1D* timeeventssig = new TH1D("timeeventssig", "timeeventssig", 20, tinit, tfin);
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


	Nlc[num] = LCphotonssig+LCphotonsbkg;
	Nsp[num] = Fitphotonssig+Fitphotonsbkg;

	//cout << "Nlc: " << Nlc[num] << endl;
	//cout << "Nsp: " << Nsp[num] << endl;

	//Write a txt file with the evetns
	FILE *fout;
	fout = fopen(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau0/flarepglc%i.txt", num),"wb");
	fprintf(fout,"%0.3lf \n", LCphotonssig+LCphotonsbkg);
	fprintf(fout,"%0.3lf \n", Fitphotonssig+Fitphotonsbkg);
	fprintf(fout,"%0.3lf \n", chipl[num]);
	for(int i = 0; i<f; i++)
	{
		fprintf(fout,"%0.3lf    %0.8lf\n", tsim[i], Esim[i]);
	}

	//Write a root file with the important plots
	TString rootoutfile(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau0/Simulationpglc%i.root", num));
	TFile* outputfile = new TFile(rootoutfile, "recreate");
	signaltimehisto->SetTitle(" Signal Time generated events");
	signaltimehisto->Write();
	signalenergyhisto->SetTitle("Signal energy generated events");
	signalenergyhisto->Write();
	laghisto->SetTitle(Form("Events with lag %1.2f", tau));
	laghisto->Write();
	smearhistosignal->SetTitle("Signal Smeared events");
	smearhistosignal->Write();
	plhisto->SetTitle("PL template: 0.4-1 TeV");
	plhisto->SetName("PLtemplate");
	plhisto->Write();
	lchisto->SetTitle("LC template: 0.3-0.4 TeV");
	lchisto->SetName("LCtemplate");
	lchisto->Write();
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
	delete lchisto;
	delete energyeventssig;
	delete timeeventssig;
	delete energyeventsbkg;
	delete timeeventsbkg;
	delete bkgtimehisto;
	delete smearhistosignal;
	delete smearhistobkg;

}


void LCsepsim()
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
	Double_t sumtmean1, sumtmean2, sumchipl, sumchilc, sumpl, sumratio, sumLC, sumPL;
	sumtmean1 = 0;
	sumtmean2 = 0;
	sumpl = 0;
	sumratio = 0;
	sumLC = 0;
	sumPL = 0;
	sumchilc = 0,
	sumchipl = 0;

	for(Int_t i=0; i<simnum; i++)
	{
		sumtmean1+=tmean1[i];
		sumtmean2+=tmean2[i];
		sumpl+=ind[i];
		sumratio+=SBratio[i];
		sumLC+=Nlc[i];
		sumPL+=Nsp[i];
		sumchipl+=chipl[i];
		sumchilc+=chilc[i];

	}

	Double_t meantmean1, meantmean2, meanpl, meanratio, meanLC, meanPL, meanchipl, meanchilc;
	meantmean1=sumtmean1/simnum;
	meantmean2=sumtmean2/simnum;
	meanpl=sumpl/simnum;
	meanratio=sumratio/simnum;
	meanLC=sumLC/simnum;
	meanPL=sumPL/simnum;
	meanchipl=sumchipl/simnum;
	meanchilc=sumchilc/simnum;


	TH1D* tmeanhisto1 = new TH1D("tmeanhisto1", "tmeanhisto1", 50, meantmean1-700, meantmean1+700);
	TH1D* tmeanhisto2 = new TH1D("tmeanhisto2", "tmeanhisto2", 50, meantmean2-700, meantmean2+700);
	TH1D* indhisto = new TH1D("indhisto", "indhisto", 50, meanpl-4., meanpl+5.); //0.6
	TH1D* chiplhisto = new TH1D("chiplhisto", "chiplhisto", 50, meanchipl-3, meanchipl+3);
	TH1D* chilchisto = new TH1D("chilchisto", "chilchisto", 50, meanchilc-2, meanchilc+2);
	TH1D* Nlchisto = new TH1D("Nlchisto", "Nlchisto", 50, meanLC-100, meanLC+100);
	TH1D* Nsphisto = new TH1D("Nsphisto", "Nsphisto", 50, meanPL-100, meanPL+100);
	TH1D* ratiohisto = new TH1D("ratiohisto", "ratiohisto", 50, meanratio-3., meanratio+4.); //1


	for(Int_t i = 0; i<=iter; i++)
	{
		tmeanhisto1->Fill(tmean1[i]);
		tmeanhisto2->Fill(tmean2[i]);
		indhisto->Fill(ind[i]);
		chiplhisto->Fill(chipl[i]);
		chilchisto->Fill(chilc[i]);
		Nlchisto->Fill(Nlc[i]);
		Nsphisto->Fill(Nsp[i]);
		ratiohisto->Fill(SBratio[i]);
	}

	gStyle->SetOptStat(111111);

	TCanvas* big = new TCanvas("tau0simlc", "tau0simlc", 1000, 1000);
	big->Divide(2,4);
	big->cd(1);
	gPad->SetTitle("LC Tmax1");
	tmeanhisto1->SetName("LC Tmax1");
	tmeanhisto1->GetXaxis()->SetTitle("s");
	tmeanhisto1->Draw("EP");
	big->cd(2);
	gPad->SetTitle("LC Tmax2");
	tmeanhisto2->SetName("LC Tmax2");
	tmeanhisto2->GetXaxis()->SetTitle("s");
	tmeanhisto2->Draw("EP");
	big->cd(3);
	gPad->SetTitle("Chi2/NDF(LC)");
	chilchisto->SetName("#Chi^{2}/NDF(LC)");
	chilchisto->GetXaxis()->SetTitle("#Chi^{2}/NDF");
	chilchisto->Draw("EP");
	big->cd(4);
	gPad->SetTitle("S/B ratio");
	ratiohisto->SetName("S/B ratio");
	ratiohisto->GetXaxis()->SetTitle("ratio(0.25-1TeV)");
	ratiohisto->Fit("gaus", "Q");
	ratiohisto->Draw("EP");
	big->cd(5);
	gPad->SetTitle("PL index");
	indhisto->SetName("PL index");
	indhisto->GetXaxis()->SetTitle("index");
	indhisto->Fit("gaus", "Q");
	indhisto->Draw("EP");
	big->cd(6);
	gPad->SetTitle("Chi2/NDF(PL)");
	chiplhisto->SetName("#Chi^{2}/NDF(PL)");
	chiplhisto->GetXaxis()->SetTitle("#Chi^{2}/NDF");
	chiplhisto->Draw("EP");
	big->cd(7);
	Nlchisto->SetTitle("Events LC");
	Nlchisto->SetName("Events LC");
	Nlchisto->GetXaxis()->SetTitle("events");
	Nlchisto->Draw("EP");
	big->cd(8);
	Nsphisto->SetTitle("Events Fit");
	Nsphisto->SetName("Events Fit");
	Nsphisto->GetXaxis()->SetTitle("events");
	Nsphisto->Draw("EP");
	big->SaveAs("./tau0/Canvastau0pglc.png");



}

