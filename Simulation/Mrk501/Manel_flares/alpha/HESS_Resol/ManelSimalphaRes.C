/*
 * ManelSim.C
 *
 *  Created on: Oct 3, 2016
*   Author: lnogues
*
*   This codes generates data for the flare of Mrk501 in 2005 following
*   the simulation developed by Manel Martinez is the fortran original
*   code.
*
 */




void ManelSimalphaRes(double alpha){

	//Variables
	const Int_t numberOfEvents = 1491;
	Double_t E[numberOfEvents],t[numberOfEvents];
	Double_t tmin, tmax, Emin, Emax;
	tmin=0;
	tmax=2750.;
	Emin=150.;
	Emax=12000;

	//Spectrum variables
	Double_t signalIndex=-2.4;
	Double_t bkgIndex=-2.7;
	Double_t specIntegralSignal=(1./(signalIndex+1))*(TMath::Power(Emax,signalIndex+1)-TMath::Power(Emin,signalIndex+1));
	Double_t specIntegralbkg=(1./(bkgIndex+1))*(TMath::Power(Emax,bkgIndex+1)-TMath::Power(Emin,bkgIndex+1));
	std::cout << "Integral Spectrum signal = " << specIntegralSignal << std::endl;
	std::cout << "Integral Spectrum bkg = " << specIntegralbkg << std::endl;

   //LC inputs
	Double_t gausMax, gausWidth, SBratio;
	gausMax=2004.9;
	gausWidth=221.51;
	SBratio=0.369;

	//QG inputs
	Double_t z=0.034;
	Double_t H0=70.;
	Double_t PlanckScale=1.2*TMath::Power(10,19);
	Double_t QGpar=alpha;
	Double_t PcTom=3.085*TMath::Power(10,19);
	Double_t QGfixedPart=z*PcTom/H0;
	std::cout << "QG fixed part = " << QGfixedPart << std::endl;

	//Separated histograms for signal and background
	TH1D* signaltime = new TH1D("signaltime", "signaltime", 20, tmin, tmax);
	TH1D* bkgtime = new TH1D("bkgtime", "bkgtime", 20, tmin, tmax);


	//Simulation
	Double_t bkgProb=(SBratio*(tmax-tmin)*specIntegralbkg)/(SBratio*(tmax-tmin)*specIntegralbkg+(1-SBratio)*50.*specIntegralSignal);
	std::cout << "BkgProb = " << bkgProb << std::endl;
	TRandom* ran1 = new TRandom();
	ran1->SetSeed(0);
	TRandom* ran11 = new TRandom();
	ran11->SetSeed(0);
	TRandom* ran12 = new TRandom();
	ran12->SetSeed();
	for(int i=0; i<numberOfEvents; i++){
		Double_t randomNum = ran11->Uniform(0,1);
		if(randomNum<bkgProb){
//			std::cout << "It's bkg" << std::endl;

			E[i] = TMath::Power((ran1->Uniform()*specIntegralbkg*(bkgIndex+1)+TMath::Power(Emin,bkgIndex+1)),(1/(bkgIndex+1)));
			t[i] = tmin+(tmax-tmin)*ran12->Uniform();
			bkgtime->Fill(t[i]);
		}
		else{
//			std::cout << "It's signal" << std::endl;
			E[i] = TMath::Power((ran1->Uniform()*specIntegralSignal*(signalIndex+1)+TMath::Power(Emin,signalIndex+1)),(1/(signalIndex+1)));
			Double_t a, b;
			gRandom->Rannor(a,b);
			t[i] = gausMax+gausWidth*a;
			signaltime->Fill(t[i]);
		}

	}


	/*TCanvas* disttime = new TCanvas();
	disttime->Divide(2,1);
	disttime->cd(1);
	bkgtime->Draw("");
	disttime->cd(2);
	signaltime->Draw("");*/

	//Injection of delay
	Double_t td[numberOfEvents];
	Double_t QGDelay= QGfixedPart*QGpar/PlanckScale;
	for(int i = 0; i<numberOfEvents; i++)
		{
			td[i]=t[i]+QGDelay*E[i];
		}

	//Smearing of energy
	Double_t Es[numberOfEvents];
	Double_t a,b;
	Double_t resolution = 0.10;
	gRandom->SetSeed();
	for(int i = 0; i<numberOfEvents; i++)
		{
			gRandom->Rannor(a,b);
			Es[i]=E[i]*(1.+b*resolution);
		}

	//Histograms of intrinsic time and energy

	TH1D* energyint = new TH1D("energyint", "Energy", 12, Emin, Emax);
	TH1D* energyrec = new TH1D("energyrec", "energyrec", 12, Emin, Emax);
	TH1D* timeint = new TH1D("timeint", "Time", 20, tmin, tmax);
	TH1D* timedel = new TH1D("timedel", "timedel", 20, tmin, tmax);
	for(int i=0; i<numberOfEvents; i++){
		energyint->Fill(E[i]);
		energyrec->Fill(Es[i]);
		timeint->Fill(t[i]);
		timedel->Fill(td[i]);
	}

	/*gStyle->SetOptStat(0);
	TCanvas* dist = new TCanvas();
	dist->Divide(2,1);
	dist->cd(1);
	energyint->Draw("EB");
	energyrec->SetLineColor(2);
	energyrec->Draw("EBsame");
	dist->cd(2);
	timeint->Draw("");
	timedel->SetLineColor(2);
	timedel->Draw("same");


	dist->cd(1);
	TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(energyint,"Intrinsic","le");
	leg1->AddEntry(energyrec,"Smeared","le");
	leg1->Draw();

	dist->cd(2);
	TLegend* leg2 = new TLegend(0.1,0.7,0.4,0.9);
	leg2->AddEntry(timeint,"Intrinsic","l");
	leg2->AddEntry(timedel,"Delayed","l");
	leg2->Draw();*/


	//Create a txt file with the events
	FILE *fout;
	fout = fopen(Form("Manel_flare_res_alpha%i.txt", (int)alpha),"wb");
	for(int i = 0; i<numberOfEvents; i++)
	{
		fprintf(fout,"%0.8lf    %0.8lf   \n", td[i], Es[i]);
	}




}
