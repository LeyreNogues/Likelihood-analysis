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




void ManelSimalpha(double alpha){

	//Variables
	const Int_t numberOfEvents = 1491;
	Double_t E[numberOfEvents],t[numberOfEvents], td[numberOfEvents];;
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
	ran12->SetSeed(0);
	Double_t QGDelay= QGfixedPart*QGpar/PlanckScale;
	for(int i=0; i<numberOfEvents; i++){
		Double_t randomNum = ran11->Uniform(0,1);
//		std::cout << "Random number " << randomNum << std::endl;

		if(randomNum<bkgProb){
//			std::cout << "It's bkg" << std::endl;
			E[i] = TMath::Power((ran1->Uniform(0,1)*specIntegralbkg*(bkgIndex+1)+TMath::Power(Emin,bkgIndex+1)),(1/(bkgIndex+1)));
			t[i] = tmin+(tmax-tmin)*ran12->Uniform(0,1);
			td[i]=t[i]+QGDelay*E[i];
			bkgtime->Fill(t[i]);
//			std::cout << "Energy = " << E[i] << " GeV"<< std::endl;
//			std::cout << "Time = " << t[i] << " s"<< std::endl;
		}
		else{
//			std::cout << "It's signal" << std::endl;
			E[i] = TMath::Power((ran1->Uniform()*specIntegralSignal*(signalIndex+1)+TMath::Power(Emin,signalIndex+1)),(1/(signalIndex+1)));
			Double_t a, b;
			gRandom->Rannor(a,b);
			t[i] = gausMax+gausWidth*a;
			td[i]=t[i]+QGDelay*E[i];
			signaltime->Fill(t[i]);
//			std::cout << "Energy = " << E[i] << " GeV"<< std::endl;
//			std::cout << "Time = " << t[i] << " s"<< std::endl;
		}

	}


	/*TCanvas* disttime = new TCanvas();
	disttime->Divide(2,1);
	disttime->cd(1);
	bkgtime->Draw("");
	disttime->cd(2);
	signaltime->Draw("");*/

	//Histograms of intrinsic time and energy

	gStyle->SetOptStat(0);
	TH1D* energyint = new TH1D("energyint", "energyint", 12, Emin, Emax);
	TH1D* timeint = new TH1D("timeint", "timeint", 20, tmin, tmax);
	TH1D* timedel = new TH1D("timedel", "timedel", 20, tmin, tmax);
	for(int i=0; i<numberOfEvents; i++){
		energyint->Fill(E[i]);
		timeint->Fill(t[i]);
		timedel->Fill(td[i]);
	}

	/*TCanvas* dist = new TCanvas();
	dist->Divide(2,1);
	dist->cd(1);
	dist->SetLogx();
	dist->SetLogy();
	energyint->Draw("");
	dist->cd(2);
	timeint->Draw("");
	timedel->SetLineColor(2);
	timedel->Draw("same");*/

	TCanvas* dist = new TCanvas();
	timeint->Draw("");
	timeint->SetTitle("");
	timeint->GetXaxis()->SetTitle("Time(s)");
	timedel->SetLineColor(2);
	timedel->Draw("same");

	//Create a txt file with the events
//	FILE *fout;
//	fout = fopen(Form("Manel_flare_alpha%i.txt", (int)alpha),"wb");
//	for(int i = 0; i<numberOfEvents; i++)
//	{
//		fprintf(fout,"%0.8lf    %0.8lf   \n", td[i], E[i]);
//	}




}
