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




void DoubleSimalpha(double alpha){

	//Variables
	const Int_t numberOfEvents_flare1 = 1491;
	const Int_t numberOfEvents_flare2 = 1491;
	const Int_t numberOfEvents = numberOfEvents_flare1+numberOfEvents_flare2;
	Double_t E1[numberOfEvents_flare1],t1[numberOfEvents_flare1], td1[numberOfEvents_flare1];
	Double_t E2[numberOfEvents_flare2],t2[numberOfEvents_flare2], td2[numberOfEvents_flare2];
	Double_t tmin1, tmax1, Emin1, Emax1;
	tmin1=0;
	tmax1=2750.;
	Emin1=150.;
	Emax1=12000;
	Double_t tmin2, tmax2, Emin2, Emax2;
	tmin2=2750;
	tmax2=5500.;
	Emin2=150.;
	Emax2=12000;

	//Spectrum variables
	Double_t signalIndex1=-2.4;
	Double_t bkgIndex1=-2.7;
	Double_t signalIndex2=-2.4;
	Double_t bkgIndex2=-2.7;
	Double_t specIntegralSignal1=(1./(signalIndex1+1))*(TMath::Power(Emax1,signalIndex1+1)-TMath::Power(Emin1,signalIndex1+1));
	Double_t specIntegralbkg1=(1./(bkgIndex1+1))*(TMath::Power(Emax1,bkgIndex1+1)-TMath::Power(Emin1,bkgIndex1+1));
	Double_t specIntegralSignal2=(1./(signalIndex2+1))*(TMath::Power(Emax2,signalIndex2+1)-TMath::Power(Emin2,signalIndex2+1));
	Double_t specIntegralbkg2=(1./(bkgIndex2+1))*(TMath::Power(Emax2,bkgIndex2+1)-TMath::Power(Emin2,bkgIndex2+1));


   //LC inputs
	Double_t gausMax1, gausWidth1, SBratio1;
	Double_t gausMax2, gausWidth2, SBratio2;
	gausMax1=2004.9;
	gausWidth1=221.51;
	SBratio1=0.369;
	gausMax2=4754.9;
	gausWidth2=221.51;
	SBratio2=0.369;

	//QG inputs
	Double_t z=0.034;
	Double_t H0=70.;
	Double_t PlanckScale=1.2*TMath::Power(10,19);
	Double_t QGpar=alpha;
	Double_t PcTom=3.085*TMath::Power(10,19);
	Double_t QGfixedPart=z*PcTom/H0;


	//Common simulation parameters
	TRandom* ran1 = new TRandom();
	ran1->SetSeed(0);
	TRandom* ran11 = new TRandom();
	ran11->SetSeed(0);
	TRandom* ran12 = new TRandom();
	ran12->SetSeed(0);
	Double_t QGDelay= QGfixedPart*QGpar/PlanckScale;



	//Simulation flare 1
	Double_t bkgProb1=(SBratio1*(tmax1-tmin1)*specIntegralbkg1)/(SBratio1*(tmax1-tmin1)*specIntegralbkg1+(1-SBratio1)*50.*specIntegralSignal1);
	std::cout << "BkgProb1 = " << bkgProb1 << std::endl;
	for(int i=0; i<numberOfEvents_flare1; i++){
		Double_t randomNum = ran11->Uniform();
		if(randomNum<bkgProb1){
			E1[i] = TMath::Power((ran1->Uniform(0,1)*specIntegralbkg1*(bkgIndex1+1)+TMath::Power(Emin1,bkgIndex1+1)),(1/(bkgIndex1+1)));
			t1[i] = tmin1+(tmax1-tmin1)*ran12->Uniform();
			td1[i]=t1[i]+QGDelay*E1[i];
		}
		else{
			E1[i] = TMath::Power((ran1->Uniform()*specIntegralSignal1*(signalIndex1+1)+TMath::Power(Emin1,signalIndex1+1)),(1/(signalIndex1+1)));
			Double_t a, b;
			gRandom->Rannor(a,b);
			t1[i] = gausMax1+gausWidth1*a;
			td1[i]=t1[i]+QGDelay*E1[i];
		}
	}

	//Simulation flare 2
	Double_t bkgProb2=(SBratio2*(tmax2-tmin2)*specIntegralbkg2)/(SBratio2*(tmax2-tmin2)*specIntegralbkg2+(1-SBratio2)*50.*specIntegralSignal2);
	std::cout << "BkgProb2 = " << bkgProb2 << std::endl;
	for(int i=0; i<numberOfEvents_flare2; i++){
		Double_t randomNum = ran11->Uniform();
		if(randomNum<bkgProb2){
			E2[i] = TMath::Power((ran1->Uniform()*specIntegralbkg2*(bkgIndex2+1)+TMath::Power(Emin2,bkgIndex2+1)),(1/(bkgIndex2+1)));
			t2[i] = tmin2+(tmax2-tmin2)*ran12->Uniform();
			td2[i]=t2[i]+QGDelay*E2[i];
		}
		else{
			E2[i] = TMath::Power((ran1->Uniform()*specIntegralSignal2*(signalIndex2+1)+TMath::Power(Emin2,signalIndex2+1)),(1/(signalIndex2+1)));
			Double_t a, b;
			gRandom->Rannor(a,b);
			t2[i] = gausMax2+gausWidth2*a;
			td2[i]=t2[i]+QGDelay*E2[i];
		}
	}


	//Put both flares together
	Double_t E[numberOfEvents],t[numberOfEvents],td[numberOfEvents],index[numberOfEvents];;
	for(int i=0; i<numberOfEvents; i++){
		if(i<numberOfEvents_flare1){
			index[i]=1;
			E[i]=E1[i];
			t[i]=t1[i];
			td[i]=td1[i];

		}
		else{
			index[i]=2;
			E[i]=E2[i-numberOfEvents_flare2];
			t[i]=t2[i-numberOfEvents_flare2];
			td[i]=td2[i-numberOfEvents_flare2];
		}

	}

	//Plot the separated flares
	TH1D* energyint1 = new TH1D("energyint1", "Energy Distribution", 12, Emin1, Emax1);
	TH1D* timeint1 = new TH1D("timeint1", "Time Distribution", 20, tmin1, tmax1);
	TH1D* energyint2 = new TH1D("energyint2", "Energy Distribution", 12, Emin2, Emax2);
	TH1D* timeint2 = new TH1D("timeint2", "Time Distribution", 20, tmin2, tmax2);
	for(int i=0; i<numberOfEvents_flare1; i++){
		energyint1->Fill(E1[i]);
		timeint1->Fill(t1[i]);
		energyint2->Fill(E2[i]);
		timeint2->Fill(t2[i]);
	}

	gStyle->SetOptStat(0);

	/*TCanvas* dist1 = new TCanvas();
	dist1->Divide(2,1);
	dist1->cd(1);
	energyint1->GetXaxis()->SetTitle("Energy (GeV)");
	energyint1->SetMarkerStyle(31);
	energyint1->Draw("P");
	dist1->cd(2);
	timeint1->GetXaxis()->SetTitle("Time(s)");
	timeint1->Draw("");

	TCanvas* dist2 = new TCanvas();
	dist2->Divide(2,1);
	dist2->cd(1);
	energyint2->GetXaxis()->SetTitle("Energy (GeV)");
	energyint2->SetMarkerStyle(31);
	energyint2->Draw("P");
	dist2->cd(2);
	timeint2->GetXaxis()->SetTitle("Time(s)");
	timeint2->Draw("");*/



	//Histograms of intrinsic time and energy (Both flares together)

	TH1D* energyint = new TH1D("energyint", "Energy Distribution", 12, Emin1, Emax2);
	TH1D* timeint = new TH1D("timeint", "Time Distribution", 40, tmin1, tmax2);
	TH1D* timedel = new TH1D("timedel", "timedel", 40, tmin1, tmax2);
	for(int i=0; i<numberOfEvents; i++){
		energyint->Fill(E[i]);
		timeint->Fill(t[i]);
		timedel->Fill(td[i]);
	}

	/*gStyle->SetOptStat(0);

	TCanvas* dist = new TCanvas();
	dist->Divide(2,1);
	dist->cd(1);
	dist->SetLogx();
	dist->SetLogy();
	energyint->GetXaxis()->SetTitle("Energy (GeV)");
	energyint->SetMarkerStyle(31);
	energyint->Draw("P");
	dist->cd(2);
	timeint->GetXaxis()->SetTitle("Time(s)");
	timeint->Draw("");
	timedel->SetLineColor(2);
	timedel->Draw("same");


	TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
	leg1->AddEntry(timeint,"Intrinsic t","l");
	leg1->AddEntry(timedel,"Delayed t","l");
	leg1->SetTextSize(0.04);
	leg1->Draw();*/


	//Create a txt file with the events
	FILE *fout;
	fout = fopen(Form("Two_flares_alpha%i.txt", (int)alpha),"wb");
	for(int i = 0; i<numberOfEvents; i++)
	{
		fprintf(fout,"%i     %0.8lf    %0.8lf   \n", index[i], td[i], E[i]);
	}




}
