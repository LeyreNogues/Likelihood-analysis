/*
 * DrawDistributions.C
 *
 *  Created on: Oct 14, 2016
 *      Author: lnogues
 */



#include<fstream>
#include<iostream>
#include<TH1.h>
#include<TCanvas.h>

void DrawDistributions(char *fname){

	const Int_t numberOfEvents = 14000;
	Double_t E[numberOfEvents],t[numberOfEvents];
	Double_t tmin, tmax, Emin, Emax;
	tmin=0;
	tmax=13000;
	Emin=100.;
	Emax=20000.;


	ifstream in;
	in.open(fname);
	double v1,v2;
	int Npoints = 0;
	while(1)
	{
	  in >> v1 >> v2;
	  if (!in.good()) break;
	  t[Npoints] = v1;
	  E[Npoints] = v2;
	  Npoints++;
	}

	TH1D* energyint = new TH1D("energyint", "energyint", 50, Emin, Emax);
	TH1D* timedel = new TH1D("timedel", "timedel", 200, tmin, tmax);
	for(int i=0; i<Npoints; i++){
		energyint->Fill(E[i]);
		timedel->Fill(t[i]);
	}

	TCanvas* dist = new TCanvas();
	dist->Divide(2,1);
	dist->cd(1);
	dist->cd(1)->SetLogx();
	dist->cd(1)->SetLogy();
	energyint->Draw("");
	dist->cd(2);
	timedel->SetLineColor(2);
	timedel->Draw("");

}
