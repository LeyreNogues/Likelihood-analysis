/*
 * Randomize_data.C
 *
 *  Created on: Jul 16, 2018
 *      Author: lnogues
 */

#include<iostream>
#include<fstream>
#include<istream>
#include<TGraph.h>
#include<TGraph2D.h>
#include<TF1.h>
#include<TFile.h>
#include<TH1D.h>
#include<TStyle.h>
#include<TMath.h>
#include<TRandom.h>
#include<TCanvas.h>
#include<TLegend.h>
#include <TVector.h>
#include <TVectorT.h>
#include<TMath.h>
#include "TTree.h"
#include <TROOT.h>
#include <TRandom3.h>
#include <algorithm>    // std::random_shuffle
#include <vector>       // std::vector
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand

using namespace std;


void Randomize_data(Int_t index){

	//Read data
	int Npoints = 0;
	TH1D* data_time;
	TH1D* data_energy;
	ifstream in;ls

	in.open("../../../ReadEvents/Thesis_EventList/Event_Kazuma_120.txt");
//	in.open("./Event_Kazuma_120.txt");//Hola
	Double_t t[20000],E[20000];
	Double_t v1,v2;

	Double_t max_time = 0.;
	while(1)
	{
		in >> v1 >> v2;
		if (!in.good()) break;
		t[Npoints] = v1;
		E[Npoints] = v2;
		Npoints++;
	}
	cout << "Npoints(Number of events): " << Npoints << endl;

	//Randomized energy
	for (unsigned i = 0; i < 20; i++)
	{
		cout << E[i] << " ";
	}
	cout << endl;
	std::srand (unsigned(std::time(0)));

	random_shuffle(&E[0],&E[Npoints-1]);

	for (unsigned i = 0; i < 20; i++)
	{
		cout << E[i] << " ";
	}
	cout << endl;
	cout << "Events randomized!" << endl;

	//Write events
	FILE *fout;
	fout = fopen(Form("./Random_Kazuma_120_%i.txt", index),"wb");
	for(int i = 0; i<Npoints; i++)
	{
		fprintf(fout,"%0.8lf    %0.8lf   \n", t[i], E[i]);

	}
	cout << "Written randomized events! " << endl;


}



