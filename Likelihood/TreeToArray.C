/*
 * TreeToArray.C
 *
 *  Created on: Sep 6, 2018
 *      Author: lnogues
 */
#include <Rtypes.h>
#include <stddef.h>
#include <TArrayD.h>
#include <TAttLine.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TString.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include <MMath.h>

using namespace std;

void TreeToArray(){

	Int_t alpha_index = 0;
	Double_t alpha_array[1000]={0};

	for(Int_t i=0; i<5; i++){ //1000
		Int_t index=i;
		TString InputfileName = Form("Likelihood_Kazuma_120_Inter_Quadratic_NoBkg%i.root",index);
		cout << "File name: " << InputfileName << endl;

		TFile* Eventfile = TFile::Open(InputfileName, "READ");
		if (!Eventfile) {
			cout << "ERROR: Eventfile " << Eventfile << " not found..." << endl;
		}
		else{
			cout << "I found the file" << endl;

			//Read the Tree
			TTree *t1 = (TTree*)Eventfile->Get("t1");
			if (!t1) {
				cout << "ERROR: tree file " << t1 << " not found..." << endl;
			}

			cout << "I read the tree" << endl;
			Double_t alpha_value;

			t1->SetBranchAddress("alpha",&alpha_value);
			t1->GetEntry(0);
			cout << "alpha_value: " << alpha_value << endl;

			alpha_array[alpha_index] = alpha_value;

			cout << "component: " << alpha_array[alpha_index] << endl;
			alpha_index++;

			Eventfile->Close();

		}
	}

	cout << "Number of found files: " << alpha_index << endl;

	//Compute MAD and Sn values
	Double_t MAD_value = MMath::MAD(alpha_index,alpha_array);
	Double_t Sn_value = MMath::Sn(alpha_index,alpha_array);
	cout << "MAD_value: " << MAD_value << endl;
	cout << "Sn_value: " << Sn_value << endl;

	//Compute standard deviation
	Double_t alpha_mean = 0;
	Double_t alpha_RMS = 0;
	for(Int_t i=0; i<alpha_index; i++){
		alpha_mean+=alpha_array[i];
		alpha_RMS+=alpha_array[i]*alpha_array[i];
	}

	alpha_mean = alpha_mean/alpha_index;
	alpha_RMS = 2*TMath::Sqrt((alpha_RMS-(alpha_mean*alpha_mean*alpha_index))/(alpha_index-1.));
	cout << "alpha mean: " << alpha_mean << endl;
	cout << "Sta. Deviation: " << alpha_RMS << endl;



}


