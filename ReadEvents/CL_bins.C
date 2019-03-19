/*
 * CL_bins.C
 *
 *  Created on: Jul 21, 2018
 *      Author: lnogues
 */

#include<iostream>
#include<fstream>
#include<istream>
#include <math.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TStyle.h>
#include<iostream>
#include<fstream>
#include<istream>
#include<TVectorD.h>
#include<TMath.h>

using namespace std;

void CL_bins(Double_t dist, Double_t en1,Double_t en2){

	const Int_t number = 11;
	Double_t C_alpha[number];
	Double_t CL_level[number];
	Double_t Kappa_alpha[number];

	CL_level[0] = 0.60;
	CL_level[1] = 0.50;
	CL_level[2] = 0.40;
	CL_level[3] = 0.20;
	CL_level[4] = 0.15;
	CL_level[5] = 0.10;
    CL_level[6] = 0.05;
	CL_level[7] = 0.025;
	CL_level[8] = 0.01;
	CL_level[9] = 0.005;
	CL_level[10] = 0.001;

	C_alpha[0] = 0.78;
	C_alpha[1] = 0.83;
	C_alpha[2] = 0.9;
	C_alpha[3] = 1.07;
	C_alpha[4] = 1.14;
	C_alpha[5] = 1.22;
	C_alpha[6] = 1.36;
	C_alpha[7] = 1.48;
	C_alpha[8] = 1.63;
	C_alpha[9] = 1.73;
	C_alpha[10] = 1.95;

	cout << "Dist: " << dist << endl;
	cout << "en1: " << en1 << " en2: " << en2 << endl;

	for(Int_t i=0; i<number; i++){
		Double_t product = en1*en2;
		Double_t sum =  en1+en1;
//		cout << "product: " << product << endl;
//		cout << "sum: " << sum << endl;
		Kappa_alpha[i] = C_alpha[i]*TMath::Sqrt(sum/product);
		cout << "CL: " << CL_level[i]  << " c(alpha): " << C_alpha[i] << " Kappa_alha: " << Kappa_alpha[i] << endl;
	}
}

	Double_t alpha_crit(Double_t dist, Double_t en1,Double_t en2){

		Double_t A = dist*dist*(-2)*en1*en2/(en1+en2);
		Double_t exp_A = TMath::Exp(A);
		Double_t alpha_crit = 2*exp_A;
		return alpha_crit;

	}




