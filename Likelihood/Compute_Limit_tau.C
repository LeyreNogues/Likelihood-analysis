/*
 * Compute_Limit_tau.C
 *
 *  Created on: Aug 5, 2018
 *      Author: lnogues
 */

#include<TFile.h>
#include<TCanvas.h>
#include<TH1.h>
#include <fstream>
#include <iostream>
#include <TMath.h>
#include <TF1.h>

using namespace std;

// type: 1 linear, 2 Quadratic
// lag: measured lag
// stat: statistical error
// sys: systematic error
// nsig: 2 for 95%CL, 3 for 99%CL, etc.


Double_t Compute_Limit_tau(int type, double lag, double stat, double sys, double nsig = 2.) {

	Double_t H_0_s_1 = 2.29e-18;
	Double_t z =0.031;

	cout << "Compute LIV tau limit for case: " << type << endl;
	cout << "The mean value is: " << lag << endl;
	cout << "The stat is: " << stat << endl;
	cout << "The syst is: " << sys << endl;
	cout << "The nsig is: " << nsig << endl;

	Double_t res; //lag_nsig = lag + nsig * sqrt(stat * stat + sys * sys);

	Double_t lag_nsig = stat ;

	cout << "Final tau: " << lag_nsig << endl;

	if(type == 1) {
		res = z/H_0_s_1/lag_nsig;
	}
	else if(type == 2) {
		res = TMath::Sqrt((3/2.)*z/H_0_s_1/lag_nsig);
	}

	return res;
}

