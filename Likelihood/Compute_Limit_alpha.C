/*
 * Compute_Limit_alpha.C
 *
 *  Created on: Apr 13, 2018
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


Double_t Compute_Limit_alpha(Int_t type, Double_t lag, Double_t stat, Double_t sys, Double_t nsig = 2.){

	Double_t PlanckScale=1.2*TMath::Power(10,19);

	cout << "Compute LIV limit for case: " << type << endl;
	cout << "The mean value is: " << lag << endl;
	cout << "The stat is: " << stat << endl;
	cout << "The syst is: " << sys << endl;
	cout << "The nsig is: " << nsig << endl;

	Double_t lag_nsig = lag + nsig * TMath::Sqrt(stat * stat + sys * sys); //Final alpha
//	Double_t lag_nsig = stat; //Final alpha


	cout << "Final alpha: " << lag_nsig << endl;

	Double_t res = 0.; //resolution
	  if(type==1){ //Linear
		res= PlanckScale/lag_nsig;
	  }
	  if(type==2){ //Quadratic
		res= TMath::Sqrt((3000/2.)*(PlanckScale/lag_nsig));
	  }
	  return res;
}


