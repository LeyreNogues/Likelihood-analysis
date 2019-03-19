/*
 * Select_Root.C
 *
 *  Created on: Apr 12, 2018
 *      Author: lnogues
 */


#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <iostream>

using namespace std;


void Select_Root(Double_t minalpha, Double_t maxalpha){

Int_t numberOfFiles = 4;
TFile* RootFile;

for(Int_t i=0; i<numberOfFiles; i++){

	RootFile = TFile::Open(Form("LikelihoodMrk421_SAMPLE_10res%i.root", i), "READ");
		if (!RootFile) {
			cout << "ERROR: file " << RootFile << " not found..."	<< endl;
		}
	TTree *t = (TTree*)RootFile->Get("t1");
	Double_t alpha;
	t->SetBranchAddress("alpha",&alpha);

	t->GetEvent(0);
	cout << "alpha: " << alpha << endl;

	if (alpha > minalpha && alpha < maxalpha){ cout << RootFile->GetName() << endl; }

	}
}
