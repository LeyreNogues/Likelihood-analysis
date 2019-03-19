/*
 * Select_Samples.C
 *
 *  Created on: Mar 22, 2018
 *      Author: lnogues
 *
 *
 *      This codes allows to select small subsamples from a tree of events
 *      with lots of statistics. The  number of events per sample can be modified.
 *
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
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>


using namespace std;

void Select_Samples(Int_t numberOfEvents, Int_t index){
	//Tree version
	cout << "I will read now" << endl;
	TString pathToReadEvents = "./";
	const TString fileName6 = pathToReadEvents + "Event_tree_Mrk421_20140_0_KDE.root";
	TFile* Eventfile = TFile::Open(fileName6, "READ");
	if (!Eventfile) {
		cout << "ERROR: file " << Eventfile << " not found..." << endl;
	}

	cout << "I found the file" << endl;

	//Read the Tree
	TTree *tevents = (TTree*)Eventfile->Get("Event_tree_Mrk421_20140_0_KDE.root");
	if (!tevents) {
		cout << "ERROR: file " << tevents << " not found..." << endl;
	}
	Double_t Es, td;
	tevents->SetBranchAddress("Es",&Es);
	tevents->SetBranchAddress("td",&td);


	//--Events in tree
	vector<double> E, t, E_sel, t_sel;
	E.clear();
	t.clear();
	E_sel.clear();
	t_sel.clear();

	Long64_t nevents = tevents->GetEntries();
	cout << "Number of events in the tree: " << nevents << endl;
	Int_t Npoints = 0;
	for(Long64_t i=0;i<nevents;i++){
		tevents->GetEntry(i);
		E.push_back(Es);
		t.push_back(td);
		Npoints++;
	}

	cout << "Number of events " << numberOfEvents << endl;
	cout << "Length of the array E: " << E.size() << endl;
	cout << "Length of the array t: " << t.size() << endl;


	//---Loop: The original array is sampled and the sampled elements are removed out of the
	//array in order not to sample them again. Even though the periodicity of the
	//TRandom 3 is much larger than the array length...(numberOfEvents: 12994)

	cout << "Going to loop..." << endl;
	TRandom3 *ran = new TRandom3();
	ran->SetSeed(0);
	for(Int_t i=0; i<numberOfEvents; i++){
		Int_t number = E.size();
//		cout << "Legth: " << number << endl;
		Int_t component = ran->Uniform(0,number);
//		cout << "component: " << component << endl;
//		cout << "E[component]: " << E[component] << endl;

		E_sel.push_back(E.at(component));
		t_sel.push_back(t.at(component));
		E.erase(E.begin()+component);
		t.erase(t.begin()+component);

//		cout << "Length of E: " << E.size() << endl;
//		cout << "Length of E_sel: " << E_sel.size() << endl;

	  }

	//Create a txt file with the results of the Selection.
	Int_t VHE_events = 0;
	FILE *fout;
	fout = fopen(Form("./Sample_10res%i_KDE.txt", index),"wb");

	for(int i = 0; i<numberOfEvents; i++)
	{
		fprintf(fout,"%0.8lf    %0.8lf   \n", t_sel[i], E_sel[i]);
		if(E_sel[i]>4000) VHE_events++;

	}
	cout << "VHE events: " << VHE_events << endl;
	cout << "Written sample events! " << endl;




}

