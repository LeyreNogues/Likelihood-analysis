/*
 * Add_Results.C
 *
 *  Created on: Feb 9, 2017
 *      Author: lnogues
 */


void Add_Results(){

	const Int_t numOfResults = 5; //Number of results to combine
	TH1D* results_histo[numOfResults];
	TList* list_of_histos = new TList;

	for(int i=0; i<numOfResults; i++){
//		TString pathToReadCode = "./PIC_alpha20_results/";
		TString pathToReadCode = "./results_alpha20res/";
		const TString fileName = pathToReadCode + Form("100Sim_alpha20res_%i_histos.root", i+1);
		TFile* file = TFile::Open(fileName, "READ");
		if (!file){
			std::cout << "ERROR: file " << file << " not found..." << std::endl;
		}

		//--Read histogram
		results_histo[i] = dynamic_cast<TH1D*>(file->Get("Alpha"));
		if (!results_histo[i]) {
			std::cout << "ERROR: results_histo object not found... " << std::endl;
		}
		list_of_histos->Add(results_histo[i]);
	}

	//Draw separately
	TCanvas* sepCanvas = new TCanvas();
	sepCanvas->Divide(3,2);
	sepCanvas->cd(1);
	results_histo[0]->Draw();
	sepCanvas->cd(2);
	results_histo[1]->Draw();
	sepCanvas->cd(3);
	results_histo[2]->Draw();
	sepCanvas->cd(4);
	results_histo[3]->Draw();
	sepCanvas->cd(5);
	results_histo[4]->Draw();

	//Add histograms
	TH1D* add_histo = (TH1D*)results_histo[0]->Clone("add_histo");
	add_histo->Reset();
	add_histo->Merge(list_of_histos);


	//Draw together
	TCanvas* allCanvas = new TCanvas();
	add_histo->Draw();

}



