/*
 * readIRFsforLIV.C
 *
 *  Created on: Jan 18, 2017
 *      Author: Leyre Nogués
 *
 *  This code allows the extraction of the Collection area and the energy bias and resolution.
 *  The latter are computed in E, fitted with one or two Gaussians, and also in LogE, fitte with
 *  a Gaussian.
 */


#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TArrayD.h>
#include <TMath.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCollection.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TGTab.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TNtupleD.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include "TLegend.h"
#include <TFile.h>

#include "MAnalysisProblems.h"
#include "MAnalysisWizard.h"
#include "MArgs.h"
#include "MAstro.h"
#include "MAverageEnergy.h"
#include "MBinning.h"
#include "MCalcExcess.h"
#include "MContinue.h"
#include "MFillH.h"
#include "MGeomCamMagic.h"
#include "MGraphicsWizard.h"
#include "MEnergyEst.h"
#include "MEnv.h"
#include "MEvtLoop.h"
#include "MF.h"
#include "MFHadAlpha.h"
#include "MH3.h"
#include "MHadAlphaCut.h"
#include "MHadronness.h"
#include "MHAlphaEnergyTheta.h"
#include "MHEffectiveOnTime.h"
#include "MHExcessEnergyTheta.h"
#include "MHFlux.h"
#include "MHillas.h"
#include "MHMcCollectionArea.h"
#include "MHMcEnergyMigration.h"
#include "MLog.h"
#include "MLogManip.h"
#include "MMcCollectionAreaCalc.h"
#include "MParList.h"
#include "MPointingPos.h"
#include "MRawRunHeader.h"
#include "MReadMarsFile.h"
#include "MReadTree.h"
#include "MSrcPosCalc.h"
#include "MStatusDisplay.h"
#include "MStereoPar.h"
#include "MTaskList.h"
#include "MTheta2vsEest.h"
#include "MStatusArray.h"

using namespace std;

const Double_t normfac = 1./TMath::Sqrt(TMath::TwoPi());

Double_t Double_Gaus(Double_t* x, Double_t*  par){

	Double_t simgaus = (normfac/par[2])*par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	Double_t simgaus2 = (normfac/par[4])*(1-par[0])*TMath::Exp(-0.5*TMath::Power((x[0]-par[3])/par[4],2));
	return simgaus + simgaus2;
}

Double_t Simple_Gaus(Double_t* x, Double_t*  par){
	Double_t simgaus = (normfac/par[2])*par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	return simgaus;
}

void readIRFsforLIV(TString sOutputFluteName){

	Bool_t plot_results = true;

	TFile* fileOutputFlute;
	fileOutputFlute = TFile::Open(sOutputFluteName, "READ");
		if (!fileOutputFlute) {
			cout << "ERROR: file " << sOutputFluteName << " not found..."	<< endl;
		}

	//-------------------------------Create a new root file----------------------------------

	TString rootoutfile("IRFsMrk421_2014flare_Kazuma.root"); //Here name of file
	TFile* IRFoutputFile = new TFile(rootoutfile, "recreate");

	//-------------------------------Effective Area----------------------------------

	cout << "Obtaining Effective Area from Flute Output file:" << endl;
	gStyle->SetOptStat(0);

	MHMcCollectionArea* mHMcCollectionAreaEtrue =
	dynamic_cast<MHMcCollectionArea*>(fileOutputFlute->Get("MHMcCollectionAreaEtrue"));


	if (!mHMcCollectionAreaEtrue) {
		cout << "ERROR: MHMcCollectionAreaEtrue object not found... "<< endl;
	}

	const TH1D* CollArEtrue = mHMcCollectionAreaEtrue->GetHistCoarse()->ProjectionX();
	Int_t numAeffBinsX = CollArEtrue->GetNbinsX();
	const Double_t* arrayCollArBinsX = CollArEtrue->GetXaxis()->GetXbins()->GetArray();

	cout << "Number or bins in Etrue(from Flute): " << numAeffBinsX << endl;

	//We get maximum of histo to extrapolate manually for VHE
	Double_t max_of_histo = 0;
	for (Int_t i = 1; i <=numAeffBinsX; i++){
		if(max_of_histo<CollArEtrue->GetBinContent(i)) max_of_histo = CollArEtrue->GetBinContent(i);
	}
	cout << "Maximum in Collection area bins: "  << max_of_histo<< endl;
	cout << "Last bin: "  << CollArEtrue->GetBinContent(numAeffBinsX) << endl;
	cout << "Pre-last bin: "  << CollArEtrue->GetBinContent(numAeffBinsX-1) << endl;

	cout << "Obtaining TGraph from the TH1 (for other LIV codes)" << endl;
	cout << "We add some extra points to the Graph to extrapolate at VHE" << endl;


	//Filling the TGraph
	TGraph* CollArEtrue_Graph = new TGraph(numAeffBinsX);
	for (Int_t i = 0; i < numAeffBinsX; i++){
		CollArEtrue_Graph->SetPoint(i,CollArEtrue->GetBinCenter(i+1), CollArEtrue->GetBinContent(i+1));
	}


	//Give names and save everything in the root file
	CollArEtrue_Graph->SetName("CollArEtrue_Graph");
	CollArEtrue->Write();
	CollArEtrue_Graph->Write();


	if(plot_results){
		//Plot collection area resulting graphs
		TCanvas* CollArCanvas = new TCanvas();
		CollArCanvas->cd(1);
		CollArCanvas->cd(1)->SetLogx();
		CollArCanvas->cd(1)->SetLogy();
		CollArEtrue->DrawCopy();
		CollArEtrue_Graph->SetLineColor(kRed);
		CollArEtrue_Graph->Draw("same");
	}



	//-----------------------Bias and Resolution vs energy-----------------------------

	//We extrac the histograms from Flute Output and take slices in Etrue and Erec.
	//The slices are fit with two Gaussiand and the resulting parameters are saved.

	MHMcEnergyMigration *m = (MHMcEnergyMigration*)fileOutputFlute->Get("MHMcEnergyMigration");
	if(!m) {
		cout << "ERROR: MHMcEnergyMigration object not found... " << endl;
	}


	TH2D* EstTrue_True_True = (TH2D*)m->GetHist2()-> Project3D("zy"); // <-- for (E_est-E_true)/E_true vs. E_true
	TH2D* EstTrue_Est_Est = (TH2D*)m->GetHist3()-> Project3D("zx");  // <-- for (E_est-E_true)/E_true vs. E_est


	if(plot_results){
		TCanvas* resCan =  new TCanvas();
		resCan->Divide(2,1);
		resCan->cd(1);
		resCan->cd(1)->SetLogx();
		EstTrue_True_True->Draw("zcol");
		resCan->cd(2);
		resCan->cd(2)->SetLogx();
		EstTrue_Est_Est->Draw("zcol");
	}

	//Vectors to save the fitted values in Etrue and Erec.
	//Etrue
	std::vector<double> Const1_Etrue;
	std::vector<double> Mean1_Etrue;
	std::vector<double> Sigma1_Etrue;
	std::vector<double> Mean2_Etrue;
	std::vector<double> Sigma2_Etrue;
	std::vector<double> Const0_Etrue;
	std::vector<double> Mean0_Etrue;
	std::vector<double> Sigma0_Etrue;
	std::vector<double> Const0Log_Etrue;
	std::vector<double> Mean0Log_Etrue;
	std::vector<double> Sigma0Log_Etrue;



	//Erec
	std::vector<double> Const1_Erec;
	std::vector<double> Mean1_Erec;
	std::vector<double> Sigma1_Erec;
	std::vector<double> Mean2_Erec;
	std::vector<double> Sigma2_Erec;
	std::vector<double> Const0_Erec;
	std::vector<double> Mean0_Erec;
	std::vector<double> Sigma0_Erec;
	std::vector<double> Const0Log_Erec;
	std::vector<double> Mean0Log_Erec;
	std::vector<double> Sigma0Log_Erec;




	//------------------------------------------Etrue---------------------------------------------
	Int_t numberOfBinsEtrue = EstTrue_True_True->GetNbinsX();
	TH1D* ProjY_Etrue;
	TH1D* ProjY_Etrue_Log;

	TF1* dobleGaus;
	TF1* simpleGaus;
	TF1* simpleGausLog;
	Double_t cutInEDown = 100.;
	Double_t cutInEUp = 60000.;
	Int_t fit_draw = 100;

	std::vector<double> new_Etrue_array;
	Double_t last_array_value;
	Int_t first_taken_bin_Etrue = 0;



	cout << "Number of Etrue bins: " << numberOfBinsEtrue << endl;

	for(Int_t i=100; i<=100; i++){
//		cout << "Bin " << i << endl;
//		cout << "Energy True " << EstTrue_True_True->GetXaxis()->GetBinCenter(i) << endl;
		ProjY_Etrue = EstTrue_True_True->ProjectionY("",i,i, "");
		if(ProjY_Etrue->GetEntries() == 0)continue;
		if(EstTrue_True_True->GetXaxis()->GetBinLowEdge(i)>cutInEDown && EstTrue_True_True->GetXaxis()->GetBinLowEdge(i)<cutInEUp){
//			cout << "Bin accepted " << i << endl;
			if(first_taken_bin_Etrue==0){
				first_taken_bin_Etrue=i-1;
				cout << "First bin take in Etrue: " << first_taken_bin_Etrue << endl;
			}
			new_Etrue_array.push_back(EstTrue_True_True->GetXaxis()->GetBinLowEdge(i));
			last_array_value = EstTrue_True_True->GetXaxis()->GetBinLowEdge(i+1);
			const TArrayD* bins = ProjY_Etrue->GetXaxis()->GetXbins();
			TArrayD* bins_Log = new TArrayD();
//			cout << "bins array size: " << bins->GetSize() << endl;

			Int_t index_log = 0;
			for(Int_t j = 0; j<bins->GetSize(); j++){
				if((bins->GetAt(j)+1)>0){
					bins_Log->Set(index_log+1);
					bins_Log->SetAt(TMath::Log10(bins->GetAt(j)+1), index_log);
					index_log++;
				}
			}
//			cout << "bins_Log size: " << bins_Log->GetSize() << endl;
			ProjY_Etrue_Log = new TH1D(Form("ProjY_Etrue_Log%i",i), Form("ProjY_Etrue_Log%i",i), bins_Log->GetSize()-1, bins_Log->GetArray());
			Int_t first_bin = (bins->GetSize()-1) - (bins_Log->GetSize()-1);
//			cout << "First bin: " << first_bin << endl;

			for(Int_t k = 1; k<=bins_Log->GetSize()-1; k++){
				ProjY_Etrue_Log->SetBinContent(k, ProjY_Etrue->GetBinContent(first_bin));
				first_bin++;
			}

			//Fit with one Gaussian
			simpleGaus =  new TF1("Simple_Gaus", Simple_Gaus, bins->At(0), bins->At(bins->GetSize()-1), 3);
			simpleGaus->SetParameters(ProjY_Etrue->GetEntries(), 0., 0.15);
			ProjY_Etrue->Fit(simpleGaus,"0Q","",-1,1);
			Double_t Gaus0_const = ProjY_Etrue->GetFunction("Simple_Gaus")->GetParameter(0);
			Double_t Gaus0_mean = ProjY_Etrue->GetFunction("Simple_Gaus")->GetParameter(1);
			Double_t Gaus0_sigma = ProjY_Etrue->GetFunction("Simple_Gaus")->GetParameter(2);
			cout << "Gaus0: const mean sigma ->" << Gaus0_const << " " << Gaus0_mean << " " << Gaus0_sigma << endl;

			//--Save Values---
			Const0_Etrue.push_back(Gaus0_const);
			Mean0_Etrue.push_back(Gaus0_mean);
			Sigma0_Etrue.push_back(Gaus0_sigma);


			//Fit with Two Gaussians
			dobleGaus =  new TF1("Doble_Gaus", Double_Gaus, bins->At(0), bins->At(bins->GetSize()-1), 5);
			dobleGaus->SetParameters(ProjY_Etrue->GetEntries(), Gaus0_mean, Gaus0_sigma, Gaus0_mean, Gaus0_sigma*2);
			dobleGaus->SetParLimits(0,0.,ProjY_Etrue->GetEntries());
			dobleGaus->SetParLimits(1,-1.,1.);
			dobleGaus->SetParLimits(2, 0.05, 0.3);
			dobleGaus->SetParLimits(3, -1.,1.);
			dobleGaus->SetParLimits(4, 0.1,0.7);
			ProjY_Etrue->Fit(dobleGaus,"0");
			Double_t Gaus1_const = ProjY_Etrue->GetFunction("Doble_Gaus")->GetParameter(0);
			Double_t Gaus1_mean = ProjY_Etrue->GetFunction("Doble_Gaus")->GetParameter(1);
			Double_t Gaus1_sigma = ProjY_Etrue->GetFunction("Doble_Gaus")->GetParameter(2);
			Double_t Gaus2_mean = ProjY_Etrue->GetFunction("Doble_Gaus")->GetParameter(3);
			Double_t Gaus2_sigma = ProjY_Etrue->GetFunction("Doble_Gaus")->GetParameter(4);

//			if(Gaus2_sigma<Gaus1_sigma){ //Gaussian exchanged in the fit
//				Double_t temp_const = Gaus1_const;
//				Double_t temp_mean = Gaus1_mean;
//				Double_t temp_sigma = Gaus1_sigma;
//
//				Gaus1_mean = Gaus2_mean;
//				Gaus1_sigma = Gaus2_sigma;
//
//				Gaus2_mean = temp_mean;
//				Gaus2_sigma = temp_sigma;
//			}
//			dobleGaus->SetParameters(Gaus1_const, Gaus1_mean, Gaus1_sigma, Gaus2_mean, Gaus2_sigma);


			//--Save Values---
			Const1_Etrue.push_back(Gaus1_const);
			Mean1_Etrue.push_back(Gaus1_mean);
			Sigma1_Etrue.push_back(Gaus1_sigma);
			Mean2_Etrue.push_back(Gaus2_mean);
			Sigma2_Etrue.push_back(Gaus2_sigma);

			//Fit with one Gaussian in LogE
			simpleGausLog =  new TF1("SimpleGaus_Log", Simple_Gaus, bins_Log->At(0), bins_Log->At(bins_Log->GetSize()-1), 3);
			simpleGausLog->SetParameters(ProjY_Etrue_Log->GetEntries(), 0., 0.15);

			ProjY_Etrue_Log->Fit(simpleGausLog,"0Q");
			Double_t Gaus0Log_const = ProjY_Etrue_Log->GetFunction("SimpleGaus_Log")->GetParameter(0);
			Double_t Gaus0Log_mean = ProjY_Etrue_Log->GetFunction("SimpleGaus_Log")->GetParameter(1);
			Double_t Gaus0Log_sigma = ProjY_Etrue_Log->GetFunction("SimpleGaus_Log")->GetParameter(2);

			//--Save Values---
			Const0Log_Etrue.push_back(Gaus0Log_const);
			Mean0Log_Etrue.push_back(Gaus0Log_mean);
			Sigma0Log_Etrue.push_back(Gaus0Log_sigma);


			if(i == (fit_draw) && plot_results){//first_taken_bin_Etrue
				cout << "------------------ Fit results --------------------------" << endl;
				cout << "Gaus1: const mean sigma ->" << Gaus1_const << " " << Gaus1_mean << " " << Gaus1_sigma << endl;
				cout << "Gaus1: const mean sigma ->" << Gaus2_mean << " " << Gaus2_sigma << endl;
				cout << "Gaus0: const mean sigma ->" << Gaus0_const << " " << Gaus0_mean << " " << Gaus0_sigma << endl;
				cout << "Gaus0Log: const mean sigma ->" << Gaus0Log_const << " " << Gaus0Log_mean << " " << Gaus0Log_sigma << endl;

				gStyle->SetOptFit(1111);
				TH1D* copy_ProjY_Etrue = (TH1D*)ProjY_Etrue->Clone("ProjY_Etrue_draw");
				copy_ProjY_Etrue->SetStats();
				TCanvas* projection = new TCanvas();
				projection->Divide(3,1);
				projection->cd(1);
				copy_ProjY_Etrue->Draw("");
				TF1* DobleGaus_draw =  (TF1*)dobleGaus->Clone("DobleGaus_draw");
				DobleGaus_draw->Draw("same");
				projection->cd(2);
				copy_ProjY_Etrue->Draw("");
				TF1* Gaus_draw = (TF1*)simpleGaus->Clone("Gaus_draw");
				Gaus_draw->Draw("same");
				projection->cd(3);
				TH1D* copy_ProjY_Etrue_Log = (TH1D*)ProjY_Etrue_Log->Clone("ProjY_Etrue_Log_draw");
				copy_ProjY_Etrue_Log->SetStats();
				copy_ProjY_Etrue_Log->GetXaxis()->SetTitle("Log10(Erec/Etrue)");
				copy_ProjY_Etrue_Log->Draw("");
				TF1* GausLog_draw =  (TF1*)simpleGausLog->Clone("Gaus_draw");
				GausLog_draw->Draw("same");


			}
		}

	}

	cout << "Number of Entries in the saved arrays: " << Const1_Etrue.size() << endl;
	cout << "Last value in Etrue array: " << last_array_value << endl;
	new_Etrue_array.push_back(last_array_value);
	cout << "Number of components in Etrue new array: " << new_Etrue_array.size() << endl;


//	//Fill histograms with results
//	Int_t numberOfNewBins = Const1_Etrue.size();
//	Double_t* array_bins_new_Etrue = new_Etrue_array.data();
//
//
//	TH1D* Const1_histo_Etrue = new TH1D("Const1_histo_Etrue", "Const1_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Mean1_histo_Etrue = new TH1D("Mean1_histo_Etrue", "Mean1_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Sigma1_histo_Etrue = new TH1D("Sigma1_histo_Etrue", "Sigma1_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Mean2_histo_Etrue = new TH1D("Mean2_histo_Etrue", "Mean2_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Sigma2_histo_Etrue = new TH1D("Sigma2_histo_Etrue", "Sigma2_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Const0_histo_Etrue = new TH1D("Const0_histo_Etrue", "Const0_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Mean0_histo_Etrue = new TH1D("Mean0_histo_Etrue", "Mean0_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Sigma0_histo_Etrue = new TH1D("Sigma0_histo_Etrue", "Sigma0_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Const0Log_histo_Etrue = new TH1D("Const0Log_histo_Etrue", "Const0Log_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Mean0Log_histo_Etrue = new TH1D("Mean0Log_histo_Etrue", "Mean0Log_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//	TH1D* Sigma0Log_histo_Etrue = new TH1D("Sigma0Log_histo_Etrue", "Sigma0Log_histo_Etrue", numberOfNewBins, array_bins_new_Etrue);
//
//	for(Int_t i = 1 ; i<=numberOfNewBins; i++){
//		Const1_histo_Etrue->SetBinContent(i, Const1_Etrue.at(i-1));
//		Const0_histo_Etrue->SetBinContent(i, Const0_Etrue.at(i-1));
//		Const0Log_histo_Etrue->SetBinContent(i, Const0Log_Etrue.at(i-1));
//		Mean1_histo_Etrue->SetBinContent(i, Mean1_Etrue.at(i-1));
//		Mean2_histo_Etrue->SetBinContent(i, Mean2_Etrue.at(i-1));
//		Mean0_histo_Etrue->SetBinContent(i, Mean0_Etrue.at(i-1));
//		Mean0Log_histo_Etrue->SetBinContent(i, Mean0Log_Etrue.at(i-1));
//		Sigma1_histo_Etrue->SetBinContent(i, Sigma1_Etrue.at(i-1));
//		Sigma2_histo_Etrue->SetBinContent(i, Sigma2_Etrue.at(i-1));
//		Sigma0_histo_Etrue->SetBinContent(i, Sigma0_Etrue.at(i-1));
//		Sigma0Log_histo_Etrue->SetBinContent(i, Sigma0Log_Etrue.at(i-1));
//	}
//
//	cout << "Bins in final histos: " << Sigma0Log_histo_Etrue->GetNbinsX() << endl;
//	cout << "Mean in final histos: " << Sigma0Log_histo_Etrue->GetMean() << endl;
//	cout << "Content of last bin: " << Sigma0Log_histo_Etrue->GetBinContent(Sigma0Log_histo_Etrue->GetNbinsX()) << endl;
//
//	if(plot_results){
//		//Plot the results
//		TCanvas* const2canvas = new TCanvas();
//		const2canvas->SetLogx();
//		Const1_histo_Etrue->SetLineColor(kBlue);
//		Const1_histo_Etrue->Draw("");
//
//		TLegend* const2legend = new TLegend(0.1,0.7,0.4,0.9);
//		const2legend->AddEntry(Const1_histo_Etrue, "Const paratermer, Gauss 1", "l");
//		const2legend->Draw("");
//
//		TCanvas* mean2canvas = new TCanvas();
//		mean2canvas->SetLogx();
//		Mean1_histo_Etrue->SetLineColor(kBlue);
//		Mean1_histo_Etrue->Draw("");
//		Mean2_histo_Etrue->SetLineColor(kRed);
//		Mean2_histo_Etrue->Draw("same");
//
//		TLegend* main2legend = new TLegend(0.1,0.7,0.4,0.9);
//		main2legend->AddEntry(Mean1_histo_Etrue, "Mean paratermer, Gauss 1", "l");
//		main2legend->AddEntry(Mean2_histo_Etrue, "Mean paratermer, Gauss 2", "l");
//		main2legend->Draw("");
//
//		TCanvas* sigma2canvas = new TCanvas();
//		sigma2canvas->SetLogx();
//		Sigma1_histo_Etrue->SetLineColor(kBlue);
//		Sigma1_histo_Etrue->Draw("");
//		Sigma2_histo_Etrue->SetLineColor(kRed);
//		Sigma2_histo_Etrue->Draw("same");
//
//		TLegend* sigma2legend = new TLegend(0.1,0.7,0.4,0.9);
//		sigma2legend->AddEntry(Sigma1_histo_Etrue, "Sigma paratermer, Gauss 1", "l");
//		sigma2legend->AddEntry(Sigma2_histo_Etrue, "Sigma paratermer, Gauss 2", "l");
//		sigma2legend->Draw("");
//	}

	//Save Results
//	Const1_histo_Etrue->Write();
//	Mean1_histo_Etrue->Write();
//	Mean2_histo_Etrue->Write();
//	Sigma1_histo_Etrue->Write();
//	Sigma2_histo_Etrue->Write();


/*
	//------------------------------------------Erec---------------------------------------------
	Int_t numberOfBinsErec = EstTrue_Est_Est->GetNbinsX();
	TH1D* ProjY_Erec;
	TH1D* ProjY_Erec_Log;

	cutInEDown = 0.;
	cutInEUp = 70000.;
	fit_draw = 28;

	std::vector<double> new_Erec_array;
	Int_t first_taken_bin_Erec = 0;


	cout << "Number of Erec bins: " << numberOfBinsErec << endl;

	for(Int_t i=1; i<=numberOfBinsErec; i++){
//		cout << "Bin " << i << endl;
//		cout << "Energy Rec " << EstTrue_True_Est->GetXaxis()->GetBinCenter(i) << endl;
		ProjY_Erec = EstTrue_Est_Est->ProjectionY("",i,i, "");
		if(ProjY_Erec->GetEntries() == 0)continue;
		if(EstTrue_Est_Est->GetXaxis()->GetBinLowEdge(i)>cutInEDown && EstTrue_Est_Est->GetXaxis()->GetBinLowEdge(i)<cutInEUp){
//			cout << "Bin accepted " << i << endl;
			if(first_taken_bin_Erec==0){
				first_taken_bin_Erec=i-1;
				cout << "First bin take in Erec: " << first_taken_bin_Erec << endl;
			}
			new_Erec_array.push_back(EstTrue_Est_Est->GetXaxis()->GetBinLowEdge(i));
			last_array_value = EstTrue_Est_Est->GetXaxis()->GetBinLowEdge(i+1);
			const TArrayD* bins = ProjY_Erec->GetXaxis()->GetXbins();
			TArrayD* bins_Log = new TArrayD();
//			cout << "bins array size: " << bins->GetSize() << endl;

			Int_t index_log = 0;
			for(Int_t j = 0; j<bins->GetSize(); j++){
//				cout << "Bin bins: " << bins->GetAt(j)+1 << endl;
				if((bins->GetAt(j)+1)>0){
					bins_Log->Set(index_log+1);
//					cout << "bins_Log length: " << bins_Log->GetSize() << endl;
					bins_Log->SetAt(TMath::Log10(bins->GetAt(j)+1), index_log);
//					cout << "Before: " << bins->GetAt(j)+1 << endl;
//					cout << "After: " << TMath::Log10(bins->GetAt(j)+1)<< endl;
//					cout << "bins_Log component: " << bins_Log->At(index_log) << endl;
					index_log++;
				}
			}
//			cout << "bins_Log length: " << bins_Log->GetSize() << endl;
			ProjY_Erec_Log = new TH1D(Form("ProjY_Erec_Log%i",i), Form("ProjY_Erec_Log%i",i), bins_Log->GetSize()-1, bins_Log->GetArray());
			Int_t first_bin = (bins->GetSize()-1) - (bins_Log->GetSize()-1);
//			cout << "First bin: " << first_bin << endl;

			for(Int_t k = 1; k<=bins_Log->GetSize()-1; k++){
				ProjY_Erec_Log->SetBinContent(k, ProjY_Erec->GetBinContent(first_bin));
				first_bin++;
			}

			//Normalize histograms before fit
			ProjY_Erec->Scale(1./ProjY_Erec->Integral());
			ProjY_Erec_Log->Scale(1./ProjY_Erec_Log->Integral());

			//Fit with one Gaussian
			ProjY_Erec->Fit("gaus","Q0");
			Double_t Gaus0_const = ProjY_Erec->GetFunction("gaus")->GetParameter(0);
			Double_t Gaus0_mean = ProjY_Erec->GetFunction("gaus")->GetParameter(1);
			Double_t Gaus0_sigma = ProjY_Erec->GetFunction("gaus")->GetParameter(2);
			//--Save Values---
			Const0_Erec.push_back(Gaus0_const);
			Mean0_Erec.push_back(Gaus0_mean);
			Sigma0_Erec.push_back(Gaus0_sigma);


			//Fit with Two Gaussians
			dobleGaus =  new TF1("Doble Gauss", Double_Gaus, bins->At(0), bins->At(bins->GetSize()-1), 6);
			dobleGaus->SetParameters(Gaus0_const, Gaus0_mean, Gaus0_sigma, Gaus0_const/2., Gaus0_mean, Gaus0_sigma*2);
			ProjY_Erec->Fit(dobleGaus,"Q0");
			Double_t Gaus1_const = ProjY_Erec->GetFunction("Doble Gauss")->GetParameter(0);
			Double_t Gaus1_mean = ProjY_Erec->GetFunction("Doble Gauss")->GetParameter(1);
			Double_t Gaus1_sigma = ProjY_Erec->GetFunction("Doble Gauss")->GetParameter(2);
			Double_t Gaus2_mean = ProjY_Erec->GetFunction("Doble Gauss")->GetParameter(4);
			Double_t Gaus2_sigma = ProjY_Erec->GetFunction("Doble Gauss")->GetParameter(5);

			if(Gaus2_sigma<Gaus1_sigma){ //Gaussian exchanged in the fit
				Double_t temp_mean = Gaus1_mean;
				Double_t temp_sigma = Gaus1_sigma;

				Gaus1_mean = Gaus2_mean;
				Gaus1_sigma = Gaus2_sigma;

				Gaus2_mean = temp_mean;
				Gaus2_sigma = temp_sigma;
			}
			//--Save Values---
			Const1_Erec.push_back(Gaus1_const);
			Mean1_Erec.push_back(Gaus1_mean);
			Sigma1_Erec.push_back(Gaus1_sigma);
			Mean2_Erec.push_back(Gaus2_mean);
			Sigma2_Erec.push_back(Gaus2_sigma);

			//Fit with one Gaussian in LogE
			ProjY_Erec_Log->Fit("gaus","Q0");
			Double_t Gaus0Log_const = ProjY_Erec_Log->GetFunction("gaus")->GetParameter(0);
			Double_t Gaus0Log_mean = ProjY_Erec_Log->GetFunction("gaus")->GetParameter(1);
			Double_t Gaus0Log_sigma = ProjY_Erec_Log->GetFunction("gaus")->GetParameter(2);
			//--Save Values---
			Const0Log_Erec.push_back(Gaus0Log_const);
			Mean0Log_Erec.push_back(Gaus0Log_mean);
			Sigma0Log_Erec.push_back(Gaus0Log_sigma);


			if(i == (fit_draw+first_taken_bin_Erec) && plot_results){
				cout << "------------------ Fit results --------------------------" << endl;
				cout << "Gaus1: const mean sigma ->" << Gaus1_const << " " << Gaus1_mean << " " << Gaus1_sigma << endl;

				gStyle->SetOptFit(1111);
				TH1D* copy_ProjY_Erec = (TH1D*)ProjY_Erec->Clone("ProjY_Erec_draw");
				copy_ProjY_Erec->SetStats();
				TCanvas* projection = new TCanvas();
				projection->Divide(3,1);
				projection->cd(1);
				copy_ProjY_Erec->Draw("");
				TF1* DobleGaus_draw =  new TF1("DobleGaus_draw", Double_Gaus, bins->At(0), bins->At(bins->GetSize()-1), 6);
				DobleGaus_draw->SetParameters(Gaus1_const, Gaus1_mean, Gaus1_sigma, Gaus2_mean, Gaus2_sigma);
				DobleGaus_draw->Draw("same");
				projection->cd(2);
				copy_ProjY_Erec->Draw("");
				TF1* Gaus_draw =  new TF1("Gaus_draw", Simple_Gaus, bins->At(0), bins->At(bins->GetSize()-1), 6);
				Gaus_draw->SetParameters(Gaus0_const, Gaus0_mean, Gaus0_sigma);
				Gaus_draw->Draw("same");
				projection->cd(3);
				TH1D* copy_ProjY_Erec_Log = (TH1D*)ProjY_Erec_Log->Clone("ProjY_Erec_Log_draw");
				copy_ProjY_Erec_Log->SetStats();
				copy_ProjY_Erec_Log->GetXaxis()->SetTitle("Log10(Erec/Etrue)");
				copy_ProjY_Erec_Log->Draw("");
				TF1* GausLog_draw =  new TF1("GausLog_draw", Simple_Gaus, bins->At(0), bins->At(bins->GetSize()-1), 6);
				GausLog_draw->SetParameters(Gaus0Log_const, Gaus0Log_mean, Gaus0Log_sigma);
				GausLog_draw->Draw("same");

			}
		}
	}

	cout << "Number of Entries in the saved arrays: " << Const1_Erec.size() << endl;
	cout << "Last value in Etrue array: " << last_array_value << endl;
	new_Erec_array.push_back(last_array_value);
	cout << "Number of components in Erec new array: " << new_Erec_array.size() << endl;


	//Fill histograms with results
	Int_t numberOfNewBins_rec = Const1_Erec.size();
	Double_t* array_bins_new_Erec = new_Erec_array.data();


	TH1D* Const1_histo_Erec = new TH1D("Const1_histo_Erec", "Const1_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Mean1_histo_Erec = new TH1D("Mean1_histo_Erec", "Mean1_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Sigma1_histo_Erec = new TH1D("Sigma1_histo_Erec", "Sigma1_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Mean2_histo_Erec = new TH1D("Mean2_histo_Erec", "Mean2_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Sigma2_histo_Erec = new TH1D("Sigma2_histo_Erec", "Sigma2_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Const0_histo_Erec = new TH1D("Const0_histo_Erec", "Const0_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Mean0_histo_Erec = new TH1D("Mean0_histo_Erec", "Mean0_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Sigma0_histo_Erec = new TH1D("Sigma0_histo_Erec", "Sigma0_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Const0Log_histo_Erec = new TH1D("Const0Log_histo_Erec", "Const0Log_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Mean0Log_histo_Erec = new TH1D("Mean0Log_histo_Erec", "Mean0Log_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);
	TH1D* Sigma0Log_histo_Erec = new TH1D("Sigma0Log_histo_Erec", "Sigma0Log_histo_Erec", numberOfNewBins_rec, array_bins_new_Erec);

	for(Int_t i = 1 ; i<=numberOfNewBins_rec; i++){
		Const1_histo_Erec->SetBinContent(i, Const1_Erec.at(i-1));
		Const0_histo_Erec->SetBinContent(i, Const0_Erec.at(i-1));
		Const0Log_histo_Erec->SetBinContent(i, Const0Log_Erec.at(i-1));
		Mean1_histo_Erec->SetBinContent(i, Mean1_Erec.at(i-1));
		Mean2_histo_Erec->SetBinContent(i, Mean2_Erec.at(i-1));
		Mean0_histo_Erec->SetBinContent(i, Mean0_Erec.at(i-1));
		Mean0Log_histo_Erec->SetBinContent(i, Mean0Log_Erec.at(i-1));
		Sigma1_histo_Erec->SetBinContent(i, Sigma1_Erec.at(i-1));
		Sigma2_histo_Erec->SetBinContent(i, Sigma2_Erec.at(i-1));
		Sigma0_histo_Erec->SetBinContent(i, Sigma0_Erec.at(i-1));
		Sigma0Log_histo_Erec->SetBinContent(i, Sigma0Log_Erec.at(i-1));
	}

	cout << "Bins in final histos: " << Sigma0Log_histo_Erec->GetNbinsX() << endl;
	cout << "Mean in final histos: " << Sigma0Log_histo_Erec->GetMean() << endl;
	cout << "Content of last bin: " << Sigma0Log_histo_Erec->GetBinContent(Sigma0Log_histo_Erec->GetNbinsX()) << endl;

	if(plot_results){
		//Plot the results
		TCanvas* const2canvas = new TCanvas();
		const2canvas->SetLogx();
		Const1_histo_Erec->SetLineColor(kBlue);
		Const1_histo_Erec->Draw("");

		TLegend* const2legend = new TLegend(0.1,0.7,0.4,0.9);
		const2legend->AddEntry(Const1_histo_Erec, "Const paratermer, Gauss 1", "l");
		const2legend->Draw("");

		TCanvas* mean2canvas = new TCanvas();
		mean2canvas->SetLogx();
		Mean1_histo_Erec->SetLineColor(kBlue);
		Mean1_histo_Erec->Draw("");
		Mean2_histo_Erec->SetLineColor(kRed);
		Mean2_histo_Erec->Draw("same");

		TLegend* main2legend = new TLegend(0.1,0.7,0.4,0.9);
		main2legend->AddEntry(Mean1_histo_Erec, "Mean paratermer, Gauss 1", "l");
		main2legend->AddEntry(Mean2_histo_Erec, "Mean paratermer, Gauss 2", "l");
		main2legend->Draw("");

		TCanvas* sigma2canvas = new TCanvas();
		sigma2canvas->SetLogx();
		Sigma1_histo_Erec->SetLineColor(kBlue);
		Sigma1_histo_Erec->Draw("");
		Sigma2_histo_Erec->SetLineColor(kRed);
		Sigma2_histo_Erec->Draw("same");

		TLegend* sigma2legend = new TLegend(0.1,0.7,0.4,0.9);
		sigma2legend->AddEntry(Sigma1_histo_Erec, "Sigma paratermer, Gauss 1", "l");
		sigma2legend->AddEntry(Sigma2_histo_Erec, "Sigma paratermer, Gauss 2", "l");
		sigma2legend->Draw("");
	}


	//Save Results
	Const1_histo_Erec->Write();
	Mean1_histo_Erec->Write();
	Mean2_histo_Erec->Write();
	Sigma1_histo_Erec->Write();
	Sigma2_histo_Erec->Write();
*/

}



