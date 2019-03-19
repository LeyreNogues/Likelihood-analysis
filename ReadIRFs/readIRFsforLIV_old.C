/*
 * readIRFsforLIV.C
 *
 *  Created on: Jan 18, 2017
 *      Author: Leyre Nogués
 *
 *  This code allows the extraction of the Migration Matrix and Collection Area from flute files.
 *  It allow also the conversion of both of them to LogE.
 *  It also extracts directly from flute the histogram to compute the energy bias a resolution
 *  and converts them to LogE before doing a Gaussian Fit.
 */


#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TArrayD.h>
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

Double_t Double_Gaus(Double_t* x, Double_t*  par){
	Double_t simgaus = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	Double_t simgaus2 = par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5],2));
	return simgaus + simgaus2;
}

void readIRFsforLIV(TString sOutputFluteName){

	TFile* fileOutputFlute;
	fileOutputFlute = TFile::Open(sOutputFluteName, "READ");
		if (!fileOutputFlute) {
			cout << "ERROR: file " << sOutputFluteName << " not found..."	<< endl;
		}

	//-------------------------------Create a new root file----------------------------------

	TString rootoutfile("IRFsMrk421_2014flare_Kazuma_LogE.root"); //Here name of file
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


	//-----Convert Effective Area to Log Scale-----
	TH1D *CollAreaLog10Etrue = new TH1D("EffectiveAreaEtrue", "Effective Area;Log_{10}(E' [GeV]);Collection Area[m^{2}]",numAeffBinsX, TMath::Log10(arrayCollArBinsX[0]),TMath::Log10(arrayCollArBinsX[numAeffBinsX]));
	for(Int_t i = 1; i < numAeffBinsX + 1; i++){
		CollAreaLog10Etrue->SetBinContent(i, CollArEtrue->GetBinContent(i));
		CollAreaLog10Etrue->SetBinError(i, CollArEtrue->GetBinError(i));
	}

	//Create a TGraph from the histogram
	TGraph* CollAreaLog10Etrue_Graph = new TGraph(numAeffBinsX);
	for (Int_t i = 0; i < numAeffBinsX; i++){
		CollAreaLog10Etrue_Graph->SetPoint(i,CollAreaLog10Etrue->GetBinCenter(i+1), CollAreaLog10Etrue->GetBinContent(i+1));
		}


	//Give names and save everything in the root file
	CollArEtrue_Graph->SetName("CollArEtrue_Graph");
	CollAreaLog10Etrue->SetName("CollAreaLog10Etrue");
	CollAreaLog10Etrue_Graph->SetName("CollAreaLog10Etrue_Graph");
	CollArEtrue->Write();
	CollArEtrue_Graph->Write();
	CollAreaLog10Etrue->Write();
	CollAreaLog10Etrue_Graph->Write();


	//Plot collection area resulting graphs
	TCanvas* CollArCanvas = new TCanvas();
	CollArCanvas->Divide(2,1);
	CollArCanvas->cd(1);
	CollArCanvas->cd(1)->SetLogx();
	CollArCanvas->cd(1)->SetLogy();
	CollArEtrue->DrawCopy();
	CollArEtrue_Graph->SetLineColor(kRed);
	CollArEtrue_Graph->Draw("same");
	CollArCanvas->cd(2);
	CollArCanvas->cd(2)->SetLogy();
	CollAreaLog10Etrue->DrawCopy();
	CollAreaLog10Etrue_Graph->SetLineColor(kRed);
	CollAreaLog10Etrue_Graph->Draw("same");

//	CollArCanvas->SaveAs("Flute_CollAr.png");
//	CollArCanvas->SaveAs("Flute_CollAr.pdf");



	//-------------------------------Migration Matrix----------------------------------

	cout << "Obtaining Migration Matrix from Flute Output file:" << endl;

	TH3D* MigMatrix3D = dynamic_cast<TH3D*>(fileOutputFlute->Get("migrmatrix"));
	if(!MigMatrix3D) {
		cout << "ERROR: MigMatrix3D object not found... " << endl;
	}


	TH2D* MigMatrix = (TH2D*) MigMatrix3D->Project3D("yx");
	delete MigMatrix3D;


	Int_t numMigBinsX = MigMatrix->GetNbinsX();
	Int_t numMigBinsY = MigMatrix->GetNbinsY();
	cout << "Bins in Erec: " << numMigBinsX << endl;
	cout << "Bins in Etrue: " << numMigBinsY << endl;
	const Double_t* arrayMigBinsX =	MigMatrix->GetXaxis()->GetXbins()->GetArray();
	const Double_t* arrayMigBinsY = MigMatrix->GetYaxis()->GetXbins()->GetArray();

	//-----Convert Migration matrix to Log Scale (Better to do gaussian fit)-----

	TH2D *LogMigMatrix = new TH2D("LogMigMatrix", ";Log_{10}(E'[GeV]); Log_{10}(E[GeV])", numMigBinsX, TMath::Log10(arrayMigBinsX[0]),TMath::Log10(arrayMigBinsX[numMigBinsX]), numMigBinsY,TMath::Log10(arrayMigBinsY[0]),TMath::Log10(arrayMigBinsY[numMigBinsY]));
	for (Int_t i = 1; i < numMigBinsX + 1; i++) {
		for (Int_t j = 1; j < numMigBinsY + 1; j++) {
			LogMigMatrix->SetBinContent(i, j,
					MigMatrix->GetBinContent(i, j));
			LogMigMatrix->SetBinError(i, j, MigMatrix->GetBinError(i, j));
			}
	}

	//Plot resulting graphs
	TCanvas* migmatcanvas = new TCanvas();
	migmatcanvas->Divide(2,1);
	migmatcanvas->cd(1);
	migmatcanvas->cd(1)->SetLogx();
	migmatcanvas->cd(1)->SetLogy();
	MigMatrix->DrawCopy("zcol");
	migmatcanvas->cd(2);
	LogMigMatrix->DrawCopy("zcol");
//	migmatcanvas->SaveAs("Flute_MigMatrix.png");
//	migmatcanvas->SaveAs("Flute_MigMatrix.pdf");



	//Give names to save everything in the root file
	MigMatrix->SetName("MigMatrix");
	LogMigMatrix->SetName("LogMigMatrix");
	MigMatrix->Write();
	LogMigMatrix->Write();



	//-----------------------Bias and Resolution vs energy-----------------------------

	gStyle->SetOptFit(1111);
	Int_t smooth_times = 30;

	//-----Options for bias and res-----
	Bool_t fit_gaus = true; //Fit or not the bins with gaussians.
	if (fit_gaus == true) cout << "I will fit gaussians" << endl;
	else cout << "I will take directly the histograms" << endl;

	//--Draw a bin fit-- (Only true if fit_gaus == true)
	Bool_t draw_fit = true;
	Int_t draw_fit_number = 17;


	//Cut the axis arrays
	Int_t numMigBinsY_cut = 0;
	Int_t numMigBinsX_cut = 0;
	Double_t Emin_cut = 100.;
	Double_t Emax_cut = 30000;

	//--Y
	Int_t binY_Emin = MigMatrix->GetYaxis()->FindBin(Emin_cut);
	Int_t binY_Emax = MigMatrix->GetYaxis()->FindBin(Emax_cut);
	cout << "The first bin is " << binY_Emin << " and the last " << binY_Emax << endl;
	numMigBinsY_cut = (binY_Emax-binY_Emin)+1;
	Int_t indexY_new = (binY_Emax-binY_Emin)+1;
	const Int_t index_newY2 = indexY_new+1;
	Double_t arrayMigBinsY_cut[index_newY2];
	cout << "The new histo will have " << numMigBinsY_cut << " Y bins and an array of " << indexY_new << " components." << endl;

	for(Int_t i = 0; i<=indexY_new; i++){
		Int_t index = (binY_Emin-1)+1*i;
		Double_t value = *(arrayMigBinsY+index);
		arrayMigBinsY_cut[i]= value;
	}

	//--X
	Int_t binX_Emin = MigMatrix->GetXaxis()->FindBin(Emin_cut);
	Int_t binX_Emax = MigMatrix->GetXaxis()->FindBin(Emax_cut);
	cout << "The first bin is " << binX_Emin << " and the last " << binX_Emax << endl;
	numMigBinsX_cut = (binX_Emax-binX_Emin)+1;
	Int_t indexX_new = (binX_Emax-binX_Emin)+1;
	const Int_t index_newX2 = indexX_new+1;
	Double_t arrayMigBinsX_cut[index_newX2];
	cout << "The new histo will have " << numMigBinsX_cut << " X bins and an array of " << indexX_new << " components." << endl;

	for(Int_t i = 0; i<=indexX_new; i++){
		Int_t index = (binX_Emin-1)+1*i;
		Double_t value = *(arrayMigBinsX+index);
		arrayMigBinsX_cut[i]= value;
	}


	//-----Plots vs Etrue-----
	cout << "Here I compute Bias and Resolution vs Etrue: for Simulation" << endl;
	cout << "For a given Etrue, I proyect in y-axis and get a histo vs Erec" << endl;

	TH1D* Etrue_projX;
	TF1* Fit_projX;
	Double_t Etrue;
	TH1D* BiasvsEtrue = new TH1D("BiasvsEtrue", "BiasvsEtrue", numMigBinsY,arrayMigBinsY);
	TH1D* ResvsEtrue = new TH1D("ResolutionVsEtrue", "ResolutionVsEtrue", numMigBinsY,arrayMigBinsY);
	TH1D* BiasvsEtrue_Cut = new TH1D("BiasvsEtrue_cut", "BiasvsEtrue_cut", numMigBinsY_cut,arrayMigBinsY_cut);
	TH1D* ResvsEtrue_Cut = new TH1D("ResvsEtrue_cut", "ResvsEtrue_cut", numMigBinsY_cut,arrayMigBinsY_cut);

	TH1D* FitvsEtrue;
	TF1* FitGausEtrue;


	for (Int_t i = 1; i <=numMigBinsY; i++) {
		Etrue = TMath::Power(10,LogMigMatrix->GetYaxis()->GetBinCenter(i));
//		cout << "Etrue: " << i << " " << Etrue << endl;
		Etrue_projX = MigMatrix->ProjectionX("Proj in Erec",i,i, "");

		if(Etrue_projX->GetMean(1)>0){ //Not-empty bin
			if(fit_gaus==true){
				Etrue_projX->Fit("gaus","Q0");
				Fit_projX = Etrue_projX->GetFunction("gaus");
				if(draw_fit == true && i == draw_fit_number){
					FitvsEtrue = (TH1D*)Etrue_projX->Clone("FitvsEtrue");
					FitGausEtrue = (TF1*)Fit_projX->Clone("FitGaus");
				}
				Double_t Erec_fit = Fit_projX->GetParameter(1);
				Double_t Sigma_Erec_fit = Fit_projX->GetParameter(2);
//				cout << "Erec_fit: " << Erec_fit << "  Sigma_Erec_fit: " << Sigma_Erec_fit << endl;
				Double_t bias= Erec_fit/Etrue;
				Double_t res= Sigma_Erec_fit/Etrue;
//				cout << "Bias: " << bias << "  Resolution: " << res << endl;
				BiasvsEtrue->SetBinContent(i, bias);
				ResvsEtrue->SetBinContent(i, res);
				if(i>=binY_Emin && i<=binY_Emax){
					BiasvsEtrue_Cut->SetBinContent(i-binY_Emin+1, bias);
					ResvsEtrue_Cut->SetBinContent(i-binY_Emin+1, res);
				}

			}
			else{
				Double_t Erec_histo = Etrue_projX->GetMean(1);
				Double_t Sigma_Erec_histo = Etrue_projX->GetRMS(1);
//				cout << "Erec_histo: " << Erec_histo << "  Sigma_Erec_histo: " << Sigma_Erec_histo << endl;
				Double_t bias= Erec_histo/Etrue;
				Double_t res= Sigma_Erec_histo/Etrue;
//				cout << "Bias: " << bias << "  Resolution: " << res << endl;
				BiasvsEtrue->SetBinContent(i, bias);
				ResvsEtrue->SetBinContent(i, res);
				if(i>=binY_Emin && i<=binY_Emax){
					BiasvsEtrue_Cut->SetBinContent(i-binY_Emin+1, bias);
					ResvsEtrue_Cut->SetBinContent(i-binY_Emin+1, res);
				}
			}
		}
	}

	//---Plot resulting graphs---
	TCanvas* ResandBiasEtrue_canvas = new TCanvas();
	ResandBiasEtrue_canvas->SetTitle("Resolution and Bias vs Etrue");
	ResandBiasEtrue_canvas->Divide(2,1);
	ResandBiasEtrue_canvas->cd(1);
	ResandBiasEtrue_canvas->cd(1)->SetLogx();
	BiasvsEtrue_Cut->Smooth(smooth_times);
	BiasvsEtrue_Cut->GetYaxis()->SetTitle("Erec/Etrue");
	BiasvsEtrue_Cut->GetXaxis()->SetTitle("Etrue(GeV)");
	BiasvsEtrue_Cut->SetLineColor(kBlue);
	BiasvsEtrue_Cut->SetLineWidth(2);
//	BiasvsEtrue->Draw("");
	BiasvsEtrue_Cut->Draw("");
	ResandBiasEtrue_canvas->cd(2);
	ResandBiasEtrue_canvas->cd(2)->SetLogx();
	ResvsEtrue_Cut->Smooth(smooth_times);
	ResvsEtrue_Cut->GetYaxis()->SetTitle("#sigma/Etrue");
	ResvsEtrue_Cut->GetXaxis()->SetTitle("Etrue(GeV)");
	ResvsEtrue_Cut->SetLineWidth(2);
	ResvsEtrue_Cut->SetLineColor(kRed);
//	ResvsEtrue->Draw("");
	ResvsEtrue_Cut->Draw("");

	if(draw_fit == true){
		TCanvas* fit_canvas = new TCanvas();
		FitvsEtrue->Draw();
		FitGausEtrue->Draw("same");
	}


	//Give names and save everything in the root file
//	BiasvsEtrue->SetName("BiasVsEtrue");
//	ResvsEtrue->SetName("ResolutionVsEtrue");
//	BiasvsEtrue->Write();
//	ResvsEtrue->Write();

//	cout << "Bias and Resolution vs Etrue computed and saved" << endl;


	//-----Plots vs Erec-----

	cout << "Here I compute Bias and Resolution vs Erec: for Likelihood" << endl;
	cout << "For a given Erec, I proyect in x-axis and get a histo vs Etrue" << endl;

	TH1D* Erec_projY;
	TF1* Fit_projY;
	Double_t Erec;
	TH1D* BiasvsErec = new TH1D("BiasVsErec","BiasVsErec", numMigBinsX, arrayMigBinsX);
	TH1D* ResvsErec = new TH1D("ResolutionVsErec","ResolutionVsErec", numMigBinsX, arrayMigBinsX);
	TH1D* BiasvsErec_Cut = new TH1D("BiasvsErec_cut", "BiasvsErec_cut", numMigBinsX_cut,arrayMigBinsX_cut);
	TH1D* ResvsErec_Cut = new TH1D("ResvsErec_cut", "ResvsErec_cut", numMigBinsX_cut,arrayMigBinsX_cut);

	TH1D* FitvsErec;
	TF1* FitGausErec;

	for (Int_t i = 1; i <=numMigBinsX; i++) {
		Erec = TMath::Power(10,LogMigMatrix->GetXaxis()->GetBinCenter(i));
//		cout << "Erec: " << i << " " << Erec << endl;
		Erec_projY = MigMatrix->ProjectionY("Proj in Etrue",i,i, "");

		if(Erec_projY->GetMean(1)>0){ //Not-empty bins
			if(fit_gaus == true){
				Erec_projY->Fit("gaus","Q0");
				Fit_projY = Erec_projY->GetFunction("gaus");
				if(draw_fit == true && i == draw_fit_number){
					FitvsErec = (TH1D*)Erec_projY->Clone("FitvsErec");
					FitGausErec = (TF1*)Fit_projY->Clone("FitGausErec");
				}
				Double_t Etrue_fit = Fit_projY->GetParameter(1);
				Double_t Sigma_Etrue_fit = Fit_projY->GetParameter(2);
//				cout << "Etrue_fit: " << Etrue_fit << "  Sigma_Etrue_fit: " << Sigma_Etrue_fit << endl;
				Double_t bias= Etrue_fit/Erec;
				Double_t res= Sigma_Etrue_fit/Erec;
//				cout << "Bias: " << bias << "  Resolution: " << res << endl;
				BiasvsErec->SetBinContent(i, bias);
				ResvsErec->SetBinContent(i, res);
				if(i>=binX_Emin && i<=binX_Emax){
					BiasvsErec_Cut->SetBinContent(i-binX_Emin+1, bias);
					ResvsErec_Cut->SetBinContent(i-binX_Emin+1, res);
				}
			}
			else{
				Double_t Etrue_histo = Erec_projY->GetMean(1);
				Double_t Sigma_Etrue_histo = Erec_projY->GetRMS(1);
//				cout << "Etrue_histo: " << Etrue_histo << "  Sigma_Etrue_histo: " << Sigma_Etrue_histo << endl;
				Double_t bias= Etrue_histo/Erec;
				Double_t res= Sigma_Etrue_histo/Erec;
//				cout << "Bias: " << bias << "  Resolution: " << res << endl;
				BiasvsErec->SetBinContent(i, bias);
				ResvsErec->SetBinContent(i, res);
				if(i>=binX_Emin && i<=binX_Emax){
					BiasvsErec_Cut->SetBinContent(i-binX_Emin+1, bias);
					ResvsErec_Cut->SetBinContent(i-binX_Emin+1, res);
				}

			}
		}
	}

	//---Plot resulting graphs---
	TCanvas* ResandBiasErec_canvas = new TCanvas();
	ResandBiasErec_canvas->SetTitle("Resolution and Bias vs Erec");
	ResandBiasErec_canvas->Divide(2,1);
	ResandBiasErec_canvas->cd(1);
	ResandBiasErec_canvas->cd(1)->SetLogx();
	BiasvsErec_Cut->Smooth(smooth_times);
	BiasvsErec_Cut->GetYaxis()->SetTitle("Etrue/Erec");
	BiasvsErec_Cut->GetXaxis()->SetTitle("Erec(GeV)");
	BiasvsErec_Cut->SetLineColor(kBlue);
	BiasvsErec_Cut->SetLineWidth(2);
	BiasvsErec_Cut->Draw("");
	ResandBiasErec_canvas->cd(2);
	ResandBiasErec_canvas->cd(2)->SetLogx();
	ResvsErec_Cut->Smooth(smooth_times);
	ResvsErec_Cut->GetYaxis()->SetTitle("#sigma/Erec");
	ResvsErec_Cut->GetXaxis()->SetTitle("Erec(GeV)");
	ResvsErec_Cut->SetLineWidth(2);
	ResvsErec_Cut->SetLineColor(kRed);
	ResvsErec_Cut->Draw("");

	if(draw_fit == true){
		TCanvas* fit_canvas = new TCanvas();
		FitvsErec->Draw();
		FitGausErec->Draw("same");
	}

	//Give names and save everything in the root file
	BiasvsErec->SetName("BiasVsErec");
	ResvsErec->SetName("ResolutionVsErec");
	BiasvsErec->Write();
	ResvsErec->Write();

	cout << "Bias and Resolution vs Erec computed and saved" << endl;



	//-----------------------Bias and Resolution vs log energy-----------------------------


	//-----Plots vs LogEtrue-----

	cout << "Here I compute Bias and Resolution vs LogEtrue: for Simulation" << endl;
	cout << "For a given LogEtrue, I proyect in y-axis and get a histo vs LogErec" << endl;
	cout << "For fit case, the bias and resultion are defined in LogE" << endl;

	TH1D* LogEtrue_projX;
	TF1* LogFit_projX;
	Double_t LogEtrue;

	TH1D* BiasvslogEtrue = new TH1D("BiasVsLogEtrue","BiasVsLogEtrue", numMigBinsY,TMath::Log10(arrayMigBinsY[0]),TMath::Log10(arrayMigBinsY[numMigBinsY]));
	TH1D* ResvslogEtrue = new TH1D("ResolutionVsLogEtrue","ResolutionVsLogEtrue", numMigBinsY,TMath::Log10(arrayMigBinsY[0]),TMath::Log10(arrayMigBinsY[numMigBinsY]));
	TH1D* BiasvslogEtrue_cut = new TH1D("BiasvslogEtrue_cut","BiasvslogEtrue_cut", numMigBinsY_cut,TMath::Log10(arrayMigBinsY_cut[0]),TMath::Log10(arrayMigBinsY_cut[numMigBinsY_cut]));
	TH1D* ResvslogEtrue_cut = new TH1D("ResvslogEtrue_cut","ResvslogEtrue_cut", numMigBinsY_cut,TMath::Log10(arrayMigBinsY_cut[0]),TMath::Log10(arrayMigBinsY_cut[numMigBinsY_cut]));

	TH1D* FitvsLogEtrue;
	TF1* FitGausLogEtrue;

	for (Int_t i = 1; i <=numMigBinsY; i++) {
		LogEtrue = LogMigMatrix->GetYaxis()->GetBinCenter(i);
//		cout << "LogEtrue: " << i << " " << LogEtrue << endl;
		LogEtrue_projX = LogMigMatrix->ProjectionX("Proj in logErec",i,i, "");
		if(LogEtrue_projX->GetMean(1)>0){//This assures bins are not empty
			if(fit_gaus==true){
				LogEtrue_projX->Fit("gaus","Q0");
				LogFit_projX = LogEtrue_projX->GetFunction("gaus");
				if(draw_fit == true && i == draw_fit_number){
					FitvsLogEtrue = (TH1D*)LogEtrue_projX->Clone("FitvsLogEtrue");
					FitGausLogEtrue = (TF1*)LogFit_projX->Clone("FitGausLogEtrue");
				}
				Double_t LogErec_fit = LogFit_projX->GetParameter(1);
				Double_t Sigma_LogErec_fit = LogFit_projX->GetParameter(2);
//				cout << "LogErec_fit: " << LogErec_fit << "  Sigma_LogErec_fit: " << Sigma_LogErec_fit << endl;
				Double_t bias = LogErec_fit/LogEtrue; //Is it ok?
				Double_t res= Sigma_LogErec_fit;
//				cout << "Bias: " << bias << "  Resolution: " << res << endl;
				BiasvslogEtrue->SetBinContent(i, bias);
				ResvslogEtrue->SetBinContent(i, res);
				if(i>=binY_Emin && i<=binY_Emax){
					BiasvslogEtrue_cut->SetBinContent(i-binY_Emin+1, bias);
					ResvslogEtrue_cut->SetBinContent(i-binY_Emin+1, res);
				}

			}
			else{
				Double_t LogErec_histo = LogEtrue_projX->GetMean(1);
				Double_t Sigma_LogErec_histo = LogEtrue_projX->GetRMS(1);
//				cout << "LogErec_histo: " << LogErec_histo << "  Sigma_LogErec_histo: " << Sigma_LogErec_histo << endl;
				Double_t bias = LogErec_histo/LogEtrue;
				Double_t res= Sigma_LogErec_histo;
//				cout << "Bias: " << bias << "  Resolution: " << res << endl;
				BiasvslogEtrue->SetBinContent(i,bias);
				ResvslogEtrue->SetBinContent(i, res);
				if(i>=binY_Emin && i<=binY_Emax){
					BiasvslogEtrue_cut->SetBinContent(i-binY_Emin+1, bias);
					ResvslogEtrue_cut->SetBinContent(i-binY_Emin+1, res);
				}
			}
		}
	}

		//---Plot resulting graphs---
		TCanvas* ResandBiasLogEtrue_canvas = new TCanvas();
		ResandBiasLogEtrue_canvas->SetTitle("Resolution and Bias vs LogEtrue");
		ResandBiasLogEtrue_canvas->Divide(2,1);
		ResandBiasLogEtrue_canvas->cd(1);
//		BiasvslogEtrue_cut->Smooth(smooth_times);
		BiasvslogEtrue_cut->GetYaxis()->SetTitle("LogErec/LogEtrue");
		BiasvslogEtrue_cut->GetXaxis()->SetTitle("LogEtrue");
		BiasvslogEtrue_cut->SetLineColor(kBlue);
		BiasvslogEtrue_cut->SetLineWidth(2);
		BiasvslogEtrue_cut->DrawCopy("");
		ResandBiasLogEtrue_canvas->cd(2);
//		ResvslogEtrue_cut->Smooth(smooth_times);
		ResvslogEtrue_cut->GetYaxis()->SetTitle("#sigma");
		ResvslogEtrue_cut->GetXaxis()->SetTitle("LogEtrue");
		ResvslogEtrue_cut->SetLineWidth(2);
		ResvslogEtrue_cut->SetLineColor(kRed);
		ResvslogEtrue_cut->DrawCopy("");

		if(draw_fit == true){
			TCanvas* fit_canvas = new TCanvas();
			FitvsLogEtrue->SetStats(1);
			FitvsLogEtrue->DrawCopy();
			FitGausLogEtrue->Draw("same");
		}

		//Give names and save everything in the root file
		BiasvslogEtrue->SetName("BiasVsLogEtrue");
		ResvslogEtrue->SetName("ResolutionVsLogEtrue");
		BiasvslogEtrue->Write();
		ResvslogEtrue->Write();
		BiasvslogEtrue_cut->SetName("BiasvslogEtrue_cut");
		ResvslogEtrue_cut->SetName("ResvslogEtrue_cut");
		BiasvslogEtrue_cut->Write();
		ResvslogEtrue_cut->Write();



		//-----Plots vs LogErec-----
		TH1D* LogErec_projY;
		TH1D* LogErec_projY_corrected = new TH1D("Proj_corr", "Proj_corr" , numMigBinsY,TMath::Log10(arrayMigBinsY[0]),TMath::Log10(arrayMigBinsY[numMigBinsY]));
		TF1* LogFit_projY;
		TF1* LogFit_projY_corrected;
		Double_t LogErec;

		TH1D* BiasvslogErec = new TH1D("BiasVsLogErec", "BiasVsLogErec", numMigBinsX, TMath::Log10(arrayMigBinsX[0]),TMath::Log10(arrayMigBinsX[numMigBinsX]));
		TH1D* ResvslogErec = new TH1D("ResolutionVsLogErec", "ResolutionVsLogErec", numMigBinsX, TMath::Log10(arrayMigBinsX[0]),TMath::Log10(arrayMigBinsX[numMigBinsX]));
		TH1D* BiasvslogErec_cut = new TH1D("BiasvslogErec_cut", "BiasvslogErec_cut", numMigBinsX_cut, TMath::Log10(arrayMigBinsX_cut[0]),TMath::Log10(arrayMigBinsX_cut[numMigBinsX_cut]));
		TH1D* ResvslogErec_cut = new TH1D("ResvslogErec_cut", "ResvslogErec_cut", numMigBinsX_cut, TMath::Log10(arrayMigBinsX_cut[0]),TMath::Log10(arrayMigBinsX_cut[numMigBinsX_cut]));

		TH1D* FitvsLogErec;
		TF1* FitGausLogErec;
		TH1D* FitvsLogErec_corr;
		TF1* FitGausLogErec_corr;

		for (Int_t i = 1; i <=numMigBinsX; i++) {
			LogErec = LogMigMatrix->GetXaxis()->GetBinCenter(i);
//			cout << "LogErec: " << i << " " << LogErec << endl;
			LogErec_projY = LogMigMatrix->ProjectionY("Proj in LogEtrue",i,i, "");
			for(Int_t j = 1; j<=numMigBinsY; j++){
				if(i == draw_fit_number){
//					cout << "Bin number: " << i << endl;
//					cout << "Old content: " << LogErec_projY->GetBinContent(j) << endl;
//					cout << "Factor: " << TMath::Power(10., 1.6*LogErec_projY->GetBinCenter(j)) << endl;
//					cout << "New content: " << LogErec_projY->GetBinContent(j)*TMath::Power(10., 0.6*LogErec_projY->GetBinCenter(j)) << endl;
				}

				LogErec_projY_corrected->SetBinContent(j, LogErec_projY->GetBinContent(j)*TMath::Power(10., 0.6*LogErec_projY->GetBinCenter(j)));
			}
			if(LogErec_projY->GetMean(1)>0){ //This assures bins are not empty
				if(fit_gaus == true){
					LogErec_projY->Fit("gaus","Q0");
					LogFit_projY = LogErec_projY->GetFunction("gaus");
					LogErec_projY_corrected->Fit("gaus","Q0");
					LogFit_projY_corrected = LogErec_projY_corrected->GetFunction("gaus");
					if(draw_fit == true && i == draw_fit_number){
						FitvsLogErec = (TH1D*)LogErec_projY->Clone("FitvsLogErec");
						FitGausLogErec = (TF1*)LogFit_projY->Clone("FitGausLogErec");
						FitvsLogErec_corr = (TH1D*)LogErec_projY_corrected->Clone("FitvsLogErec_corr");
						FitGausLogErec_corr = (TF1*)LogFit_projY_corrected->Clone("FitGausLogErec_corr");
					}
					Double_t LogEtrue_fit = LogFit_projY_corrected->GetParameter(1);
					Double_t Sigma_LogEtrue_fit = LogFit_projY_corrected->GetParameter(2);
//					cout << "LogEtrue_fit: " << LogEtrue_fit << "  Sigma_LogEtrue_fit: " << Sigma_LogEtrue_fit << endl;
					Double_t bias = LogEtrue_fit/LogErec;
					Double_t res= Sigma_LogEtrue_fit;
//					cout << "Bias: " << bias << "  Resolution: " << res << endl;
					BiasvslogErec->SetBinContent(i, bias);
					ResvslogErec->SetBinContent(i, res);
					if(i>=binX_Emin && i<=binX_Emax){
						BiasvslogErec_cut->SetBinContent(i-binX_Emin+1, bias);
						ResvslogErec_cut->SetBinContent(i-binX_Emin+1, res);
					}

				}
				else{
					Double_t LogEtrue_histo = LogErec_projY_corrected->GetMean(1);
					Double_t Sigma_LogEtrue_histo = LogErec_projY_corrected->GetRMS(1);
//					cout << "LogEtrue_histo: " << LogEtrue_histo << "  Sigma_LogEtrue_histo: " << Sigma_LogEtrue_histo << endl;
					Double_t bias = LogEtrue_histo/LogErec;
					Double_t res= Sigma_LogEtrue_histo;
//					cout << "Bias: " << bias << "  Resolution: " << res << endl;
					BiasvslogErec->SetBinContent(i, bias);
					ResvslogErec->SetBinContent(i, res);
					if(i>=binX_Emin && i<=binX_Emax){
						BiasvslogErec_cut->SetBinContent(i-binX_Emin+1, bias);
						ResvslogErec_cut->SetBinContent(i-binX_Emin+1, res);
					}
				}
			}
		}



		//---Plot resulting graphs---
		TCanvas* ResandBiasLogErec_canvas = new TCanvas();
		ResandBiasLogErec_canvas->SetTitle("Resolution and Bias vs LogErec");
		ResandBiasLogErec_canvas->Divide(2,1);
		ResandBiasLogErec_canvas->cd(1);
//		BiasvslogErec_cut->Smooth(smooth_times);
		BiasvslogErec_cut->GetYaxis()->SetTitle("LogEtrue/LogErec");
		BiasvslogErec_cut->GetXaxis()->SetTitle("LogErec");
		BiasvslogErec_cut->SetLineColor(kBlue);
		BiasvslogErec_cut->SetLineWidth(2);
		BiasvslogErec_cut->DrawCopy("");
		ResandBiasLogErec_canvas->cd(2);
//		ResvslogErec_cut->Smooth(smooth_times);
		ResvslogErec_cut->GetYaxis()->SetTitle("#sigma");
		ResvslogErec_cut->GetXaxis()->SetTitle("LogErec");
		ResvslogErec_cut->SetLineWidth(2);
		ResvslogErec_cut->SetLineColor(kRed);
		ResvslogErec_cut->DrawCopy("");


		if(draw_fit == true){
			TCanvas* fit_canvas = new TCanvas();
			fit_canvas->Divide(2,1);
			fit_canvas->cd(1);
			FitvsLogErec->SetStats(1);
			FitvsLogErec->DrawCopy();
			FitGausLogErec->DrawCopy("same");
			fit_canvas->cd(2);
			FitvsLogErec_corr->SetStats(1);
			FitvsLogErec_corr->DrawCopy("EB");
			FitGausLogErec_corr->DrawCopy("same");


//			TCanvas* normhistos = new TCanvas();
//			FitvsLogErec->Scale(1./FitvsLogErec->Integral());
//			FitvsLogErec->SetLineColor(kBlue);
//			FitvsLogErec->Draw("HIST");
//			FitvsLogErec_corr->Scale(1./FitvsLogErec_corr->Integral());
//			FitvsLogErec_corr->SetLineColor(kRed);
//			FitvsLogErec_corr->Draw("same");
		}

		//Give names to save everything in the root file
		BiasvslogErec->SetName("BiasVsLogErec");
		ResvslogErec->SetName("ResolutionVsLogErec");
		BiasvslogErec->Write();
		ResvslogErec->Write();
		BiasvslogErec_cut->SetName("BiasvslogErec_cut");
		ResvslogErec_cut->SetName("ResvslogErec_cut");
		BiasvslogErec_cut->Write();
		ResvslogErec_cut->Write();
		IRFoutputFile->Close();




}



