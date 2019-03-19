/*
 * Histo_Interpolate.C
 *
 *  Created on: Oct 16, 2017
 *      Author: Leyre Nogués
 *
 *      Code to play with the time histogram of the events:
 *      Interpolation, fit and delay.
 *
 *      The codes also fills with fake events the time gaps to apply the Kernel Density method.
 */


#include <Math/InterpolationTypes.h>
#include <Math/Interpolator.h>
#include <Rtypes.h>
#include <TArrayD.h>
#include <TAttLine.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TVector.h>
#include <iostream>
#include "TTree.h"
#include <TRandom3.h>
#include<TStyle.h>


using namespace std;

Double_t LCGauss(Double_t* x, Double_t* par){
//	par[0] = 31.54;
//	par[1] = 84.50;
//	par[2] = 4294.06;
//	par[3] = 2350.84;
//	par[4] = 1561.74;
	Double_t gaus_value = 0;

	//Compute values of the gaussian
	if(x[0]<par[2])
		gaus_value = par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[3],2));
	else
		gaus_value = par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[4],2));

	return gaus_value;
}

static Double_t One_Asymgaus(Double_t* x, Double_t*  par){

	if(x[0]<par[2]) return par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[3],2));
	else return par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[4],2));
}

void Histo_Interpolate(){

	//Read Events and histogram.
	//--Time histogram--
	TFile* HistoFile = TFile::Open("./Paper_Histos/Histo_36ON_Paper.root", "READ");

	if (!HistoFile) {
			cout << "ERROR: file " << HistoFile << " not found..."	<< endl;
		}
	TH1D* lightcurve = dynamic_cast<TH1D*>(HistoFile->Get("ONvariation"));
	if (!lightcurve) {
		cout << "LC_insecs not found" << endl;
	}

	lightcurve->SetStats(1);
	Int_t numberOfBins = lightcurve->GetNbinsX();
	const Double_t* arrayBins = lightcurve->GetXaxis()->GetXbins()->GetArray();
	cout << "Number of bins in LC: " << numberOfBins << endl;

	TCanvas* plothistocanvas =  new TCanvas();
	lightcurve->SetTitle("Akima Interpolation");
	lightcurve ->Draw();

	//Interpolation of the histogram
	std::vector<double> x_points;
	std::vector<double> y_points;
	std::vector<double> bin_width_points;
	std::vector<double> x_holes;

	// Create a TTree
	TTree t1("t1","Tree with simple variables");
	Double_t x_point, y_point, bin_width;
	t1.Branch("x",&x_point,"px/D");
	t1.Branch("y",&y_point,"py/D");
	t1.Branch("width",&bin_width,"pz/D");


	Int_t bins_index = 0, holes_index = 0;

	//Set first point from histogram limit
	x_points.push_back(0);
	y_points.push_back(lightcurve->GetBinContent(1));
	bin_width_points.push_back(lightcurve->GetBinWidth(1));

	bins_index++;
	x_point = 0;
	y_point = lightcurve->GetBinContent(1);
	bin_width = lightcurve->GetBinWidth(1);
	t1.Fill();


	for(Int_t bin = 1; bin <= numberOfBins; bin++){
		if(lightcurve->GetBinContent(bin)!= 0){
			x_points.push_back(lightcurve->GetBinCenter(bin));
			x_point = lightcurve->GetBinCenter(bin);
			y_points.push_back(lightcurve->GetBinContent(bin));
			y_point = lightcurve->GetBinContent(bin);
			bin_width_points.push_back(lightcurve->GetBinWidth(bin));
			bin_width = lightcurve->GetBinWidth(bin);
			bins_index++;
			t1.Fill();
		}
		else{
			x_holes.push_back(lightcurve->GetBinCenter(bin));
			holes_index++;
		}
	}

	//Set last point from histogram limit
	x_points.push_back(lightcurve->GetBinLowEdge(numberOfBins+1));
	y_points.push_back(lightcurve->GetBinContent(numberOfBins));
	bin_width_points.push_back(lightcurve->GetBinWidth(numberOfBins));
	cout << "Last point x: " << lightcurve->GetBinLowEdge(numberOfBins+1) << endl;
	cout << "Last point y: " << lightcurve->GetBinContent(numberOfBins) << endl;
	cout << "Last point width: " << lightcurve->GetBinWidth(numberOfBins) << endl;
	x_point = lightcurve->GetBinLowEdge(numberOfBins+1);
	y_point = lightcurve->GetBinContent(numberOfBins);
	bin_width = lightcurve->GetBinWidth(numberOfBins);
	t1.Fill();
	bins_index++;


	cout << "High Edge: " << lightcurve->GetBinLowEdge(numberOfBins+1) << endl;

	cout << "Full bins: " << bins_index << endl;
	cout << "Empty bin: " << holes_index << endl;

	cout << "First width: " << bin_width_points[0] << endl;
	cout << "Second width: " << bin_width_points[1] << endl;
	cout << "Last width: " << bin_width_points[bin_width_points.size()-1] << endl;
	cout << "Pre-last width: " << bin_width_points[bin_width_points.size()-2] << endl;





	//---------------------------Linear interpolation---------------------------------------
	ROOT::Math::Interpolator int_linear(x_points, y_points, ROOT::Math::Interpolation::kLINEAR);
	TH1D* linear_histo = new TH1D("TH1_linear", "TH1_linear", numberOfBins, arrayBins);
	TGraph* linear_graph = new TGraph(numberOfBins-1);
	for(Int_t bin = 1; bin < numberOfBins; bin++){
		linear_histo->SetBinContent(bin, int_linear.Eval(lightcurve->GetBinCenter(bin)));
		linear_graph->SetPoint(bin-1, lightcurve->GetBinCenter(bin), int_linear.Eval(lightcurve->GetBinCenter(bin)));
	}

//	linear_histo->SetLineColor(kRed);
//	linear_histo->Draw("same");
//	linear_graph->SetLineColor(kBlack);
//	linear_graph->Draw("csame");





	//--------------------------------Cubic Spline interpolation---------------------------------
	ROOT::Math::Interpolator int_cspline(x_points, y_points, ROOT::Math::Interpolation::kCSPLINE);
	TH1D* cspline_histo = new TH1D("TH1_cspline", "TH1_cspline", numberOfBins, arrayBins);
	TGraph* cspline_graph = new TGraph(numberOfBins-1);
	for(Int_t bin = 1; bin < numberOfBins; bin++){
		cspline_histo->SetBinContent(bin, int_cspline.Eval(lightcurve->GetBinCenter(bin)));
		cspline_graph->SetPoint(bin-1, lightcurve->GetBinCenter(bin), int_cspline.Eval(lightcurve->GetBinCenter(bin)));
	}

//	cspline_histo->SetLineColor(kRed);
//	cspline_histo->Draw("same");
//	cspline_graph->SetLineColor(kBlack);
//	cspline_graph->Draw("csame");





	//-------------------Akima Spline interpolation------------------------- (Stable to outliers)
	ROOT::Math::Interpolator int_aspline(x_points, y_points, ROOT::Math::Interpolation::kAKIMA);
	TH1D* aspline_histo = new TH1D("TH1_aspline", "TH1_aspline", numberOfBins, arrayBins);
	TGraph* aspline_graph = new TGraph(numberOfBins-1);
	for(Int_t bin = 1; bin <=numberOfBins; bin++){
		aspline_histo->SetBinContent(bin, int_aspline.Eval(lightcurve->GetBinCenter(bin)));
		aspline_graph->SetPoint(bin-1, lightcurve->GetBinCenter(bin), int_aspline.Eval(lightcurve->GetBinCenter(bin)));
	}


//	TCanvas* hello = new TCanvas();
//	aspline_histo->SetLineColor(kRed);
//	aspline_histo->Draw("same");
//	lightcurve->Draw("same");
//	aspline_graph->SetLineColor(kBlack);
//	aspline_graph->Draw("csame");







	//------------------------------------For Kernel Density method------------------------------

	//Create list with extra events in wobble gaps
//	std::vector<double> extra_time_events;
//	TRandom3* rant = new TRandom3();
//	rant->SetSeed(0);
//	for(Int_t bin = 1; bin <= numberOfBins; bin++){
//		if(lightcurve->GetBinContent(bin)==0){
//			Int_t numberGenEvents = round(aspline_histo->GetBinContent(bin)*lightcurve->GetBinWidth(bin)/lightcurve->GetBinWidth(bin+1));
//			for(Int_t i= 0; i<numberGenEvents; i++){
//				Double_t time = rant->Uniform(aspline_histo->GetBinLowEdge(bin)+1, aspline_histo->GetBinLowEdge(bin+1)-1);
//				extra_time_events.push_back(time);
//				lightcurve->Fill(time);
//			}
//			cout << "After bin: " << bin << " we had " << extra_time_events.size() << " extra events" << endl;
//		}
//	}
//
//	lightcurve->Draw("same");
//
//	FILE *fout;
//	fout = fopen("./ExtraEventsTime_OFF_Kazuma_120.txt","wb");
//	Double_t energy = 0.;
//	for(Int_t i = 0; i<extra_time_events.size(); i++)
//	{
//		fprintf(fout,"%0.8lf    %0.8lf   \n", extra_time_events[i], energy);
//	}



	//Create also extra events out of the window limits.
//	std::vector<double> extra_time_events;
//	Double_t extra_time = 8000.; //s
//	Int_t extra_events = 20000;
//	TF1* lc_function = new TF1("lc_function", LCGauss, arrayBins[0]-extra_time, arrayBins[numberOfBins-1]+extra_time,5);
//	lc_function->SetParameters(29.93, 72.86, 4407.01, 2307.10, 1375.21);
//	TRandom3* ranx = new TRandom3();
//	ranx->SetSeed(0);
//	TRandom3* rany = new TRandom3();
//	rany->SetSeed(1);
//	Int_t number_of_bins = (int)((arrayBins[numberOfBins-1]+extra_time) - (arrayBins[0]-extra_time))/60; //60s bins
//	cout << "Bins: " << number_of_bins << endl;
//	TH1D* outOfLimits_histo = new TH1D("hi","hi", number_of_bins, arrayBins[0]-extra_time, arrayBins[numberOfBins-1]+extra_time);
//	Int_t simulated_events = 0;
//	Int_t accepted_events = 0;
//
//
//
//	for(Int_t i = 0; i<extra_events; i++){
//		Double_t x = ranx->Uniform(arrayBins[0]-extra_time, arrayBins[numberOfBins-1]+extra_time);
//		Double_t y = rany->Uniform(0, 120);
//		simulated_events++;
//		if(x > arrayBins[0] && x<arrayBins[numberOfBins-1]) continue;
//		if(y < lc_function->Eval(x)){
//			accepted_events++;
//			extra_time_events.push_back(x);
//			outOfLimits_histo->Fill(x);
//		}
//	}
//
//	lc_function->Draw("");
//	outOfLimits_histo->Draw("same");
//	cout << "Sim events: " <<  simulated_events << endl;
//	cout << "Acep events: " <<  accepted_events << endl;
//	cout << "Array size: " <<  extra_time_events.size() << endl;
//
//
//	FILE *fout;
//	fout = fopen("./ExtraLimitsEventsTime_ON_Kazuma_1min.txt","wb");
//	Double_t energy = 0.;
//	for(Int_t i = 0; i<extra_time_events.size(); i++)
//	{
//		fprintf(fout,"%0.8lf    %0.8lf   \n", extra_time_events[i], energy);
//	}



	 //Extra events out of the window limits for the background
//		std::vector<double> extra_time_events;
//		Double_t extra_time = 8000.; //s
//		Int_t extra_events = 1800;
//		TRandom3* ranx = new TRandom3();
//		ranx->SetSeed(0);
//		TRandom3* rany = new TRandom3();
//		rany->SetSeed(1);
//		Int_t number_of_bins = (int)((arrayBins[numberOfBins-1]+extra_time) - (arrayBins[0]-extra_time))/60; //60s bins
//		cout << "Bins: " << number_of_bins << endl;
//		TH1D* outOfLimits_histo = new TH1D("hi","hi", number_of_bins, arrayBins[0]-extra_time, arrayBins[numberOfBins-1]+extra_time);
//
//
//		for(Int_t i = 0; i<extra_events; i++){
//			Double_t x = ranx->Uniform(arrayBins[0]-extra_time, arrayBins[0]);
//			Double_t y = rany->Uniform(arrayBins[numberOfBins-1], arrayBins[numberOfBins-1]+extra_time);
//			extra_time_events.push_back(x);
//			extra_time_events.push_back(y);
//			outOfLimits_histo->Fill(x);
//			outOfLimits_histo->Fill(y);
//
//		}


//		outOfLimits_histo->Draw("");
//		lightcurve->Draw("same");
//		cout << "Array size: " <<  extra_time_events.size() << endl;
//
//
//		FILE *fout;
//		fout = fopen("./ExtraLimitsEventsTime_OFF_Kazuma_1min.txt","wb");
//		Double_t energy = 0.;
//		for(Int_t i = 0; i<extra_time_events.size(); i++)
//		{
//			fprintf(fout,"%0.8lf    %0.8lf   \n", extra_time_events[i], energy);
//		}




//------------------------------------------------------------------------------------------






	//Interpolator offer integration and derivation methods. TGraphs are only used to "see" them.
	Int_t number_of_points = 10000;
	TGraph* general =  new TGraph(number_of_points);
	Double_t step = 12500/number_of_points;
	Double_t point = x_points[0];
	Int_t iter=0;

	cout << "First point: " << x_points[0] << " Last point: " << x_points[bins_index-1] << endl;

	while(point+iter*step<=x_points[bins_index-1]){
//		cout << "Step: " << point+iter*step << endl;
		general->SetPoint(iter, point+iter*step, int_aspline.Eval(point+iter*step));
		iter++;
	}

//	general->SetTitle("Delayed LC, alpha=60, E=20 TeV");
//	general->SetLineColor(kBlack);
	general->SetLineWidth(2);
	general->Draw("lsame");




//------------------------------Test to delay in time the interpolation------------------------


//	//Moving the TGraph
//
//	//--Amount of delay--
//	Double_t alpha = 60;
//	Double_t z=0.031;
//	Double_t H0=70.;
//	Double_t PlanckScale=1.2*TMath::Power(10,19);
//	Double_t QGpar=alpha;
//	Double_t PcTom=3.085*TMath::Power(10,19);
//	Double_t QGfixedPart=z*PcTom/H0;
//	std::cout << "alpha: " << alpha << std::endl;
//
//	Double_t QGDelay= QGfixedPart*QGpar/PlanckScale;
//	std::cout << "Injected tau: " << QGDelay << std::endl;
//	Double_t energy = 20000;
//
//	Double_t final_delay = QGDelay*energy; //Energy in GeV
//	std::cout << "Final Delay: " << QGDelay << std::endl;
//
//
//
//
//	//--Delayed TGraph--
//
//	TGraph* general_delayed =  new TGraph(number_of_points);
//	iter=0;
//
//	while(point+iter*step<=x_points[207]){
//		if((point+iter*step)-final_delay > x_points[0])general_delayed->SetPoint(iter, point+iter*step, int_aspline.Eval((point+iter*step)-final_delay));
//		else general_delayed->SetPoint(iter, point+iter*step, int_aspline.Eval(x_points[0]));
//		iter++;
//	}
//	general_delayed->SetLineColor(kRed);
//	general_delayed->SetLineWidth(2);
//	general_delayed->Draw("same");



//-------------------------------------------------------------------------------------------




//Create a root file to save the interpolation
//We save a tree with the points and the widths (Needed to apply later Poissonian fluctuations)

	 TFile* inter_outputfile = new TFile("Inter_Paper_36ONevents_width.root","recreate");
//	 x_vector.Write();
//	 y_vector.Write();
//	 error_vector.Write();
	 t1.Write();
//	 general->Write();




//--------------Extrapolation out of the observation window with an asymmetric Gauss----------
//Only for the ON interpolation!!

	gStyle->SetOptFit(111111);
	TF1* function = new TF1("function", One_Asymgaus, 0, 13000, 5);
	function->SetNpx(1000);
	function->SetParameters(0.5, 0.35, 4000, 2590, 1500);
	lightcurve->Fit(function);
	function->Draw("same");
	plothistocanvas->SetName("Interpolation_canvas");
	plothistocanvas->Write();
	inter_outputfile->Close();



}

