/*
 * Peak_Search.C
 *
 *  Created on: Jan 10, 2018
 *      Author: lnogues
 *
 *
 *      Code to look for peaks in the time histogram of the events before trying to fit it
 *      with LC_Like fit.
 *
 *
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
#include <TMath.h>
#include <iostream>
#include <TSpectrum.h>
#include "TTree.h"

using namespace std;

void Peak_Search(){

	//Read the histogram of the lightcurve
	TFile* HistoFile = TFile::Open("./FinalEventLists/LC_01-20TeV_ONhistograms_Mrk421_2014.root", "READ");
			if (!HistoFile) {
				cout << "ERROR: file " << HistoFile << " not found..."	<< endl;
			}
		TH1D* lightcurve = dynamic_cast<TH1D*>(HistoFile->Get("LC_insecs"));
		if (!lightcurve) {
			cout << "LC_insecs not found" << endl;
		}

	Int_t numberOfBins = lightcurve->GetNbinsX();
	const Double_t* arrayBins = lightcurve->GetXaxis()->GetXbins()->GetArray();
	cout << "Number of bins in LC: " << numberOfBins << endl;

//	TCanvas* LC_canvas = new TCanvas();
//	lightcurve->Draw();

	//Interpolate the histogram
	std::vector<double> x_points;
	std::vector<double> y_points;

	//----Set first point from histogram limit
	x_points.push_back(0);
	y_points.push_back(lightcurve->GetBinContent(1));


	for(Int_t bin = 1; bin < numberOfBins; bin++){
		if(lightcurve->GetBinContent(bin)!= 0){
			x_points.push_back(lightcurve->GetBinCenter(bin));
			y_points.push_back(lightcurve->GetBinContent(bin));
		}
	}

	//----Set last point from histogram limit
	x_points.push_back(lightcurve->GetBinLowEdge(numberOfBins+1));
	y_points.push_back(lightcurve->GetBinContent(numberOfBins));

	ROOT::Math::Interpolator int_aspline(x_points, y_points, ROOT::Math::Interpolation::kAKIMA);

	//----Do interpolated histogram
	TH1D* aspline_histo = new TH1D("TH1_aspline", "TH1_aspline", numberOfBins, arrayBins);
		for(Int_t bin = 1; bin < numberOfBins; bin++){
			aspline_histo->SetBinContent(bin, int_aspline.Eval(lightcurve->GetBinCenter(bin)));
		}

	TCanvas* histo_canvas = new TCanvas();
	aspline_histo->Draw();

	//Crear el TSpectrum
	TSpectrum* spec_lc = new TSpectrum();
	Int_t numberOfPeaks = spec_lc->Search(aspline_histo,5, "", 0.2);
	cout << "Number of peaks: " << numberOfPeaks << endl;

	//Get Results
	Int_t peaknumber = spec_lc->GetNPeaks();
	Float_t* x_pos = spec_lc->GetPositionX();
	Float_t* y_pos = spec_lc->GetPositionY();
	cout << "We found " << peaknumber << " peaks at locations: " << endl;
	for(Int_t i=0; i<peaknumber; i++){
		cout << "x pos: " << x_pos[i] << " y pos: " << y_pos[i] << endl;
	}


	//Draw polymarker
	TList *functions = aspline_histo->GetListOfFunctions();
	TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");


























}


