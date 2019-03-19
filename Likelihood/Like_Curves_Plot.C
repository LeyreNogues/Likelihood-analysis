/*
 * Like_Curves_Plot.C
 *
 *  Created on: Sep 14, 2018
 *      Author: lnogues
 */

#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include<iostream>
#include<fstream>
#include<istream>

void Like_Curves_Plot(){

	TFile* HistoFile = TFile::Open("./Likelihood_Kazuma_120_36Inter_Quadratic_DATA2.root", "READ");
	TGraph *curve = (TGraph*)HistoFile->Get("");

	TCanvas* canvas = new TCanvas();
	curve->GetXaxis()->SetTitle("#alpha");
	curve->GetYaxis()->SetTitle("-2#Delta L");
	curve->SetLineColor(kRed);
	curve->SetLineWidth(2);
	curve->Draw();


}



