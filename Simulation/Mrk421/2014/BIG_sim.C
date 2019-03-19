/*
 * BIG_sim.C
 *
 *  Created on: Jan 20, 2017
 *      Author: Leyre Nogués
 *
 *      Modified version to simulate using a Interpolation or a function. Oct2017.
 *      Modified version to generate a big number of events. Febrary 2018.
 *      Modified version to use a KDE function for event simulation. May 2018.
 *      Modified version to include EBL, two Gausians. June 2018.
 *      Modified version to include Background events. August 2018.
 *      Modified version to include or not the Background with a Boolean. August 2018.
 *
 */

#include <Math/InterpolationTypes.h>
#include <Math/Interpolator.h>
#include <Rtypes.h>
#include <stddef.h>
#include <TArrayD.h>
#include <TAttLine.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>


using namespace std;

//Functions for the PDF building
TGraph* CollAr_Graph;
TH1D *energy_dist_bkg;
TF1* flare_spectrum;
TF1* folded_spectrum_flare;
Double_t inter_min = 0.;
Double_t inter_max = 0.;
Double_t inter_min_bkg = 0.;
Double_t inter_max_bkg = 0.;
Double_t KDE_max = 0.;
Double_t KDE_max_bkg = 0.;
ROOT::Math::Interpolator *int_aspline = NULL;
ROOT::Math::Interpolator *int_aspline_bkg = NULL;
TF1* KDE_function;
TF1* KDE_function_bkg;
TF1* KDE_function_scaled;
TF1* energy_template;
TF1* energy_template_scaled;
TGraph* general;
TGraph* general_bkg;
TH1D* sim_flare_time;
TH1D* sim_flare_time_var;
TH1D* sim_bkg_time;
TH1D* sim_bkg_time_var;
TH1D* timedel_bkg;
TH1D* timedel_bkg_var;


//Graph for EBL application
TGraph *transvse;

//Spectrum variables
const Double_t ln10 = 2.302585093;
const Double_t E0 = 364.;
Double_t bin_value;
Double_t factor_scale;

//----------------------------Wobble information---------------------------------------

static const  Int_t numberOfWobbles = 14; //Number of wobbles
static const  Int_t numberOfPoints = 2*numberOfWobbles;

static const  Double_t mintimeMJD = 1.77293506944444380e+03; //From wobble info
static const Double_t maxtimeMJD = 1.77308576388889196e+03; //From wobble info

const static Double_t hole[numberOfPoints] = { 0.,  //Beggining of first wobble
 (1.77294540509259241e+03 - mintimeMJD)*24.*60.*60.,
 (1.77294579861110833e+03 - mintimeMJD)*24.*60.*60., //Wobble2
(1.77295609953703388e+03 - mintimeMJD)*24.*60.*60.,
(1.77295817129629722e+03 - mintimeMJD)*24.*60.*60., //Wobble3
(1.77296849537036906e+03 - mintimeMJD)*24.*60.*60.,
(1.77296908564814657e+03 - mintimeMJD)*24.*60.*60., //Wobble4
(1.77297917824074102e+03 - mintimeMJD)*24.*60.*60.,
(1.77297956018518744e+03 - mintimeMJD)*24.*60.*60., //Wobble 5
(1.77298983796295943e+03 - mintimeMJD)*24.*60.*60.,
(1.77299023148148262e+03 - mintimeMJD)*24.*60.*60., //Wobble6
(1.77300049768518511e+03 - mintimeMJD)*24.*60.*60.,
(1.77300089120370103e+03 - mintimeMJD)*24.*60.*60., //Wobble7
(1.77301114583333401e+03 - mintimeMJD)*24.*60.*60.,
(1.77301156250000349e+03 - mintimeMJD)*24.*60.*60., //Wobble8
(1.77302184027777548e+03 - mintimeMJD)*24.*60.*60.,
(1.77302224537036818e+03 - mintimeMJD)*24.*60.*60., //Wobble9
(1.77303248842592438e+03 - mintimeMJD)*24.*60.*60.,
(1.77303289351851708e+03 - mintimeMJD)*24.*60.*60., //Wobble10
(1.77304315972221957e+03 - mintimeMJD)*24.*60.*60.,
(1.77304371527778130e+03 - mintimeMJD)*24.*60.*60., //Wobble11
(1.77305380787036847e+03 - mintimeMJD)*24.*60.*60.,
(1.77305427083333052e+03 - mintimeMJD)*24.*60.*60., //Wobble12
(1.77306445601851738e+03 - mintimeMJD)*24.*60.*60.,
(1.77306501157407183e+03 - mintimeMJD)*24.*60.*60., //Wobble 13
(1.77307512731481256e+03 - mintimeMJD)*24.*60.*60.,
(1.77307554398148204e+03 - mintimeMJD)*24.*60.*60., //Wobble14
(maxtimeMJD - mintimeMJD)*24.*60.*60.}; //End of last wobble (not last event)


Double_t FoldingFlare(Double_t* x, Double_t* par){
	Double_t flareValue = flare_spectrum->Eval(x[0])*10;
	Double_t areaValue = CollAr_Graph->Eval(x[0]);
	if(areaValue>0)return flareValue*areaValue;
	else return 0;
}

Double_t EBLFlare(Double_t* x, Double_t* par){
	Double_t flareValue = folded_spectrum_flare->Eval(x[0]);
	Double_t transValue = transvse->Eval(log10(x[0]));
//	cout << "x[0]: " << x[0] << endl;
//	cout << "CollAr spectrum: " << flareValue << endl;
//	cout << "Trans value: " << transValue << endl;
	if(flareValue>0 && transValue>0) return flareValue*transValue;
	else return 0;
}

Double_t ScaleKDE(Double_t* x, Double_t* par){
	Double_t function_value = KDE_function->Eval(x[0]);
	Double_t factor = sim_flare_time->GetBinContent(sim_flare_time->FindBin(4000.))/KDE_function->Eval(4000.);
	return function_value*factor;
}

Double_t ScaleEnergyFunction(Double_t* x, Double_t* par){
	Double_t originalvalue = energy_template->Eval(x[0]); //a
	return originalvalue*factor_scale;
}

Double_t lc_function(Double_t t){ //CAREFUL! THIS HAS TO BE UPDATED WITH THE FIT VALUES FROM Histo_Interpolate.C

	Double_t cons_value = 0.5365;
	Double_t weight_value = 1.354;
	Double_t peak_value = 4358.44;
	Double_t sigma_left = 2342.71;
	Double_t sigma_right = 1481.20;
	Double_t gaus_value = 0;


	//Compute values of the gaussian
	if(t<peak_value)
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_left,2));
	else
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_right,2));


	//Return Interpolation or Gaussian Extrapolation
    if(t<inter_min || t>inter_max)
    	return gaus_value;
    else
		return int_aspline->Eval(t);
}

Double_t lc_function_background(Double_t t){
	return int_aspline_bkg->Eval(t);
}

Double_t gaus_function(Double_t t){

	Double_t cons_value = 0.5365;
	Double_t weight_value = 1.354;
	Double_t peak_value = 4358.44;
	Double_t sigma_left = 2342.71;
	Double_t sigma_right = 1481.20;
	Double_t gaus_value = 0;

	//Compute value of the gaussian
	if(t<peak_value)
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_left,2));
	else
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_right,2));

	return gaus_value;
}


Double_t GetEnergyFromPowerLaw(const Double_t emin, const Double_t emax, const Double_t index)
{
  const Double_t x = gRandom->Rndm();
  const Double_t np1 = index;
  const Double_t x0 = TMath::Power(emin,np1);
  const Double_t x1 = TMath::Power(emax,np1);
  return TMath::Power((x1 - x0) * x + x0,1./np1);
}

Double_t get_E(Double_t y)
{
  return exp(y)*E0;
}


Double_t Double_Gaus(Double_t* x, Double_t*  par){

	Double_t simgaus = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	Double_t simgaus2 = par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5],2));
	return simgaus + simgaus2;
}

Double_t Simple_Gaus(Double_t* x, Double_t*  par){

	Double_t simgaus = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	return simgaus;
}

Double_t evaluate_wobble_condition(Double_t t){
	Double_t wobble_condition = 0.;
	if(t>hole[1] && t<hole[2]) wobble_condition = 1.;
	if(t>hole[3] && t<hole[4]) wobble_condition = 1.;
	if(t>hole[5] && t<hole[6]) wobble_condition = 1.;
	if(t>hole[7] && t<hole[8]) wobble_condition = 1.;
	if(t>hole[9] && t<hole[10]) wobble_condition = 1.;
	if(t>hole[11] && t<hole[12]) wobble_condition = 1.;
	if(t>hole[13] && t<hole[14]) wobble_condition = 1.;
	if(t>hole[15] && t<hole[16]) wobble_condition = 1.;
	if(t>hole[17] && t<hole[18]) wobble_condition = 1.;
	if(t>hole[19] && t<hole[20]) wobble_condition = 1.;
	if(t>hole[21] && t<hole[22]) wobble_condition = 1.;
	if(t>hole[23] && t<hole[24]) wobble_condition = 1.;
	if(t>hole[25] && t<hole[26]) wobble_condition = 1.;

	return wobble_condition;
}

void BIG_sim(double alpha, int simIndex){


	//----------------------Simulation options--------------------

	Bool_t two_gaus = true; //True for 2 Gaus resolution fit, False for 1 Gaus fit
	Bool_t KDE_mode = false; // True for KDE as LC template, False for Interpolation as LC template.
	Bool_t include_CollAr = true; //Include collection area
	Bool_t include_EBL = true; //Include EBL absorption
	Bool_t PIC_job = true ; //Changes in path and plots to run the code in PIC
	Bool_t plot_results = false; //True to plot simulated distributions
	Bool_t include_bkg = true; //True to simulate background events

	if(PIC_job == true) plot_results == false; //No plots in PIC



	//---------------------Flare properties (from reading code)----------------------

	const Int_t numberOfTotalEvents = 12235; //(>120 GeV, from data)
	Int_t numberOfBkgEvents = 0.;
	if(include_bkg==1) numberOfBkgEvents = 1466; //(>120 GeV, from data)
	Int_t numberOfEvents = numberOfTotalEvents-numberOfBkgEvents; // Excess events

	cout << "Number of signal events: " << numberOfEvents << endl;
	cout << "Number of background events: " << numberOfBkgEvents << endl;

	Double_t E, t, Es, td;
	Double_t tmin, tmax, Emin, Emax;
	tmin = 0.;
	tmax = hole[numberOfPoints-1];
	Emin=120.; //Change this depending on the Emin cut!
	Emax=50000.;
	Int_t QGorder=1; //Linear or Quadratic case
	Double_t max_inter = 0.;
	Double_t max_inter_bkg = 0.;

	// Create a root File and a TTree
	TTree tevents("tevents","Tree with all events");
	tevents.Branch("E",&E,"E/D"); //True energy
	tevents.Branch("t",&t,"t/D"); //Non-delayed time
	tevents.Branch("Es",&Es,"Es/D"); //Reconstructed energy
	tevents.Branch("td",&td,"td/D"); //Delayed time

	cout << "Emin: " << Emin << " Emax: " << Emax << endl;
	cout << "tmin: " << tmin << " tmax: " << tmax << endl;




	//----------------------------  ENERGY SPECTRUM   ---------------------------

	//Read intrinsic spectrum from Fold ouput file
	TString pathToReadSpec;
	if(PIC_job == false) pathToReadSpec = "./Fold_outputs/";
	else pathToReadSpec = "./";

	const TString fileName = pathToReadSpec + "Output_fold_Zd50_LC120_EPWL.root";
	TFile* Specfile = TFile::Open(fileName, "READ");
	if (!Specfile) {
		cout << "ERROR: file " << Specfile << " not found..." << endl;
	}

	flare_spectrum = dynamic_cast<TF1*>(Specfile->Get("SpectralModel"));
	if (!flare_spectrum) {
		cout << "ERROR: SpectralModel object not found... " << endl;
	}

	Double_t specIntegralSignal=flare_spectrum->Integral(Emin, Emax);
	std::cout << "Integral Spectrum signal = " << specIntegralSignal << std::endl;




	//------------------------Spectrum systematics----------------------------------

	//Change the parameters in the spectrum (Study of systematics)
	//UNCOMMENT WHEN NEEDED

	Int_t numPar = flare_spectrum->GetNpar();
	Double_t parameter[3]= {0.};
	Double_t error_parameter[3] = {0.};
	cout << "Number of parameters: " << numPar << endl;
	for(Int_t i=0; i<numPar; i++){
		parameter[i] = flare_spectrum->GetParameter(i);
		error_parameter[i] = flare_spectrum->GetParError(i);
		cout << "Par " <<i << " value is: " << parameter[i] << " with error: " << error_parameter[i] << endl;
	}

	Int_t numOfPar = 0;
	Double_t new_value = parameter[numOfPar]+error_parameter[numOfPar];
	flare_spectrum->SetParameter(numOfPar, new_value);
	specIntegralSignal=flare_spectrum->Integral(Emin, Emax);
	std::cout << "New Integral Spectrum signal = " << specIntegralSignal << std::endl;

	std::cout << "Intrinsic spectrum loaded" << std::endl;

	//----------------------------------------------------------------------------------

	if(plot_results){
		TCanvas* canvasflare = new TCanvas();
		canvasflare->SetLogx();
		canvasflare->SetLogy();
		flare_spectrum->SetTitle("Intrinsic Spectrum (from Fold)");
		flare_spectrum->SetLineColor(2);
		flare_spectrum->Draw();

		TLegend* leg_spec = new TLegend(0.1,0.7,0.4,0.9);
		leg_spec->AddEntry(flare_spectrum, "Flare (Log parabola)", "l");
		leg_spec->Draw();
	}




	//---------------------------------  LIGHTCURVE  ------------------------------------


	if(KDE_mode == false){

	//-------------------------Read and prepare interpolation-----------------------------

		cout << " The time template is taken from an INTERPOLATION of the data histogram" << endl;

		TString pathToReadInter;
		if(PIC_job == false) pathToReadInter = "../../../ReadEvents/Paper_Inter/";
		else pathToReadInter = "./";
		const TString fileName2 = pathToReadInter + "Inter_Paper_36ONevents_width.root";
		TFile* Vectorfile = TFile::Open(fileName2, "READ");
		if (!Vectorfile) {
				cout << "ERROR: file " << Vectorfile << " not found..." << endl;
			}

		//Read the Tree
		TTree *t1 = (TTree*)Vectorfile->Get("t1");
		Double_t x_point, y_point;
		t1->SetBranchAddress("x",&x_point);
		t1->SetBranchAddress("y",&y_point);

		//Fill vectors
		std::vector<double> x_points;
		std::vector<double> y_points;

		//--First point (Link with the Gaussian)
		Double_t x_value = -90.;
		x_points.push_back(x_value);
		y_points.push_back(gaus_function(x_value));

		//--Interpolation points
		Long64_t nentries = t1->GetEntries();
		for(Long64_t i=0;i<nentries;i++){
			t1->GetEntry(i);
			x_points.push_back(x_point);
			y_points.push_back(y_point);
		}

		//--Last point (Link with the Gaussian)
		Double_t x_value2 = x_points[nentries]-x_value;
		x_points.push_back(x_value2);
		y_points.push_back(gaus_function(x_value2));

		inter_min = x_points[0];
		inter_max = x_points[nentries+1];
		tmin = x_points[1];
		tmax = x_points[nentries];
		cout << "Number of entries: " << nentries << endl;
		cout << "First time point is (tmin): " << x_points[1] << " s" << endl;
		cout << "Last time point is (tmax): " << x_points[nentries] << " s" << endl;
		cout << "New First time point is: " << x_points[0] << " s" << endl;
		cout << "New Last time point is: " << x_points[nentries+1] << " s" << endl;

		//Do and draw Akima interpolation
		int_aspline = new ROOT::Math::Interpolator(x_points, y_points, ROOT::Math::Interpolation::kAKIMA);
		Int_t number_of_points = 8000;
		general =  new TGraph(number_of_points);
		Double_t step = tmax/number_of_points;
		Double_t point = x_points[1];
		Int_t iter=0;
		Double_t scale_factor = (double)numberOfEvents/numberOfTotalEvents;
		cout << "Scale factor for ON interpolation: " << scale_factor << endl;
		while(point+iter*step<=x_points[nentries]){
				if(int_aspline->Eval(point+iter*step)>max_inter)max_inter = int_aspline->Eval(point+iter*step);
				general->SetPoint(iter, point+iter*step, int_aspline->Eval(point+iter*step));
				iter++;
		}

		if(plot_results){
			TCanvas* canvasintersignal =  new TCanvas();
			general->SetLineColor(kBlack);
			general->SetLineWidth(2);
			general->Draw("AL");
		}

		cout << "max_inter for signal: " << max_inter << endl;
		cout << "Interpolation is ready" << endl;
	}


	if(KDE_mode == true){

		//--------------------------Read and prepare KDE function--------------------------

		cout << " The time template is taken from a KDE of the data events" << endl;

		TString pathToReadKDE;
		if(PIC_job == false) pathToReadKDE = "../../../ReadEvents/Thesis_KDE/";
		else pathToReadKDE = "./";
		const TString fileName2 = pathToReadKDE + "TKDE_Kazuma_120_0.5.root";
		TFile* KDEfile = TFile::Open(fileName2, "READ");
		if (!KDEfile) {
				cout << "ERROR: file " << KDEfile << " not found..." << endl;
			}

		KDE_function = (TF1*)KDEfile->Get("TKDE_model");

		if(plot_results){
			TCanvas* canvasinter =  new TCanvas();
			KDE_function->SetLineColor(kRed);
			KDE_function->SetLineWidth(2);
			KDE_function->Draw("");
		}

		Int_t pointKDE= 100;
		Double_t stepKDE= (tmax-tmin)/pointKDE;
		for(Int_t i = 0; i<pointKDE; i++){
			Double_t eval_point = tmin+stepKDE*i;
			if(KDE_function->Eval(eval_point)>KDE_max) KDE_max = KDE_function->Eval(eval_point);
		}

		cout << "Max of KDE signal is: " << KDE_max << endl;
		cout << "KDE function is ready" << endl;

	}


	//----------------------------------BACKGROUND EVENTS-----------------------------------

	Int_t numberOfBinsBkg = 0;
	if(include_bkg==1){

		//--------------------   Background energy distribution   ------------------------
		TString pathToReadBkg;
		if(PIC_job == false) pathToReadBkg = "../../../ReadEvents/Paper_Histos/";
		else pathToReadBkg = "./";
		const TString fileNameBkg = pathToReadBkg + "Histo_36OFF_Paper.root";
		TFile* BKgfile = TFile::Open(fileNameBkg, "READ");
		if (!BKgfile) {
			cout << "ERROR: file " << BKgfile << " not found..." << endl;
		}

		energy_dist_bkg = dynamic_cast<TH1D*>(BKgfile->Get("energy_dist"));
		if (!energy_dist_bkg) {
				cout << "ERROR: energy_dist_bkg object not found... " << endl;
		}

		cout << "----- Reading the OFF histogram... ------" << endl;

		sim_bkg_time = dynamic_cast<TH1D*>(BKgfile->Get("ONvariation"));
		if (!energy_dist_bkg) {
				cout << "ERROR: ONvariation bkg object not found... " << endl;
		}

		Int_t numberEBinsBkg = energy_dist_bkg->GetNbinsX();
		cout << "Number of energy bins in bkg: " << numberEBinsBkg << endl;

		numberOfBinsBkg = sim_bkg_time->GetNbinsX();
		cout << "Number of bins in bkg time template: " << numberOfBinsBkg << endl;
		const Double_t* arrayBinsBkg = sim_bkg_time->GetXaxis()->GetXbins()->GetArray();

		if(plot_results){
			TCanvas* canvasEBkg =  new TCanvas();
			canvasEBkg->SetLogx();
			canvasEBkg->SetLogy();
			energy_dist_bkg->SetTitle("Background energy distribution");
			energy_dist_bkg->SetLineColor(kRed);
			energy_dist_bkg->SetLineWidth(2);
			energy_dist_bkg->Draw("");
		}

		//Generate histos for simulated bkg events
		timedel_bkg = new TH1D("timedel_bkg", "timedel_bkg", numberOfBinsBkg, arrayBinsBkg);
		timedel_bkg_var = new TH1D("timedel_bkg_var", "timedel_bkg_var", numberOfBinsBkg, arrayBinsBkg);


		cout << "Background Energy is ready" << endl;



		//--------------------   Background time distribution   ------------------------


		if(KDE_mode == false){

				//------------------------Background interpolation----------------------

				TString pathToReadInter;
				if(PIC_job == false) pathToReadInter = "../../../ReadEvents/Paper_Inter/";
				else pathToReadInter = "./";
				const TString fileName3 = pathToReadInter + "Inter_Paper_36OFFevents_width.root";
				TFile* Vectorfile = TFile::Open(fileName3, "READ");
				if (!Vectorfile) {
					cout << "ERROR: file " << Vectorfile << " not found..." << endl;
				}

				//Fill vectors
				std::vector<double> x_points_bkg;
				std::vector<double> y_points_bkg;
				//Read the Tree
				TTree *t1 = (TTree*)Vectorfile->Get("t1");
				Double_t x_point, y_point, bin_width;
				t1->SetBranchAddress("x",&x_point);
				t1->SetBranchAddress("y",&y_point);
				t1->SetBranchAddress("width",&bin_width);

				  //--Interpolation points
				  Long64_t nentries = t1->GetEntries();
				  for(Long64_t i=0;i<nentries;i++){
					  t1->GetEntry(i);
					  x_points_bkg.push_back(x_point);
					  y_points_bkg.push_back(y_point);
				  }

				  cout << "Number of entries: " << nentries << endl;
				  cout << "First time point is: " << x_points_bkg[0] << " s" << endl;
				  cout << "First time point value: " << y_points_bkg[0] << " s" << endl;
				  cout << "Last time point is: " << x_points_bkg[nentries-1] << " s" << endl;
				  cout << "Last time point value: " << y_points_bkg[nentries-1] << " s" << endl;



				  //Do and draw Akima interpolation
				  int_aspline_bkg = new ROOT::Math::Interpolator(x_points_bkg, y_points_bkg, ROOT::Math::Interpolation::kAKIMA);

				  Int_t number_of_points = 8000;
				  general_bkg =  new TGraph(number_of_points);
				  general_bkg->SetName("Akima interpolation for Background");
				  general_bkg->SetTitle("Akima interpolation for Background");
				  Double_t step = tmax/number_of_points;
				  Double_t point = x_points_bkg[0];
				  Int_t iter=0;
				  while(point+iter*step<=x_points_bkg[nentries-1]){
					if(int_aspline_bkg->Eval(point+iter*step)>max_inter_bkg)max_inter_bkg = int_aspline_bkg->Eval(point+iter*step);
					general_bkg->SetPoint(iter, point+iter*step, int_aspline_bkg->Eval(point+iter*step)*1./3);
					  iter++;
				  }

				  if(plot_results){
					  TCanvas* canvasinter =  new TCanvas();
					  general_bkg->SetLineColor(kBlack);
					  general_bkg->SetLineWidth(2);
					  general_bkg->Draw("AL");
				  }

				  cout << "max_inter for background: " << max_inter_bkg << endl;
				  cout << "Background Interpolation is ready" << endl;
				}


		if(KDE_mode == true){
			//---------------------------Background KDE-------------------------------

			TString pathToReadKDE;
			if(PIC_job == false) pathToReadKDE = "../../../ReadEvents/Thesis_KDE/";
			else pathToReadKDE = "./";
			const TString fileName2 = pathToReadKDE + "TKDE_Kazuma_120_0.5_OFF.root";
			TFile* KDEfile = TFile::Open(fileName2, "READ");
			if (!KDEfile) {
				cout << "ERROR: file " << KDEfile << " not found..." << endl;
			}

			KDE_function_bkg = (TF1*)KDEfile->Get("TKDE_model");

			if(plot_results){
				TCanvas* canvasinter =  new TCanvas();
				KDE_function_bkg->SetLineColor(kRed);
				KDE_function_bkg->SetLineWidth(2);
				KDE_function_bkg->Draw("");
			}

			Int_t pointKDE= 100;
			Double_t stepKDE= (tmax-tmin)/pointKDE;
			for(Int_t i = 0; i<pointKDE; i++){
				Double_t eval_point = tmin+stepKDE*i;
				if(KDE_function_bkg->Eval(eval_point)>KDE_max_bkg) KDE_max_bkg = KDE_function_bkg->Eval(eval_point);
			}

			cout << "Max of KDE signal is: " << KDE_max_bkg << endl;
			cout << "KDE Background function is ready" << endl;

		}
	}



	//--------------------------INSTRUMENT RESPONSE FUNCTIONS------------------------------

	gRandom->SetSeed(0);

	TString pathToReadIRF;
	if(PIC_job == false) pathToReadIRF = "../../../ReadIRFs/";
	else pathToReadIRF = "./";
	const TString fileName3 = pathToReadIRF + "IRFsMrk421_2014flare_Kazuma_Smooth.root";
	TFile* IRFfile = TFile::Open(fileName3, "READ");
	if (!IRFfile) {
		cout << "ERROR: file " << IRFfile << " not found..." << endl;
	}

	//---------Collection Area---------
	CollAr_Graph = dynamic_cast<TGraph*>(IRFfile->Get("CollArEtrue_Graph"));
	if (!CollAr_Graph) {
			cout << "ERROR: CollArEtrue_Graph object not found... " << endl;
	}

	if(plot_results){
		TCanvas* canvascollar = new TCanvas();
		CollAr_Graph->Draw();
	}


	//----------Double/Simple Gaussian Fit for Resolution----------

	//Const parameter
	TH1D* Const1Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Const1_histo_Etrue"));
	if (!Const1Etrue) {
			cout << "ERROR: Const1Etrue object not found... " << endl;
	}
	TH1D* Const2Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Const2_histo_Etrue"));
	if (!Const1Etrue) {
			cout << "ERROR: Const2Etrue object not found... " << endl;
	}
	TH1D* Const0Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Const0_histo_Etrue"));
	if (!Const0Etrue) {
			cout << "ERROR: Const0Etrue object not found... " << endl;
	}

	if(plot_results){
		TCanvas* const2canvas = new TCanvas();
		const2canvas->SetLogx();
		Const1Etrue->SetLineColor(kBlue);
		Const1Etrue->Draw("");
		Const2Etrue->SetLineColor(kRed);
		Const2Etrue->Draw("same");
		Const0Etrue->SetLineColor(kBlack);
		Const0Etrue->Draw("same");

		TLegend* const2legend = new TLegend(0.1,0.7,0.4,0.9);
		const2legend->AddEntry(Const1Etrue, "Const paratermer, Gauss 1", "l");
		const2legend->AddEntry(Const2Etrue, "Const paratermer, Gauss 2", "l");
		const2legend->AddEntry(Const0Etrue, "Const paratermer, 1 Gauss fit", "l");
		const2legend->Draw("");
	}

	//Mean parameter
	TH1D* Mean1Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Mean1_histo_Etrue"));
	if (!Mean1Etrue) {
			cout << "ERROR: Mean1Etrue object not found... " << endl;
	}
	TH1D* Mean2Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Mean2_histo_Etrue"));
	if (!Mean2Etrue) {
			cout << "ERROR: Mean2Etrue object not found... " << endl;
	}
	TH1D* Mean0Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Mean0_histo_Etrue"));
	if (!Mean0Etrue) {
			cout << "ERROR: Mean0Etrue object not found... " << endl;
	}

	if(plot_results){
		TCanvas* mean2canvas = new TCanvas();
		mean2canvas->SetLogx();
		Mean2Etrue->SetLineColor(kRed);
		Mean2Etrue->Draw("");
		Mean1Etrue->SetLineColor(kBlue);
		Mean1Etrue->Draw("same");
		Mean0Etrue->SetLineColor(kBlack);
		Mean0Etrue->Draw("same");

		TLegend* main2legend = new TLegend(0.1,0.7,0.4,0.9);
		main2legend->AddEntry(Mean1Etrue, "Mean paratermer, Gauss 1", "l");
		main2legend->AddEntry(Mean2Etrue, "Mean paratermer, Gauss 2", "l");
		main2legend->AddEntry(Mean0Etrue, "Mean paratermer, 1 Gauss fit", "l");
		main2legend->Draw("");
	}

	//Sigma parameter
	TH1D* Sigma1Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Sigma1_histo_Etrue"));
	if (!Sigma1Etrue) {
			cout << "ERROR: Sigma1Etrue object not found... " << endl;
	}
	TH1D* Sigma2Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Sigma2_histo_Etrue"));
	if (!Sigma2Etrue) {
			cout << "ERROR: Sigma2Etrue object not found... " << endl;
	}
	TH1D* Sigma0Etrue = dynamic_cast<TH1D*>(IRFfile->Get("Sigma0_histo_Etrue"));
	if (!Sigma0Etrue) {
			cout << "ERROR: Sigma0Etrue object not found... " << endl;
	}

	if(plot_results){
		TCanvas* sigma2canvas = new TCanvas();
		sigma2canvas->SetLogx();
		Sigma2Etrue->SetLineColor(kRed);
		Sigma2Etrue->Draw("");
		Sigma1Etrue->SetLineColor(kBlue);
		Sigma1Etrue->Draw("same");
		Sigma0Etrue->SetLineColor(kBlack);
		Sigma0Etrue->Draw("same");

		TLegend* sigma2legend = new TLegend(0.1,0.7,0.4,0.9);
		sigma2legend->AddEntry(Sigma1Etrue, "Sigma paratermer, Gauss 1", "l");
		sigma2legend->AddEntry(Sigma2Etrue, "Sigma paratermer, Gauss 2", "l");
		sigma2legend->AddEntry(Sigma0Etrue, "Sigma paratermer, 1 Gauss fit", "l");
		sigma2legend->Draw("");
	}

	if(plot_results){

		//-----Resolution test for a given energy------
		Double_t test_energy = 1000; //In GeV

		Double_t Const = Const0Etrue->GetBinContent(Const0Etrue->FindBin(test_energy));
		Double_t mean = Mean0Etrue->GetBinContent(Mean0Etrue->FindBin(test_energy));
		Double_t sigma = Sigma0Etrue->GetBinContent(Sigma0Etrue->FindBin(test_energy));
		Double_t Erec_down_simple = (mean+1)*test_energy-4*test_energy*sigma; //Limits of the integral
		if(Erec_down_simple<0)Erec_down_simple = 0;
		Double_t Erec_up_simple = (mean+1)*test_energy+4*test_energy*sigma;
		cout << "Limits simple Gaus" << endl;
		cout << "UP: " << Erec_down_simple << " DOWN: " << Erec_up_simple << endl;
		TF1* simpleGaus =  new TF1("Simple Gauss", Simple_Gaus, Erec_down_simple, Erec_up_simple, test_energy);
		simpleGaus->SetNpx(10000);

		simpleGaus->SetParameters(Const, (mean+1)*test_energy, sigma*test_energy);

		Double_t const1 = Const1Etrue->GetBinContent(Const1Etrue->FindBin(test_energy));
		Double_t const2 = Const2Etrue->GetBinContent(Const2Etrue->FindBin(test_energy));
		Double_t mean1 = Mean1Etrue->GetBinContent(Mean1Etrue->FindBin(test_energy));
		Double_t mean2 = Mean2Etrue->GetBinContent(Mean2Etrue->FindBin(test_energy));
		Double_t sigma1 = Sigma1Etrue->GetBinContent(Sigma1Etrue->FindBin(test_energy));
		Double_t sigma2 = Sigma2Etrue->GetBinContent(Sigma2Etrue->FindBin(test_energy));
		Double_t Erec_down_doble = (mean1+1)*test_energy-4*test_energy*sigma1; //Limits of the integral
		if(Erec_down_doble<0)Erec_down_doble = 0;
		Double_t Erec_up_doble = (mean1+1)*test_energy+4*test_energy*sigma1;
		cout << "Limits doble Gaus" << endl;
		cout << "UP: " << Erec_down_doble << " DOWN: " << Erec_up_doble << endl;
		TF1* dobleGaus =  new TF1("Doble Gauss", Double_Gaus, Erec_down_doble, Erec_up_doble, 6);
		dobleGaus->SetNpx(10000);

		dobleGaus->SetParameters(const1, (mean1+1)*test_energy, sigma1*test_energy, const2, (mean2+1)*test_energy, sigma2*test_energy);

		TCanvas* energytestcanvas =  new TCanvas();
		energytestcanvas->SetTitle("Resolution and Bias (Etrue = 1 TeV)");
		dobleGaus->SetLineColor(kRed);
		dobleGaus->GetXaxis()->SetTitle("E(GeV)");
		dobleGaus->GetYaxis()->SetTitle("a.u.");
		dobleGaus->Draw("");
		simpleGaus->SetLineColor(kBlue);
		simpleGaus->Draw("same");

		Double_t Es_simple = simpleGaus->GetRandom();
		Double_t Es_doble = dobleGaus->GetRandom();

		TLegend* sigma2legend = new TLegend(0.1,0.7,0.4,0.9);
		sigma2legend->AddEntry(dobleGaus, "Double Gauss", "l");
		sigma2legend->AddEntry(simpleGaus, "Simple Gauss", "l");
		sigma2legend->Draw("");

		cout << "Es_simple: " << Es_simple << " Es_doble: " << Es_doble << endl;


	}


	cout << "IRF are loaded and ready" << endl;



	//--------------------------  SPECTRA + COLLECTION AREA--------------------------------

	folded_spectrum_flare = new TF1("folded_spectrum_flare", FoldingFlare, Emin, Emax, 0);
	folded_spectrum_flare->SetNpx(10000);


	cout << "Spectrum folded with CollAr" << endl;
	cout << "Integral of folded spectrum: " << folded_spectrum_flare->Integral(Emin,Emax) << endl;
	cout << "Expected events: " << folded_spectrum_flare->Integral(Emin,Emax)*(tmax-tmin) << endl;


	//----------------------------------EBL ABSORPTION------------------------------------

	Double_t redshift = 0.031;
	Float_t z[39] = {0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895, 0.11684211, 0.13210526, 0.14736842, 0.16263158, 0.17789474, 0.19315789, 0.20842105, 0.22368421, 0.23894737, 0.25421053, 0.26947368, 0.28473684, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.2, 1.4, 1.6, 1.8, 2.};

	//Read dominguez table
	ifstream f("tau_dominguez11.dat");
	if (!f.is_open())
	{
		cout << "Cannot open input file" << endl;
		return;
	}
	Char_t ch[1000];
	Double_t data[50];
	TGraph2D* tau_vs_z_vs_e = new TGraph2D;
	Int_t ipoint = 0;

	// Set optical depth = 0 at all redshifts for E=0GeV
	tau_vs_z_vs_e->SetPoint(ipoint++, 0., 0., 0.);
	for (Int_t iz = 0; iz < 39; iz++) tau_vs_z_vs_e->SetPoint(ipoint++, 0., z[iz], 0.);
	Int_t i_energy = 1;
	for (;;)  // LOOP OVER ENERGIES (over lines)
	{
		f.getline(ch, 1000);
		if (f.eof())break; //End of file
		if (ch[0] == '#')continue; //Comment
		TString str = ch;
		TObjArray* arr;
		arr = str.Tokenize(" "); //Cut line in "words"
		for (Int_t i = 0; i < arr->GetEntries(); i++)  // Get all optical depths for this energy (Reading one line)
		{
			data[i] = ((TObjString*)arr->At(i))->GetString().Atof();
		}
		// Fill the optical depth (0) for redshift 0
		//data[0] is the energy (TeV) and the rest data[i] are the optical depth.
		tau_vs_z_vs_e->SetPoint(ipoint++, log10(data[0]*1.e3), 0., 0.); // tau=0 for z=0
		for (Int_t i = 1; i < arr->GetEntries(); i++)
		{
			tau_vs_z_vs_e->SetPoint(ipoint++, log10(data[0]*1.e3), z[i-1], data[i]);
			// E in GeV!!
			//	  cout << log10(data[0]*1.e3) << " = log10(E/GeV),  z = " << z[i-1] << ",  tau = " << data[i] << endl;
		}
		i_energy++;
	}

	tau_vs_z_vs_e->SetNpx(i_energy);
	tau_vs_z_vs_e->SetNpy(40);

	//Extract the optical depth by interpolation
	transvse = new TGraph();
	ipoint=0;
	for(Double_t log10e = 0.; log10e < log10(3.e4); log10e += 0.05)
	{
		transvse->SetPoint(ipoint, log10e, exp(-1.*tau_vs_z_vs_e->Interpolate(log10e, redshift)));
		ipoint++;
	}

	if(plot_results){
		TCanvas* transcanvas = new TCanvas;
		transcanvas->SetLogy();
		transvse->GetXaxis()->SetTitle("Log10(E/GeV)");
		transvse->GetYaxis()->SetTitle("Transmittance");
		transvse->SetLineColor(2);
		transvse->SetLineWidth(2);
		transvse->Draw("alp");
	}

	//--------------------  SPECTRA + COLLECTION AREA + EBL  ------------------------------

	TF1* EBL_spectrum_flare = new TF1("EBL_spectrum_flare", EBLFlare, Emin, Emax, 0);
	EBL_spectrum_flare->SetNpx(10000);

	if(plot_results){
		TCanvas* foldedcanvas = new TCanvas();
		foldedcanvas->SetLogx();
		foldedcanvas->SetLogy();
		folded_spectrum_flare->GetXaxis()->SetTitle("Energy (GeV)");
		folded_spectrum_flare->GetYaxis()->SetTitle("GeV^{-1}s^{-2}");
		folded_spectrum_flare->SetTitle("Spectra");
		folded_spectrum_flare->SetTitle("Spectra");
		folded_spectrum_flare->SetLineColor(2);
		folded_spectrum_flare->SetLineStyle(2);
		folded_spectrum_flare->Draw();
		flare_spectrum->SetLineColor(1);
		flare_spectrum->SetLineStyle(1);
		flare_spectrum->Draw("same");
		EBL_spectrum_flare->SetLineColor(3);
		EBL_spectrum_flare->SetLineStyle(2);
		EBL_spectrum_flare->Draw("same");


		TLegend* leg_spec_fold = new TLegend(0.1,0.7,0.4,0.9);
		leg_spec_fold->AddEntry(flare_spectrum, "Flare (Log parabola)", "l");
		leg_spec_fold->AddEntry(folded_spectrum_flare, "Flare+CollAr", "l");
		leg_spec_fold->AddEntry(EBL_spectrum_flare, "Flare+CollAr+EBL", "l");
		leg_spec_fold->Draw();
	}

	cout << "EBL applied to the folded spectrum" << endl;
	cout << "Integral after EBL: " << EBL_spectrum_flare->Integral(Emin,Emax) << endl;
	cout << "Expected events: " << EBL_spectrum_flare->Integral(Emin,Emax)*(tmax-tmin) << endl;


	//------------------------------  QG LIV DELAY INPUTS  ------------------------------
	Double_t z_source=0.031;
	Double_t H0=70.;
	Double_t PlanckScale=1.2*TMath::Power(10,19);
	Double_t PcTom=3.085*TMath::Power(10,19);
	Double_t QGfixedPart=z_source*PcTom/H0;
	std::cout << "alpha: " << alpha << std::endl;
	std::cout << "QG fixed part = " << QGfixedPart << std::endl;


	//------------------Histograms for simulated data (time and energy)--------------------

	//Equal binning in LogE
	double XX_step = 0.1;
	int n = ceil((log10(Emax)-log10(Emin))/XX_step);
	double *XX = new double[n+1];
	for(int i=0;i<n+1;i++)
	{
		XX[i]= pow(10,(log10(Emin)+i*XX_step));
	}

	//Read ON time histogram
	cout << "----- Reading the ON histogram... ------" << endl;
	TFile* HistoFile;
	if(PIC_job == false) HistoFile = TFile::Open("../../../ReadEvents/Paper_Histos/Histo_36ON_Paper.root", "READ");
	else HistoFile = TFile::Open("./Histo_36ON_Paper.root", "READ");
	if (!HistoFile) {
		cout << "ERROR: file " << HistoFile << " not found..."	<< endl;
	}
	TH1D* lightcurve = dynamic_cast<TH1D*>(HistoFile->Get("ONvariation"));
	if (!lightcurve) {
		cout << "ONvariation not found" << endl;
	}

	Int_t numberOfBins = lightcurve->GetNbinsX();
	cout << "Number of bins in initial histo: " << numberOfBins << endl;
	const Double_t* arrayBins = lightcurve->GetXaxis()->GetXbins()->GetArray();


	//Intrinsic time/energy histos
	sim_flare_time = new TH1D("sim_flare_time", "sim_flare_time", numberOfBins, arrayBins);
	sim_flare_time_var = new TH1D("sim_flare_time_var", "sim_flare_time_var", numberOfBins, arrayBins);
	TH1D* sim_flare_energy = new TH1D("sim_flare_energy", "sim_flare_energy", 2000, Emin, Emax);

	//Reconstructed time/energy histos
	TH1D* timedel = new TH1D("timedel", "timedel", numberOfBins, arrayBins);
	TH1D* timedel_var = new TH1D("timedel_var", "timedel_var", numberOfBins, arrayBins);
	TH1D* energyrec = new TH1D("energyrec", "energyrec", 2000, Emin, Emax);
	TH1D* energyrec_bkg = new TH1D("energyrec_bkg", "energyrec_bkg", 2000, Emin, Emax);



	//------------------------------Read Real Data (for comparison)-----------------------------------------

	int Npoints = 0;
	TH1D* data_time;
	TH1D* data_time_var;
	TH1D* data_energy;
	if(PIC_job == false){ //no comparison in PIC
		ifstream in;
		in.open("../../../ReadEvents/Paper_EventLists/Event_Paper_Zd50_ON.txt");
		Double_t tp[20000],Ep[20000];
		Double_t v1,v2;

		Double_t max_time = 0.;
		while(1)
		{
			in >> v1 >> v2;
			if (!in.good()) break;
			tp[Npoints] = v1;
			Ep[Npoints] = v2;
			if (tp[Npoints] > max_time) max_time = tp[Npoints];
			Npoints++;
		}
		cout << "Npoints(Number of events): " << Npoints << endl;
		cout << "Max_time: " << max_time << endl;
		if(KDE_mode == true) tmax = max_time;

		//Histos for data
		data_time = new TH1D("Time distribution", "Time distribution", numberOfBins, arrayBins);
		data_time_var = new TH1D("Time  variation distribution", "Time variation distribution", numberOfBins, arrayBins);
		data_energy = new TH1D("Energy distribution", "Energy distribution", 2000, Emin, Emax);

		for(int i = 0; i<Npoints; i++){
			data_time->Fill(tp[i]);
			data_energy->Fill(Ep[i]);
		}
		for(Int_t bin=1; bin<numberOfBins; bin++){
			data_time_var->SetBinContent(bin, data_time->GetBinContent(bin)/data_time->GetBinWidth(bin));
			data_time_var->SetBinError(bin, data_time->GetBinError(bin)/data_time->GetBinWidth(bin));
		}
	}



	//-------------------------------Signal Simulation (Energy and time together) ------------------------------------------------

	//LIV delay injection
	Double_t QGDelay = 0.;
	  if(QGorder==1){ //Linear
		QGDelay= QGfixedPart*alpha/PlanckScale;
	  }
	  if(QGorder==2){ //Quadratic
		QGDelay= QGfixedPart*alpha/PlanckScale/1000.;
	  }
	std::cout << "Injected tau: " << QGDelay << std::endl;
	std::cout << "Injected delay (in s) for 1 TeV (Linear): " << QGDelay*1000 << std::endl;
	std::cout << "Injected delay (in s) for 1 TeV (Quadratic): " << QGDelay*1000*1000 << std::endl;


	//File to save the simulated events
	TFile* inter_outputfile = new TFile(Form("Sys_Spec_%i_Paper_Par0Up_%i.root", (int)alpha, simIndex),"recreate");


	//ENERGY INPUTS (Leyre version: Sim in E, Including CollAr and EBL)
	Int_t sim_events = 0;
	Int_t acceptEvents = 0;
	Int_t plotEvents = 0;

	//TIME INPUTS (Hit or miss with Interpolation/KDE and Gaussian smearing)
	Double_t extra_time_edges = 1000.; //in seconds
	Double_t x_min = tmin-extra_time_edges;
	Double_t x_max = tmax+extra_time_edges;
	Double_t y_min = 0.;
	Double_t y_max = 0.; //Check the corresponding function (KDE or Inter)
	if(KDE_mode == false) y_max = max_inter;
	else y_max = KDE_max;
	cout << "Time: xmin: " << x_min << " xmax: " << x_max << endl;
	cout << "Time: ymin: " << y_min << " ymax: " << y_max << endl;
	TRandom3* ranx = new TRandom3();
	ranx->SetSeed(0);
	TRandom3* rany = new TRandom3();
	rany->SetSeed(0);
	Double_t random_x, random_y;




	//---------------------------- SIMULATION STARTS HERE ---------------------------------
	cout << "Starting simulation... " << endl;

	while(plotEvents<numberOfEvents){

		//----- Time simulation ----

		random_x = ranx->Uniform(x_min,x_max); //Coordinate in x (t)
		random_y = rany->Uniform(y_min,y_max); //Coordinate in y (t)

		sim_events++;
		Bool_t y_condition_t= 0;
		Bool_t wobble_condition= 0;

		if(KDE_mode == false) y_condition_t = random_y > lc_function(random_x);
		else y_condition_t = random_y > KDE_function->Eval(random_x);

		wobble_condition = evaluate_wobble_condition(random_x);

//		cout << "t_condition: " << y_condition_t << endl;
//		cout << "wobble_condition: " << wobble_condition << endl;

		if(y_condition_t || wobble_condition)continue; //Discard event
		else{

			//----- Energy simulation -----
			if(include_CollAr == false && include_EBL == false) E = flare_spectrum->GetRandom(Emin,Emax);
			if(include_CollAr == true && include_EBL == false) E = folded_spectrum_flare->GetRandom(Emin,Emax);
			if(include_CollAr == true && include_EBL == true) E = EBL_spectrum_flare->GetRandom(Emin,Emax);
			t = random_x;
			acceptEvents++;

			//Smear Energy (with the resolution function)
			if(two_gaus == false){
				Double_t Const = Const0Etrue->GetBinContent(Const0Etrue->FindBin(E));
				Double_t mean = Mean0Etrue->GetBinContent(Mean0Etrue->FindBin(E));
				Double_t sigma = Sigma0Etrue->GetBinContent(Sigma0Etrue->FindBin(E));
				Double_t Erec_down_simple = (mean+1)*E-4*E*sigma; //Limits of the integral
				if(Erec_down_simple<0)Erec_down_simple = 0;
				Double_t Erec_up_simple = (mean+1)*E+4*E*sigma;
				TF1* simpleGaus =  new TF1("Simple Gauss", Simple_Gaus, Erec_down_simple, Erec_up_simple, 3);
				simpleGaus->SetNpx(10000);

				simpleGaus->SetParameters(Const, (mean+1)*E, sigma*E);

				Es = simpleGaus->GetRandom(Erec_down_simple, Erec_up_simple);

				if(acceptEvents == 200 && plot_results){ //Plot a case for cross check
					cout << "E: " << E << " Es: " << Es << endl;
					TCanvas* simpleGausCanvas = new TCanvas();
					simpleGaus->DrawCopy();
				}
				simpleGaus->Delete();

			}
			else{
				Double_t const1 = Const1Etrue->GetBinContent(Const1Etrue->FindBin(E));
				Double_t const2 = Const2Etrue->GetBinContent(Const2Etrue->FindBin(E));
				Double_t mean1 = Mean1Etrue->GetBinContent(Mean1Etrue->FindBin(E));
				Double_t mean2 = Mean2Etrue->GetBinContent(Mean2Etrue->FindBin(E));
				Double_t sigma1 = Sigma1Etrue->GetBinContent(Sigma1Etrue->FindBin(E));
				Double_t sigma2 = Sigma2Etrue->GetBinContent(Sigma2Etrue->FindBin(E));
				Double_t Erec_down_doble = (mean1+1)*E-4*E*sigma1; //Limits of the integral
				if(Erec_down_doble<0)Erec_down_doble = 0;
				Double_t Erec_up_doble = (mean1+1)*E+4*E*sigma1;

				TF1* dobleGaus =  new TF1("Doble Gauss", Double_Gaus, Erec_down_doble, Erec_up_doble, 6);
				dobleGaus->SetNpx(10000);

				dobleGaus->SetParameters(const1, (mean1+1)*E, sigma1*E, const2, (mean2+1)*E, sigma2*E);

				Es = dobleGaus->GetRandom(Erec_down_doble, Erec_up_doble);

				if(acceptEvents == 20 && plot_results){ //Plot a case for cross check
					cout << "E: " << E << " Es: " << Es << endl;
					TCanvas* dobleGausCanvas = new TCanvas();
					dobleGaus->DrawCopy();
				}
				dobleGaus->Delete();

			}



			//----- Add LIV delay -----
			if(QGorder == 1) td=t+QGDelay*E;
			else if(QGorder == 2) td=t+QGDelay*E*E;

			//Check if the delayed event is inside a wobble hole.
			wobble_condition = evaluate_wobble_condition(td);

			//Selection of events
			if(Es>Emin && Es<Emax && td>tmin && td<tmax && wobble_condition==0){
				tevents.Fill();
				sim_flare_energy->Fill(E);
				sim_flare_time->Fill(t);
				energyrec->Fill(Es);
				timedel->Fill(td);
				plotEvents++;
			}
		}
	}

	cout << "Simulated events: " << sim_events << endl;
	cout << "Accepted events: " << acceptEvents << endl;
	cout << "Plot events: " << plotEvents << endl;
	cout << "Tree length: " << tevents.GetEntries() << endl;

	//Fill variation histogram
	for(Int_t bin=1; bin<numberOfBins; bin++){
		sim_flare_time_var->SetBinContent(bin, sim_flare_time->GetBinContent(bin)/sim_flare_time->GetBinWidth(bin));
		sim_flare_time_var->SetBinError(bin, sim_flare_time->GetBinError(bin)/sim_flare_time->GetBinWidth(bin));
		timedel_var->SetBinContent(bin, timedel->GetBinContent(bin)/timedel->GetBinWidth(bin));
		timedel_var->SetBinError(bin, timedel->GetBinError(bin)/timedel->GetBinWidth(bin));

	}



	//-------------------------------Bkg Simulation (Energy and time together) ------------------------------------------------

	Int_t sim_bkg = 0;
	Int_t accept_bkg = 0;
	Int_t plot_bkg = 0;

	if(include_bkg ==1){

		cout << "Starting BACKGROUND simulation... " << endl;

		Double_t x_min_bkg = tmin;
		Double_t x_max_bkg = tmax;
		Double_t y_min_bkg = 0.;
		Double_t y_max_bkg = 0.; //Check the corresponding function (KDE or Inter)
		if(KDE_mode == false) y_max_bkg = max_inter_bkg;
		else y_max_bkg = KDE_max_bkg;
		cout << "BKG Time: xmin: " << x_min_bkg << " xmax: " << x_max_bkg << endl;
		cout << "BKG Time: ymin: " << y_min_bkg << " ymax: " << y_max_bkg << endl;
		TRandom3* ranx_bkg = new TRandom3();
		ranx_bkg->SetSeed(0);
		TRandom3* rany_bkg = new TRandom3();
		rany_bkg->SetSeed(0);
		Double_t random_x_bkg, random_y_bkg;
		cout << "Plot events: " << plotEvents << endl;
		cout << "Total events: " << numberOfTotalEvents << endl;
		while(plotEvents<numberOfTotalEvents){
			E = energy_dist_bkg->GetRandom();
			sim_events++;
			sim_bkg++;
			Bool_t y_condition_t= 0;
			Bool_t wobble_condition= 0;
			random_x_bkg = ranx_bkg->Uniform(x_min_bkg,x_max_bkg); //Coordinate in x (t)
			random_y_bkg = rany_bkg->Uniform(y_min_bkg,y_max_bkg); //Coordinate in y (t)

			if(KDE_mode == false) y_condition_t = random_y_bkg > lc_function_background(random_x_bkg);
			else y_condition_t = random_y_bkg > KDE_function->Eval(random_x_bkg);
			wobble_condition = evaluate_wobble_condition(random_x_bkg);

	//		cout << "t_condition: " << y_condition_t << endl;
	//		cout << "wobble_condition: " << wobble_condition << endl;

			if(y_condition_t || wobble_condition)continue;
			else{
				t=td=random_x_bkg;
				Es=E; //Background directly simulated from the measured energy distribution-> No smearing
				acceptEvents++;
				accept_bkg++;
				if(Es>Emin){
					tevents.Fill();
//					sim_flare_energy->Fill(Es); //The total signal includes the bkg
//					sim_flare_time->Fill(td);
//					energyrec->Fill(Es);
					energyrec_bkg->Fill(Es);
//					timedel->Fill(td);
					timedel_bkg->Fill(td);
					plotEvents++;
					plot_bkg++;
				}
			}
		}

		//Fill variation histogram
		for(Int_t bin=1; bin<numberOfBinsBkg; bin++){
			timedel_bkg_var->SetBinContent(bin, timedel_bkg->GetBinContent(bin)/timedel_bkg->GetBinWidth(bin));
			timedel_bkg_var->SetBinError(bin, timedel_bkg->GetBinError(bin)/timedel_bkg->GetBinWidth(bin));

		}
	}

	cout << energyrec_bkg->GetEntries() << endl;

	cout << "After BACKGROUND " << endl;
	cout << "Simulated bkg: " << sim_bkg << endl;
	cout << "Accepted bkg: " << accept_bkg << endl;
	cout << "Plot bkg: " << plot_bkg << endl;
	cout << "Simulated total events: " << sim_events << endl;
	cout << "Accepted total events: " << acceptEvents << endl;
	cout << "Plot total events: " << plotEvents << endl;
	cout << "Tree length: " << tevents.GetEntries() << endl;





	if(plot_results){

		//PLOTS ENERGY AND TIME
		KDE_function_scaled = new TF1("KDE_scaled", ScaleKDE, tmin, tmax);
		data_energy->SetStats(0);
		TCanvas* energysimcanvas = new TCanvas();
		energysimcanvas->SetLogx();
		energysimcanvas->SetLogy();
		data_energy->SetLineColor(1);
		data_energy->GetXaxis()->SetTitle("Energy(GeV)");
		data_energy->SetLineWidth(2);
		data_energy->Draw("EP");
		sim_flare_energy->SetStats(0);
		sim_flare_energy->SetLineColor(2);
		sim_flare_energy->Draw("sameEP");
		energyrec->SetLineColor(4);
		energyrec->Draw("sameEP");
		bin_value = sim_flare_energy->GetBinContent(sim_flare_energy->FindBin(E0));
		cout << "bin value: " << bin_value << endl;

		//This is needed to re-scale the energy template before plotting it or using it in the chi-squared test
		if(include_CollAr == false && include_EBL == false) energy_template = (TF1*)flare_spectrum->Clone("energy_template");
		if(include_CollAr == true && include_EBL == false) energy_template = (TF1*)folded_spectrum_flare->Clone("energy_template");
		if(include_CollAr == true && include_EBL == true) energy_template = (TF1*)EBL_spectrum_flare->Clone("energy_template");
		factor_scale = bin_value/energy_template->Eval(E0);
		cout << "factor scale" << factor_scale << endl;

		energy_template_scaled = new TF1("energy_template_scaled", ScaleEnergyFunction, Emin, Emax);
		energy_template_scaled->SetLineColor(3);
		energy_template_scaled->SetLineStyle(1);
		energy_template_scaled->SetLineWidth(2);
		energy_template_scaled->Draw("same");


		TLegend* leg_energysim = new TLegend(0.7,0.7,0.9,0.9);
		leg_energysim->AddEntry(data_energy, "Data energy dist.", "l");
		leg_energysim->AddEntry(sim_flare_energy, "Simulated energy dist.", "l");
		leg_energysim->AddEntry(energyrec, "Smeared energy dist.", "l");
		leg_energysim->Draw();


		TCanvas* timesim = new TCanvas();
		data_time_var->SetStats(0);
		data_time_var->SetLineColor(1);
		data_time_var->SetLineWidth(2);
		data_time_var->Draw("EP");
		sim_flare_time_var->SetStats(0);
		sim_flare_time_var->SetLineColor(2);
		sim_flare_time_var->Draw("sameEP");
		timedel_var->SetLineColor(4);
		timedel_var->Draw("sameEP");
		if(KDE_mode == false){
			general->SetLineColor(3);
			general->Draw("same");
		}
		else{
			KDE_function_scaled->SetLineColor(3);
			KDE_function_scaled->Draw("same");
		}


		TLegend* leg_timesim = new TLegend(0.7,0.7,0.9,0.9);
		leg_timesim->AddEntry(data_time_var, "Data time dist.", "l");
		leg_timesim->AddEntry(sim_flare_time_var, "Simulated time dist.", "l");
		leg_timesim->AddEntry(timedel_var, "Delayed time dist.", "l");
		if(KDE_mode == false)leg_timesim->AddEntry(general, "Time Interpolation", "l");
		else leg_timesim->AddEntry(KDE_function_scaled, "KDE function", "l");
		leg_timesim->Draw();

		if(include_bkg==1){
			TCanvas* energysimcanvas2 = new TCanvas();
			energysimcanvas2->SetLogx();
			energysimcanvas2->SetLogy();
			energyrec_bkg->SetLineColor(2);
			energyrec_bkg->Draw("EP");
			energy_dist_bkg->Draw("same");

			TCanvas* timesim2 = new TCanvas();
			timedel_bkg_var->SetLineColor(4);
			timedel_bkg_var->Draw("EP");
			general_bkg->SetLineColor(3);
			general_bkg->Draw("same");
		}


	}




	//--------------------------------ChiSquare tests (Time and energy)-----------------------

	if(PIC_job == false && plot_results == true){

		//TIME
		Double_t chiSquare = 0;
		TH1D* 	residualhisto = new TH1D("Residuals", "Residuals", numberOfBins, arrayBins);
		Double_t residual[500]={0};
		cout << "Number of Bins: " << numberOfBins << endl;
		Int_t chi_counter = 0;

		for(int i = 1; i<=numberOfBins; i++){
			Double_t observedValue = sim_flare_time_var->GetBinContent(i);
			Double_t startBin = sim_flare_time_var->GetBinLowEdge(i);
			Double_t endBin = sim_flare_time_var->GetBinLowEdge(i+1);
			Double_t errorBin = sim_flare_time_var->GetBinError(i);
			Double_t widthBin = sim_flare_time_var->GetBinWidth(i);

			if(observedValue!=0){
			Double_t expectedValue = 0.;
			if(KDE_mode == false) expectedValue = int_aspline->Integ(startBin, endBin)*0.9/widthBin;
			if(KDE_mode == true) expectedValue = KDE_function_scaled->Integral(startBin, endBin)/widthBin;
			chiSquare+= (TMath::Power(observedValue-expectedValue,2))/expectedValue;
			residual[i-1] = (observedValue-expectedValue)/errorBin;
			residualhisto->SetBinContent(i, residual[i-1]);
//			cout << "observedValue: " << observedValue << endl;
//			cout << "expectedValue: " << expectedValue << endl;
//			cout << "chiSquare: " << chiSquare << endl;
			residual[i-1] = 0;
			chi_counter++;
			}
		}

		cout << "Chi counter time: " << chi_counter << endl;

		if(plot_results){
			TCanvas* canvas = new TCanvas("canvas", "canvas", 2000,2000);
			TPad* upper =  new TPad("upper", "upper",0.005, 0.3525, 0.995, 0.995);
			TPad* lower =  new TPad("lower", "lower", 0.005, 0.005, 0.995, 0.3475);
			canvas->cd();
			upper->Draw();
			lower->Draw();
			upper->cd();
			sim_flare_time_var->SetTitle("Chi-Squared Time");
			sim_flare_time_var->Draw("ehist");
		//	general->SetLineColor(kBlue);
			if(KDE_mode == false)general->Draw("same");
			else KDE_function_scaled->Draw("same");
			lower->cd();
			residualhisto->SetStats(0);
			residualhisto->GetYaxis()->SetTitle("Sigma");
			residualhisto->GetXaxis()->SetTitle("Time(s)");
			residualhisto->SetLineColor(4);
			residualhisto->SetLineWidth(2);
			residualhisto->Draw("EB");
		}

		//ENERGY
		chiSquare = 0;
		TH1D* residualhisto2 = new TH1D("Residuals2", "Residuals2", 2000, Emin, Emax);
		chi_counter=0;

		for(int i = 1; i<=2000; i++){
			Double_t observedValue = sim_flare_energy->GetBinContent(i);
			Double_t startBin = sim_flare_energy->GetBinLowEdge(i);
			Double_t endBin = sim_flare_energy->GetBinLowEdge(i+1);
			Double_t errorBin = sim_flare_energy->GetBinError(i);
			Double_t widthBin = sim_flare_energy->GetBinWidth(i);
	//		cout << "Low Edge energy: " << startBin << endl;
	//		cout << "High Edge energy: " << endBin << endl;
	//		cout << "Medium energy: " << midBin << endl;
	//		cout << "ErrorBin: " << errorBin << endl;

			if(observedValue!=0){
				Double_t expectedValue = energy_template_scaled->Integral(startBin,endBin)/widthBin;
				chiSquare+= (TMath::Power(observedValue-expectedValue,2))/expectedValue;
				residual[i-1] = (observedValue-expectedValue)/errorBin;
				residualhisto2->SetBinContent(i, residual[i-1]);
	//			cout << "residual: " << residual[i-1] << endl;
//				cout << "observedValue: " << observedValue << endl;
//				cout << "expectedValue: " << expectedValue << endl;
//				cout << "chiSquare: " << chiSquare << endl;
				residual[i-1] = 0;
				chi_counter++;
			}
		}

		cout << "Chi counter energy: " << chi_counter << endl;


		if(plot_results){
			TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 2000,2000);
			TPad* upper2 =  new TPad("upper2", "upper2",0.005, 0.3525, 0.995, 0.995);
			TPad* lower2 =  new TPad("lower2", "lower2", 0.005, 0.005, 0.995, 0.3475);
			canvas2->cd();
			upper2->Draw();
			lower2->Draw();
			upper2->cd();
			upper2->SetLogx();
			upper2->SetLogy();
			sim_flare_energy->SetTitle("Chi-Squared Energy");
			sim_flare_energy->Draw("ehist");
			energy_template_scaled->Draw("same");
			lower2->cd();
			lower2->SetLogx();
			residualhisto2->SetStats(0);
			residualhisto2->GetYaxis()->SetTitle("Sigma");
			residualhisto2->GetXaxis()->SetTitle("Energy(GeV)");
			residualhisto2->SetLineColor(4);
			residualhisto2->SetLineWidth(2);
			residualhisto2->Draw("EB");
		}
	}


	//--------------------------------Write root file with the tree of events-----------------------


	cout << "I write the tree of events" << endl;
	tevents.Write("Event_tree.root");


}

void DrawLightcurve(Double_t timemin, Double_t timemax){
	Int_t numberOfPoints = 100;
	Double_t step = (timemax-timemin)/numberOfPoints;
	TGraph* LC = new TGraph (numberOfPoints);
	Double_t time = 0.;

	for(Int_t i=0; i<numberOfPoints; i++){
		time = timemin + i*step;
		LC->SetPoint(i, time, lc_function(time));
	}

	TCanvas* canvas = new TCanvas();
	canvas->cd(1);
	LC->Draw();
}


