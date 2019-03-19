/*
 * Sim_alpha_Mrk421.C
 *
 *  Created on: Jan 20, 2017
 *      Author: lnogues
 *
 *      Modified version to simulate using a Interpolation or a function. Oct2017.
 *      Modified version to be launched with a job.
 *
 */

#include<iostream>
#include<fstream>
#include<istream>
#include<TGraph.h>
#include<TGraph2D.h>
#include<TF1.h>
#include<TFile.h>
#include<TH1D.h>
#include<TStyle.h>
#include<TMath.h>
#include<TRandom.h>
#include<TCanvas.h>
#include<TLegend.h>
#include <TVector.h>
#include <TVectorT.h>
#include<TMath.h>
#include <Math/InterpolationTypes.h>
#include <Math/Interpolator.h>
#include "TTree.h"
#include <TROOT.h>
#include <TRandom3.h>


using namespace std;

TGraph* CollAr_Graph;
TF1* flare_spectrum;
TF1* folded_spectrum_flare;
//Graph for EBL application
TGraph *transvse;
Double_t inter_min = 0.;
Double_t inter_max = 0.;
ROOT::Math::Interpolator *int_aspline = NULL;
bool fixed_resolution = true;
//Spectrum variables
const Double_t ln10 = 2.302585093;
const Double_t E0 = 364.;


Double_t FoldingFlare(Double_t* x, Double_t* par){
	Double_t flareValue = flare_spectrum->Eval(x[0]);
	Double_t areaValue = CollAr_Graph->Eval(x[0]);
	return flareValue*areaValue;
}

Double_t EBLFlare(Double_t* x, Double_t* par){
	Double_t flareValue = folded_spectrum_flare->Eval(x[0]);
	Double_t transValue = transvse->Eval(log10(x[0]));
	return flareValue*transValue;
}

Double_t lc_function(Double_t t){

	Double_t cons_value = 31.54;
	Double_t weight_value = 84.50;
	Double_t peak_value = 4294.06;
	Double_t sigma_left = 2350.84;
	Double_t sigma_right = 1561.74;
	Double_t gaus_value = 0;

	//Compute values of the gaussian
	if(t<peak_value)
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_left,2));
	else
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_right,2));


	//Return Interpolation or Gaussian
    if(t<inter_min || t>inter_max)
    	return gaus_value;
    else
		return int_aspline->Eval(t);
}

Double_t gaus_function(Double_t t){

	Double_t cons_value = 31.54;
	Double_t weight_value = 84.50;
	Double_t peak_value = 4294.06;
	Double_t sigma_left = 2350.84;
	Double_t sigma_right = 1561.74;
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

void Sim_alpha_Mrk421(double alpha, int simIndex){



	//---------------------Flare properties (from reading code)----------------------
	const Int_t numberOfEvents = 12994; //(>100 GeV, from Mireia data)
	Double_t E[numberOfEvents],t[numberOfEvents], Es[numberOfEvents],td[numberOfEvents];
	Double_t tmin, tmax, Emin, Emax;
	Emin=100.;
	Emax=60000.;
	Int_t QGorder=1; //Linear or Quadratic case





	//----------------------Spectrum reading and integrating---------------------------

	//Read intrinsic spectrum from Fold ouput file
	TString pathToReadSpec = "./Fold_outputs/";
//	TString pathToReadSpec = "./";
	const TString fileName = pathToReadSpec + "Output_fold_LP.root";
	TFile* Specfile = TFile::Open(fileName, "READ");
	if (!Specfile) {
		cout << "ERROR: file " << Specfile << " not found..." << endl;
	}


	flare_spectrum = dynamic_cast<TF1*>(Specfile->Get("SpectralModel"));
	if (!flare_spectrum) {
		cout << "ERROR: SpectralModel object not found... " << endl;
	}

	flare_spectrum->SetRange(Emin,Emax);
	Double_t specIntegralSignal=flare_spectrum->Integral(Emin, Emax);
	std::cout << "Integral Spectrum signal = " << specIntegralSignal << std::endl;
	std::cout << "Intrinsic spectrum loaded" << std::endl;



//	TCanvas* canvasflare = new TCanvas();
//	canvasflare->SetLogx();
//	canvasflare->SetLogy();
//	flare_spectrum->SetTitle("Intrinsic Spectrum (from Fold)");
//	flare_spectrum->SetLineColor(2);
//	flare_spectrum->Draw();
//
//	TLegend* leg_spec = new TLegend(0.1,0.7,0.4,0.9);
//	leg_spec->AddEntry(flare_spectrum, "Flare (Log parabola)", "l");
//	leg_spec->Draw();




	//-----------------------------------Read and prepare interpolation--------------------------------------
	TString pathToReadInter = "../../../ReadEvents/LC_fits/";
//	TString pathToReadInter = "./";
	const TString fileName2 = pathToReadInter + "LC_Inter_Tree_Mireia.root";
	TFile* Vectorfile = TFile::Open(fileName2, "READ");
	if (!Vectorfile) {
			cout << "ERROR: file " << Vectorfile << " not found..." << endl;
		}

	//Read the Tree
	TTree *t1 = (TTree*)Vectorfile->Get("t1");
	Double_t x_point, y_point, error_point;
   t1->SetBranchAddress("x",&x_point);
   t1->SetBranchAddress("y",&y_point);
   t1->SetBranchAddress("error",&error_point);

	//Fill vectors
	std::vector<double> x_points;
	std::vector<double> y_points;
	std::vector<double> error_points;
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
	  error_points.push_back(error_point);
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
  cout << "First time point is: " << x_points[1] << " s" << endl;
  cout << "Last time point is: " << x_points[nentries] << " s" << endl;
  cout << "New First time point is: " << x_points[0] << " s" << endl;
  cout << "New Last time point is: " << x_points[nentries+1] << " s" << endl;

  //Do and draw Akima interpolation
  	int_aspline = new ROOT::Math::Interpolator(x_points, y_points, ROOT::Math::Interpolation::kAKIMA);
	Int_t number_of_points = 4000;
	TGraph* general =  new TGraph(number_of_points);
	Double_t step = 12500/number_of_points;
	Double_t point = x_points[0];
	Int_t iter=0;
	while(point+iter*step<=x_points[nentries-1]){
			general->SetPoint(iter, point+iter*step, int_aspline->Eval(point+iter*step));
			iter++;
		}

//	TCanvas* canvasinter =  new TCanvas();
//	general->SetLineColor(kBlack);
//	general->SetLineWidth(2);
//	general->Draw("AL");

	cout << "Interpolation is ready" << endl;

	//-----------------------------------IRF reading--------------------------------------
	TString pathToReadIRF = "../../../ReadIRFs/";
//	TString pathToReadIRF = "./";
	const TString fileName3 = pathToReadIRF + "IRFsMrk421_2014flare_new.root";
	TFile* IRFfile = TFile::Open(fileName3, "READ");
	if (!IRFfile) {
		cout << "ERROR: file " << IRFfile << " not found..." << endl;
	}

	//Collection Area
	CollAr_Graph = dynamic_cast<TGraph*>(IRFfile->Get("CollArEtrue_Graph"));
	if (!CollAr_Graph) {
			cout << "ERROR: CollArEtrue_Graph object not found... " << endl;
	}

//	TCanvas* canvascollar = new TCanvas();
//	CollAr_Graph->Draw();


	//Bias and Resolution vs Etrue
	TGraph* ResolutionVsEtrue = dynamic_cast<TGraph*>(IRFfile->Get("ResolutionVsEtrue"));
	if (!ResolutionVsEtrue) {
			cout << "ERROR: ResolutionVsEtrue object not found... " << endl;
	}

	TGraph* BiasVsEtrue = dynamic_cast<TGraph*>(IRFfile->Get("BiasVsEtrue"));
	if (!BiasVsEtrue) {
			cout << "ERROR: BiasVsEtrue object not found... " << endl;
	}


//	TCanvas* canvasenerbiasresol = new TCanvas();
//	canvasenerbiasresol->Divide(2,1);
//	canvasenerbiasresol->cd(1);
//	canvasenerbiasresol->cd(1)->SetLogx();
//	ResolutionVsEtrue->Draw("ACP");
//	canvasenerbiasresol->cd(2)->SetLogx();
//	BiasVsEtrue->Draw("ACP");

	cout << "IRF are loaded and ready" << endl;





	//--------------------------Spectra folded with Coll Ar--------------------------------------

	folded_spectrum_flare = new TF1("folded_spectrum_flare", FoldingFlare, Emin, Emax, 0);
	folded_spectrum_flare->SetNpx(10000);

	cout << "Spectrum folded with CollAr" << endl;


	//----------------------------------EBL absorption------------------------------------------

	Double_t redshift = 0.03;
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

//	TCanvas* transcanvas = new TCanvas;
//	transcanvas->SetLogy();
//	transvse->GetXaxis()->SetTitle("Log10(E/GeV");
//	transvse->GetYaxis()->SetTitle("Transmitance");
//	transvse->SetLineColor(2);
//	transvse->SetLineWidth(2);
//	transvse->Draw("alp");

	//Apply the transmitance to the folded spectrum

	TF1* EBL_spectrum_flare = new TF1("EBL_spectrum_flare", EBLFlare, Emin, Emax, 0);
	EBL_spectrum_flare->SetNpx(10000);

//	TCanvas* foldedcanvas = new TCanvas();
//	foldedcanvas->SetLogx();
//	foldedcanvas->SetLogy();
//	folded_spectrum_flare->GetXaxis()->SetTitle("Energy (GeV)");
//	folded_spectrum_flare->SetTitle("Spectra");
//	folded_spectrum_flare->SetTitle("Spectra");
//	folded_spectrum_flare->SetLineColor(2);
//	folded_spectrum_flare->SetLineStyle(2);
//	folded_spectrum_flare->Draw();
//	flare_spectrum->SetLineColor(1);
//	flare_spectrum->SetLineStyle(1);
//	flare_spectrum->Draw("same");
//	EBL_spectrum_flare->SetLineColor(3);
//	EBL_spectrum_flare->SetLineStyle(2);
//	EBL_spectrum_flare->Draw("same");
//
//
//	TLegend* leg_spec_fold = new TLegend(0.1,0.7,0.4,0.9);
//	leg_spec_fold->AddEntry(flare_spectrum, "Flare (Log parabola)", "l");
//	leg_spec_fold->AddEntry(folded_spectrum_flare, "Flare+CollAr", "l");
//	leg_spec_fold->AddEntry(EBL_spectrum_flare, "Flare+CollAr+EBL", "l");
//	leg_spec_fold->Draw();

	cout << "EBL applied to the folded spectrum" << endl;


	//------------------------------------------QG inputs---------------------------------------
	Double_t z_source=0.031;
	Double_t H0=70.;
	Double_t PlanckScale=1.2*TMath::Power(10,19);
	Double_t PcTom=3.085*TMath::Power(10,19);
	Double_t QGfixedPart=z_source*PcTom/H0;
	std::cout << "alpha: " << alpha << std::endl;
	std::cout << "QG fixed part = " << QGfixedPart << std::endl;


	//------------------------Histograms for simulated data (time and energy)------------------------------

	double XX_step = 0.1;
	int n = ceil((log10(Emax)-log10(Emin))/XX_step);
	double *XX = new double[n+1];
	for(int i=0;i<n+1;i++)
	{
		XX[i]= pow(10,(log10(Emin)+i*XX_step));
	}

	//Read initial histogram to get the bin distribution
	TFile* HistoFile = TFile::Open("../../../ReadEvents/FinalEventLists/All_Events/AllEvents_ONhistograms_Mrk421_2014_Mireia.root", "READ");
//	TFile* HistoFile = TFile::Open("./LC_01-20TeV_ONhistograms_Mrk421_2014.root", "READ");
	if (!HistoFile) {
		cout << "ERROR: file " << HistoFile << " not found..."	<< endl;
	}
	TH1D* lightcurve = dynamic_cast<TH1D*>(HistoFile->Get("LC_insecs"));
	if (!lightcurve) {
		cout << "LC_insecs not found" << endl;
	}

	Int_t numberOfBins = lightcurve->GetNbinsX();
	cout << "Number of bins in initial histo: " << numberOfBins << endl;
	const Double_t* arrayBins = lightcurve->GetXaxis()->GetXbins()->GetArray();


	TH1D* sim_flare_time = new TH1D("sim_flare_time", "sim_flare_time", numberOfBins, arrayBins);
	TH1D* sim_flare_energy = new TH1D("sim_flare_energy", "sim_flare_energy", n, XX);
	sim_flare_energy->Sumw2();


	//-------------------------------------------Read Real Data (for comparison)------------------------------------------------
	ifstream in;
	in.open("../../../ReadEvents/FinalEventLists/All_Events/Mrk421_2014flare_ONEvents_Mireia_noHECut.txt");
	Double_t tp[20000],Ep[20000];
	Double_t v1,v2;
	int Npoints = 0;
	while(1)
	{
		in >> v1 >> v2;
		if (!in.good()) break;
		tp[Npoints] = v1;
		Ep[Npoints] = v2;
		Npoints++;
	}
	cout << "Npoints(Number of events): " << Npoints << endl;

	//Histos for data
	TH1D* data_time = new TH1D("Time distribution", "Time distribution", numberOfBins, arrayBins);
	TH1D* data_energy = new TH1D("Energy distribution", "Energy distribution", n, XX);
	data_energy->Sumw2();
	for(int i = 0; i<Npoints; i++){
		data_time->Fill(tp[i]);
		data_energy->Fill(Ep[i], 1/Ep[i]);
	}


	//----------------------------------------Simulation ------------------------------------------------

	//Histrogram for simulated data.
	Double_t QGDelay = 0.;
	  if(QGorder==1){ //Linear
		QGDelay= QGfixedPart*alpha/PlanckScale;
	  }
	  if(QGorder==2){ //Quadratic
		QGDelay= QGfixedPart*alpha/PlanckScale/1000.;
	  }
	std::cout << "Injected tau: " << QGDelay << std::endl;
	TH1D* timedel = new TH1D("timedel", "timedel", numberOfBins, arrayBins);
	TH1D* energyrec = new TH1D("energyrec", "energyrec", n, XX);
	energyrec->Sumw2();



	//ENERGY (Markus' version: Hit or miss with in LogE of the spectrum (CollAr, EBL later))
	Int_t index = 0;
	Int_t sim_events = 0;
	Int_t acceptEvents = 0;
	Int_t plotEvents = 0;
	Double_t a,b;
	Double_t resolution;

	TString form = flare_spectrum->GetExpFormula();
	form.ReplaceAll("[p0]*","");

	Double_t ymin = log(Emin/E0);
	Double_t ymax = log(Emax/E0);
	cout << "ymin: " << ymin << " ymax: " << ymax << endl;

	TString newform = (Form("%.0f*exp(x*(1.+%f-%f*(%f*x)))",E0,flare_spectrum->GetParameter(1),flare_spectrum->GetParameter(2),flare_spectrum->GetParameter(2)/ln10));
	TF1 *fnew = new TF1("fnew",newform,ymin,ymax);


	TRandom3 *ran = new TRandom3();
	ran->SetSeed(0);
	// initialize simulation:
	Double_t deltax = ymax-ymin;
	Double_t deltay = fnew->Eval(ymin)-fnew->Eval(ymax);
	Double_t pmin   = fnew->Eval(ymax);

	while(plotEvents<numberOfEvents){
		Double_t y = ymin + ran->Rndm()*deltax; //Coordinate in x
		Double_t p = pmin + deltay * ran->Rndm(); //Coordinate in y
		sim_events++;
		if(fnew->Eval(y) < p)continue;
		else{
			E[index] = get_E(y);
			acceptEvents++;

				if(fixed_resolution == true) resolution = 0.1;
				else resolution = ResolutionVsEtrue->Eval(E[index]);

			gRandom->Rannor(a,b);
			Es[index]=E[index]*(1.+b*resolution);
//			cout << "Smeared energy: " << Es[i] << endl;

				if(Es[index]>Emin && Es[index]<Emax){ //This selection is only for the plot!
					sim_flare_energy->Fill(E[index],1/E[index]);
					energyrec->Fill(Es[index], 1/Es[index]);
					plotEvents++;
					index++;
				}
			}
		}


	cout << "Energy distribution is simulated" << endl;
	cout << "Simulated events (Energy): " << sim_events << endl;
	cout << "Accepted events (Energy): " << acceptEvents << endl;
	cout << "Selected events for plot (Energy): " << plotEvents << endl;


	data_energy->SetStats(0);
	TCanvas* energysimcanvas = new TCanvas();
	energysimcanvas->SetLogx();
	energysimcanvas->SetLogy();
	data_energy->SetLineColor(1);
	data_energy->GetXaxis()->SetTitle("Energy(GeV)");
	data_energy->SetLineWidth(2);
	data_energy->Draw();
	sim_flare_energy->SetLineColor(3);
	sim_flare_energy->Draw("sameEP");
	energyrec->SetLineColor(2);
	energyrec->Draw("sameEP");

	TLegend* leg_energysim = new TLegend(0.1,0.7,0.4,0.9);
	leg_energysim->AddEntry(data_energy, "Data energy", "l");
	leg_energysim->AddEntry(sim_flare_energy, "Sim energy", "l");
	leg_energysim->AddEntry(energyrec, "Smeared energy", "l");
	leg_energysim->Draw();




	//TIME (Hit or miss with Interpolation and Gaussian)
	Double_t extra_time_edges = 4000.;
	Double_t x_min = tmin-extra_time_edges;
	Double_t x_max = tmax+extra_time_edges;
	Double_t y_min = 0.;
	Double_t y_max = 165;

//	cout << "x_min: " << x_min << " x_max: " << x_max << endl;
//	cout << "y_min: " << y_min << " y_max: " << y_max << endl;


	//Time simulation
	TRandom3* ranx = new TRandom3();
	ranx->SetSeed(0);
	TRandom3* rany = new TRandom3();
	rany->SetSeed(0);
	Double_t random_x, random_y;
	index = 0;
	acceptEvents = 0;
	plotEvents=0;
	while(plotEvents<numberOfEvents){
		random_x = ranx->Uniform(x_min,x_max);
		random_y = rany->Uniform(y_min,y_max);
		if(random_y <= lc_function(random_x)){
			acceptEvents++;
			t[index] = random_x;
			sim_flare_time->Fill(t[index],sim_flare_time->GetBinWidth(1)/sim_flare_time->GetBinWidth(sim_flare_time->FindBin(t[index])));
			if(QGorder == 1) td[index]=t[index]+QGDelay*E[index];
			if(QGorder == 2) td[index]=t[index]+QGDelay*E[index]*E[index];
			if(td[index]>tmin && td[index]<tmax){
				timedel->Fill(td[index],timedel->GetBinWidth(1)/timedel->GetBinWidth(timedel->FindBin(td[index])));
				plotEvents++;
				index++;
			}
		}
	}

	cout << "t[0]: " << t[0] << endl;
	cout << "td[0]: " << td[0] << endl;
	cout << "Time distribution is simulated" << endl;
	cout << "Acepted events (Time): " << acceptEvents << endl;
	cout << "Plot events (Time): " << plotEvents << endl;




	TCanvas* timesim = new TCanvas();
	data_time->SetStats(0);
	data_time->SetLineColor(1);
	data_time->SetLineWidth(2);
	data_time->Draw("EP");
	sim_flare_time->SetLineColor(2);
	sim_flare_time->Draw("sameEP");
	timedel->SetLineColor(4);
	timedel->Draw("sameEP");
	general->SetLineColor(3);
	general->Draw("same");


	TLegend* leg_timesim = new TLegend(0.1,0.7,0.4,0.9);
	leg_timesim->AddEntry(data_time, "Data time", "l");
	leg_timesim->AddEntry(sim_flare_time, "Sim time", "l");
	leg_timesim->AddEntry(timedel, "Delayed time", "l");
	leg_timesim->AddEntry(general, "Interpolation", "l");
	leg_timesim->Draw();

	//--------------------------------Select and order the events in increasing energy-----------------------

	cout << "Now I order the events in increasing Energy..."  << endl;
	Double_t temp_E;
	Double_t temp_t;
	for(Int_t i=0;i<numberOfEvents;i++)
	 {
		cout << i << endl;
		for(Int_t j=0;j<numberOfEvents-i;j++)
		{
			if(Es[j]>Es[j+1])
			{
				temp_E = Es[j];
				temp_t = td[j];
				Es[j] = Es[j+1];
				td[j] = td[j+1];
				Es[j+1]=temp_E;
				td[j+1]=temp_t;
			}
		}
	 }

	//--------------------------------Create a txt file with the events-----------------------

	//We have to select again after smearing!! IMPORTANT!!
	Int_t VHE_selected = 0;
	FILE *fout;
	fout = fopen(Form("./Mrk421_2014_flare_alpha%i_%i_qua_allfix_MarkusLP.txt", (int)alpha, simIndex),"wb");
//	fout = fopen(Form("./Mrk421_2014_flare_alpha%i_%i_fix.txt", (int)alpha, simIndex),"wb");
	FILE *fout2;
	fout2 = fopen(Form("./Mrk421_2014_flare_alpha%i_%i_qua_allfix_MarkusLP_VHE4.txt", (int)alpha, simIndex),"wb");
//	fout2 = fopen(Form("./Mrk421_2014_flare_alpha%i_%i_fix.txt", (int)alpha, simIndex),"wb");


	for(int i = 0; i<numberOfEvents; i++)
	{
		fprintf(fout,"%0.8lf    %0.8lf   \n", td[i], Es[i]);
		if(Es[i]>4000){
				fprintf(fout2,"%0.8lf    %0.8lf   \n", td[i], Es[i]);
				VHE_selected++;
			}
	}
	cout << "Written VHE events: " << VHE_selected << endl;





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


