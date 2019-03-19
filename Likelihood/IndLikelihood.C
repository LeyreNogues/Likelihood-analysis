/*
 * IndLikelihood.C
 *
 *  Created on: Jul 4, 2016
 *      Author: Leyre Nogués
 */


/* Parts of the code:
1. Read of events: Already cut and prepare list.
2. IRFs: Read CollAr and res/bias as function of energy (1 or 2 Gauss option).
3. EBL table: Read value from Dominguez table.
4. Spectrum: Intrinsic assumption from Fold Output.
5. Intrinsic time: Compute delay and extract intrinsic time of events.
6. Lightcurve at the source: Fit events and let some parameters free.
7. Normalization

Update: as time PDF, you can use an interpolation or KDE of the data events.
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
#include <TMinuit.h>
#include <TString.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>

using namespace std;



//---------------------------Likelihood options and log file----------------------------------
const Int_t enhance = 1;

Int_t QGorder= 2; //Linear or quadratic case
Bool_t ResolutionInc = true; //Include or not resolution
Bool_t include_CollAr = true;
Bool_t include_EBL = true;
Bool_t cout_PDF = false; //Print PDF values
Bool_t KDE_mode = false; //Use KDE or interpolation
Bool_t plot_results = true;
Bool_t two_gaus = true; //True for 2 Gaus resolution git, False for 1 Gaus fit
Bool_t PIC_job = false;
Bool_t poisson_fluctuations = false; //True to propagate Poisson fluctuations in the Interpolation.
Bool_t include_bkg = true;

ofstream log_file("Log_IndLike.txt", ios_base::out|ios_base::app);


//--------------------------------Functions for the PDF-----------------------------------

//IRF functions
TGraph* CollAr_Graph;
TH1D* Const1Erec;
TH1D* Const2Erec;
TH1D* Const0Erec;
TH1D* Mean1Erec;
TH1D* Mean2Erec;
TH1D* Mean0Erec;
TH1D* Sigma1Erec;
TH1D* Sigma2Erec;
TH1D* Sigma0Erec;

//Time and energy distributions
TF1* flare_spectrum;
TGraph* general;
TGraph* general_original;
TGraph* general_original_bkg;
TGraph *energy_dist_bkg;
TH1D *energy_dist_bkg_histo;
const Double_t E0 = 364.;
TH1D* energy_events;
TH1D* time_events;
Double_t Emax_bkg = 0.;;
Double_t integral_inter_signal = 0.;
Double_t integral_inter_bkg = 0.;
Double_t integral_spectrum_bkg;


//Vectors for the LC interpolation
ROOT::Math::Interpolator *int_aspline = NULL;
ROOT::Math::Interpolator *int_aspline_original = NULL;
ROOT::Math::Interpolator *int_aspline_background_original = NULL;

std::vector<double> x_points;
std::vector<double> y_points;
std::vector<double> y_mod_points;
std::vector<double> bin_widths;
Double_t inter_min = 0.;
Double_t inter_max = 0.;
Double_t new_inter_min = 0.;
Double_t new_inter_max = 0.;

std::vector<double> x_points_bkg;
std::vector<double> y_points_bkg;
std::vector<double> bin_widths_bkg;
Double_t inter_min_bkg = 0.;
Double_t inter_max_bkg = 0.;


//KDE function for LC
TF1* KDE_function;
TF1* KDE_function_bkg;

//Graph for EBL application
TGraph *transvse;

//Bins for integrals
Int_t nbinsE = 80; //100
Int_t nbinsT = 1000; //2000 //6000 for 9ON
Int_t nbinsEtrue = 80; //100
Int_t numberOfScapedEvents = 0;


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


//------------------------Flare properties (from reading code)---------------------------
Int_t numberOfEvents; //Real one
Double_t E[15000*enhance],t[15000*enhance];
Double_t tmin=0.;
Double_t tmax = hole[numberOfPoints-1];
Double_t Emin=120.;
Double_t Emax=25000.;




Double_t wobble_condition(Double_t t){
	Double_t wobble_condition = 1.;
	if(t>hole[1] && t<hole[2]) wobble_condition = 0.;
	if(t>hole[3] && t<hole[4]) wobble_condition = 0.;
	if(t>hole[5] && t<hole[6]) wobble_condition = 0.;
	if(t>hole[7] && t<hole[8]) wobble_condition = 0.;
	if(t>hole[9] && t<hole[10]) wobble_condition = 0.;
	if(t>hole[11] && t<hole[12]) wobble_condition = 0.;
	if(t>hole[13] && t<hole[14]) wobble_condition = 0.;
	if(t>hole[15] && t<hole[16]) wobble_condition = 0.;
	if(t>hole[17] && t<hole[18]) wobble_condition = 0.;
	if(t>hole[19] && t<hole[20]) wobble_condition = 0.;
	if(t>hole[21] && t<hole[22]) wobble_condition = 0.;
	if(t>hole[23] && t<hole[24]) wobble_condition = 0.;
	if(t>hole[25] && t<hole[26]) wobble_condition = 0.;

	return wobble_condition;
}



//------------------------------------------QG inputs-----------------------------------
Double_t z=0.031;
Double_t H0=70.;
//Double_t PlanckScale=1.2*TMath::Power(10,19);
Double_t PlanckScale=1.2;
Double_t QGpar=0;
//Double_t PcTom=3.085*TMath::Power(10,19);
Double_t PcTom=3.085;
Double_t QGfixedPart=z*PcTom/H0;
Double_t result_alpha;



//------------------------------Functions needed for the Likelihood--------------------------
Double_t CollArEvaluation(Double_t E){
	return CollAr_Graph->Eval(E);
}

Double_t BkgEnergyEvaluation(Double_t E){
	if(E<Emax_bkg) return energy_dist_bkg_histo->GetBinContent(energy_dist_bkg_histo->FindBin(E));
	else return 0;
}

Double_t lc_function(Double_t t){

	Double_t cons_value = 0.5363;
	Double_t weight_value = 1.353;
	Double_t peak_value = 4360.06;
	Double_t sigma_left = 2343.32;
	Double_t sigma_right = 1480.8;
	Double_t gaus_value = 0;

	//Compute values of the gaussian
	if(t<peak_value)
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_left,2));
	else
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_right,2));

	//Return Interpolation or Gaussian (For systematics study)
//    if(t<inter_min) return int_aspline_original->Eval(inter_min);
//    if(t>inter_max) return int_aspline_original->Eval(inter_max);


	//Return Interpolation or Gaussian
    if(t<inter_min || t>inter_max)
    	return gaus_value;
    else
    	if(poisson_fluctuations == 1) return int_aspline->Eval(t);
    	else if(poisson_fluctuations == 0) return int_aspline_original->Eval(t);
}

Double_t gaus_function(Double_t t){

	Double_t cons_value = 0.5363;
	Double_t weight_value = 1.353;
	Double_t peak_value = 4360.06;
	Double_t sigma_left = 2343.32;
	Double_t sigma_right = 1480.8;
	Double_t gaus_value = 0;

	//Compute value of the gaussian
	if(t<peak_value)
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_left,2));
	else
		gaus_value = cons_value + weight_value*TMath::Exp(-0.5*TMath::Power((t-peak_value)/sigma_right,2));

	return gaus_value;
}


Double_t GetDelay(Double_t energy, Double_t alpha)
{
	Double_t QGDelay = 0.;
	  if(QGorder==1){ //Linear
		QGDelay= QGfixedPart*alpha*energy/PlanckScale;
	  }

	  if(QGorder==2){ //Quadratic
		QGDelay= QGfixedPart*alpha*energy*energy/PlanckScale/1000.;
	  }

	return QGDelay;
}

Double_t Gaussian(Double_t x, Double_t mu, Double_t sigma)
{
	return (1./(sigma*TMath::Sqrt(2*TMath::Pi())))* TMath::Exp(-0.5*(x-mu)*(x-mu)/sigma/sigma);
}

Double_t GaussianE(Double_t* x, Double_t* par){

	return (1./par[0])*TMath::Exp(-0.5*(TMath::Log10(x[0])-par[1])*(TMath::Log10(x[0])-par[1])/par[2]/par[2])*(x[0]/TMath::Log(10.));
}


Double_t getPowerLawValue(Double_t x)
{
	Double_t index=-1.7;
	return TMath::Power(x,index);
}


Double_t Double_Gaus(Double_t* x, Double_t*  par){

	Double_t simgaus = par[0]*TMath::Exp(-0.5*TMath::Power((x[0]-par[1])/par[2],2));
	Double_t simgaus2 = par[3]*TMath::Exp(-0.5*TMath::Power((x[0]-par[4])/par[5],2));
	return simgaus + simgaus2;
}

Double_t Double_Gaussian(Double_t low_limit, Double_t up_limit, Double_t x, Double_t const1 , Double_t mean1, Double_t sigma1 ,Double_t const2, Double_t mean2, Double_t sigma2){

//	TF1* doble_gaus = new TF1("resol_gaus", Double_Gaus, low_limit, up_limit,6);
//	doble_gaus->SetParameters(const1, mean1, sigma1, const2, mean2, sigma2);
//	Double_t integral_value = doble_gaus->Integral(low_limit, up_limit);
//	cout << "Intergral in function: " << integral_value << endl;
//	doble_gaus->Delete();

	Double_t integral_loop = 0.;
	Double_t dEtrue = (up_limit-low_limit)/nbinsEtrue;
	for (Int_t j=0;j<nbinsEtrue;j++){
		Double_t Etrue = low_limit + (j+0.5)*dEtrue; //Half bin
//		cout << "Etrue: " << Etrue << endl;
		Double_t simgaus = const1*TMath::Exp(-0.5*TMath::Power((Etrue-mean1)/sigma1,2));
		Double_t simgaus2 = const2*TMath::Exp(-0.5*TMath::Power((Etrue-mean2)/sigma2,2));
		Double_t n=simgaus+simgaus2;
		integral_loop += n*dEtrue;
	}

//	cout << "Intergral in loop: " << integral_loop << endl;

	Double_t simgaus = const1*TMath::Exp(-0.5*TMath::Power((x-mean1)/sigma1,2));
	Double_t simgaus2 = const2*TMath::Exp(-0.5*TMath::Power((x-mean2)/sigma2,2));
	return (simgaus+simgaus2)/integral_loop;
}

//-------------------------------------PDF Building-------------------------------------


Double_t PDFevent(Double_t tr, Double_t Er, Double_t *par)
{
	Double_t QGParameter = par[0];

	Double_t spec_flare_value = 0.;
	Double_t lc_flare_value = 0.;

	Double_t PDF=0;
	spec_flare_value = flare_spectrum->Eval(Er);
	Double_t lc_point = tr-GetDelay(Er,QGParameter);
	lc_flare_value = KDE_function->Eval(lc_point);

	PDF = lc_flare_value*spec_flare_value; //Here only signal, bkg for resolution

	if(include_CollAr==true){
		Double_t collar_value = CollArEvaluation(Er);
		PDF=PDF*collar_value;
	}
	if(include_EBL==true){
		Double_t ebl_value = transvse->Eval(TMath::Log10(Er));
		PDF=PDF*ebl_value;
	}
	  return PDF;


}

Double_t PDFeventBkg(Double_t tr, Double_t Er)
{

	Double_t spec_flare_value = 0.;
	Double_t lc_flare_value = 0.;

	Double_t PDF_bkg=0;
	spec_flare_value = BkgEnergyEvaluation(Er);
	lc_flare_value = KDE_function_bkg->Eval(tr);

	PDF_bkg = lc_flare_value*spec_flare_value;

	return PDF_bkg;
}

Double_t PDFevent_Inter(Double_t tr, Double_t Er, Double_t *par)
{
	Double_t QGParameter = par[0];

	Double_t spec_flare_value = 0.;
	Double_t lc_flare_value = 0.;

	Double_t PDF=0;
	spec_flare_value = flare_spectrum->Eval(Er)*10; //x10 to have correct units (GeV, m^2)
	if(include_bkg==1) spec_flare_value = spec_flare_value*(tmax-tmin);
	Double_t lc_point = tr-GetDelay(Er,QGParameter);
	lc_flare_value = lc_function(lc_point);
	if(include_bkg==1) lc_flare_value = lc_flare_value/integral_inter_signal;

	PDF = lc_flare_value*spec_flare_value;

	if(include_CollAr==true){
		Double_t collar_value = CollArEvaluation(Er);
		PDF=PDF*collar_value;
	}
	if(include_EBL==true){
		Double_t ebl_value = transvse->Eval(TMath::Log10(Er));
		PDF=PDF*ebl_value;
	}

	  return PDF;
}

Double_t PDFeventBkg_Inter(Double_t tr, Double_t Er)
{

	Double_t spec_flare_value = 0.;
	Double_t lc_flare_value = 0.;

	Double_t PDF_bkg=0;
	spec_flare_value = BkgEnergyEvaluation(Er);
	lc_flare_value = int_aspline_background_original->Eval(tr);
	lc_flare_value = lc_flare_value/integral_inter_bkg;

	PDF_bkg = lc_flare_value*spec_flare_value;

	if(cout_PDF){
		cout << "Bkg spectrum" << spec_flare_value << endl;
		cout << "Bkg LC: " << lc_flare_value << endl;
	}

	return PDF_bkg;
}

Double_t PDFevent_Res(Double_t tr, Double_t Er, Double_t *par)
{

	//Parameters from the double Gaussian fit
	Double_t const1 = Const1Erec->GetBinContent(Const1Erec->FindBin(Er));
	Double_t const2 = Const2Erec->GetBinContent(Const2Erec->FindBin(Er));
	Double_t const0 = Const2Erec->GetBinContent(Const0Erec->FindBin(Er));
	Double_t mean1 = Mean1Erec->GetBinContent(Mean1Erec->FindBin(Er));
	Double_t mean2 = Mean2Erec->GetBinContent(Mean2Erec->FindBin(Er));
	Double_t mean0 = Mean2Erec->GetBinContent(Mean0Erec->FindBin(Er));
	Double_t sigma1 = Sigma1Erec->GetBinContent(Sigma1Erec->FindBin(Er));
	Double_t sigma2 = Sigma2Erec->GetBinContent(Sigma2Erec->FindBin(Er));
	Double_t sigma0 = Sigma0Erec->GetBinContent(Sigma0Erec->FindBin(Er));



	Double_t res = 0.;
	if (two_gaus == false) res = sigma0;
	else res = sigma2;
	Double_t Etrue_down = Er-5*Er*res; //Limits of the integral
	Double_t Etrue_up = Er+5*Er*res;
	if(Etrue_down<0) Etrue_down = 0;


	//PDF for background -> No need to integrate in Etrue
	Double_t PDF = 0.;

	//Integral over Etrue
	for (Int_t j=0;j<nbinsEtrue;j++){
		Double_t dEtrue = (Etrue_up-Etrue_down)/nbinsEtrue;
		Double_t Etrue = Etrue_down + (j+0.5)*dEtrue; //Half bin
		Double_t n = 0.;
		if (two_gaus == false) n = PDFevent(tr,Etrue, par)*Gaussian(Etrue,(1-mean0)*Er,sigma0*Er);
		else n = PDFevent(tr,Etrue, par)*Double_Gaussian(Etrue_down, Etrue_up, Etrue, const1,(1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er);
		PDF += n*dEtrue;
	}

	return PDF;
}



Double_t PDFeventRes_Inter(Double_t tr, Double_t Er, Double_t *par)
{

	//Parameters from the double Gaussian fit
	Double_t const1 = Const1Erec->GetBinContent(Const1Erec->FindBin(Er));
	Double_t const2 = Const2Erec->GetBinContent(Const2Erec->FindBin(Er));
	Double_t const0 = Const2Erec->GetBinContent(Const0Erec->FindBin(Er));
	Double_t mean1 = Mean1Erec->GetBinContent(Mean1Erec->FindBin(Er));
	Double_t mean2 = Mean2Erec->GetBinContent(Mean2Erec->FindBin(Er));
	Double_t mean0 = Mean2Erec->GetBinContent(Mean0Erec->FindBin(Er));
	Double_t sigma1 = Sigma1Erec->GetBinContent(Sigma1Erec->FindBin(Er));
	Double_t sigma2 = Sigma2Erec->GetBinContent(Sigma2Erec->FindBin(Er));
	Double_t sigma0 = Sigma0Erec->GetBinContent(Sigma0Erec->FindBin(Er));


	Double_t res = 0.;
	if (two_gaus == false) res = sigma0;
	else res = sigma2;
	Double_t Etrue_down = Er-5*Er*res; //Limits of the integral
	Double_t Etrue_up = Er+5*Er*res;
	if(Etrue_down<0) Etrue_down = 0;


	//PDF for background -> No need to integrate in Etrue
	Double_t PDF = 0.;

	//Integral over Etrue
	for (Int_t j=0;j<nbinsEtrue;j++){
		Double_t dEtrue = (Etrue_up-Etrue_down)/nbinsEtrue;
		Double_t Etrue = Etrue_down + (j+0.5)*dEtrue; //Half bin
		Double_t n = 0.;
		if (two_gaus == false) n = PDFevent_Inter(tr,Etrue, par)*Gaussian(Etrue,(1-mean0)*Er,sigma0*Er);
		else n = PDFevent_Inter(tr,Etrue, par)*Double_Gaussian(Etrue_down, Etrue_up, Etrue, const1,(1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er);
		PDF += n*dEtrue;
	}

	return PDF;
}


Double_t PDFtest(Double_t tr, Double_t Er)
{
	Double_t QGParameter=QGpar;

	Double_t spec_flare_value = 0.;
	Double_t lc_flare_value = 0.;

	Double_t PDF=0;
	spec_flare_value = flare_spectrum->Eval(Er);

	Double_t lc_point = tr-GetDelay(Er,QGParameter);

	lc_flare_value = KDE_function->Eval(lc_point);
	PDF = lc_flare_value*spec_flare_value;

	if(include_CollAr==true){
		Double_t collar_value = CollArEvaluation(Er);
		PDF=PDF*collar_value;
		if(cout_PDF == true) cout << "CollAr value: " << collar_value << endl;
	}

	if(include_EBL==true){
		Double_t ebl_value = transvse->Eval(TMath::Log10(Er));
		PDF=PDF*ebl_value;
		if(cout_PDF == true) cout << "EBL value: " << ebl_value << endl;

	}

	if(cout_PDF == true){
		cout << "tr Er " << tr << " " << Er << endl;
		cout << "QGpar " << QGParameter << endl;
		cout << "GetDelay: " << GetDelay(Er,QGParameter) << endl;
		cout << "spec_flare_value: " << spec_flare_value << endl;
		cout << "lc_point: " << lc_point << endl;
		cout << "lc_flare_value: " << lc_flare_value << endl;
		cout << "PDF(aft CA and EBL) " << PDF << endl;

	}

	return PDF;
}


Double_t PDFtest_Inter(Double_t tr, Double_t Er){

	Double_t QGParameter=QGpar;

	Double_t spec_flare_value = 0.;
	Double_t lc_flare_value = 0.;

	Double_t PDF=0;
	spec_flare_value = flare_spectrum->Eval(Er)*10;
	if(include_bkg==1) spec_flare_value = spec_flare_value*(tmax-tmin);

	Double_t lc_point = tr-GetDelay(Er,QGParameter);
	lc_flare_value = lc_function(lc_point);
	if(include_bkg==1) lc_flare_value = lc_flare_value/integral_inter_signal;

	PDF = lc_flare_value*spec_flare_value;

	if(include_CollAr==true){
		Double_t collar_value = CollArEvaluation(Er);
		PDF=PDF*collar_value;
		if(cout_PDF == true) cout << "CollAr value: " << collar_value << endl;

	}
	if(include_EBL==true){
		Double_t ebl_value = transvse->Eval(TMath::Log10(Er));
		PDF=PDF*ebl_value;
		if(cout_PDF == true) cout << "EBL value: " << ebl_value << endl;
	}

	if(cout_PDF == true){
		cout << "tr Er " << tr << " " << Er << endl;
		cout << "QGpar " << QGParameter << endl;
		cout << "GetDelay: " << GetDelay(Er,QGParameter) << endl;
		cout << "spec_flare_value: " << spec_flare_value << endl;
		cout << "lc_point: " << lc_point << endl;
		cout << "lc_flare_value: " << lc_flare_value << endl;
		cout << "PDF(aft CA and EBL) " << PDF << endl;

	}

	return PDF;

}


Double_t PDFtestRes(Double_t tr, Double_t Er)
{
	//Parameters from the double Gaussian fit
	Double_t const1 = Const1Erec->GetBinContent(Const1Erec->FindBin(Er));
	Double_t const2 = Const2Erec->GetBinContent(Const2Erec->FindBin(Er));
	Double_t const0 = Const0Erec->GetBinContent(Const0Erec->FindBin(Er));
	Double_t mean1 = Mean1Erec->GetBinContent(Mean1Erec->FindBin(Er));
	Double_t mean2 = Mean2Erec->GetBinContent(Mean2Erec->FindBin(Er));
	Double_t mean0 = Mean0Erec->GetBinContent(Mean0Erec->FindBin(Er));
	Double_t sigma1 = Sigma1Erec->GetBinContent(Sigma1Erec->FindBin(Er));
	Double_t sigma2 = Sigma2Erec->GetBinContent(Sigma2Erec->FindBin(Er));
	Double_t sigma0 = Sigma0Erec->GetBinContent(Sigma0Erec->FindBin(Er));


	Double_t res = 0.;
	if (two_gaus == false) res = sigma0;
	else res = sigma2;
	Double_t Etrue_down = Er-5*Er*res; //Limits of the integral
	Double_t Etrue_up = Er+5*Er*res;
	if(Etrue_down<0) Etrue_down = 0;

	if(cout_PDF){
		cout << "Const0 Const1 Const2 " << const0 << " " << const1 << " " << const2 << endl;
		cout << "mean0 mean1 mean2 " << mean0 << " " << mean1 << " " << mean2 << endl;
		cout << "sigma0 sigma1 sigma2 " << sigma0 << " " << sigma1 << " " << sigma2 << endl;

		cout << "tr Er " << tr << " " << Er << endl;
		cout << "Resolution " << res << endl;
		cout << "Etrue_down Etrue_up " << Etrue_down << " " << Etrue_up << endl;
		if (two_gaus == false) cout << "I use a single Gaus" << endl;
		if (two_gaus == true) cout << "I use a doble Gaus" << endl;
	}

//		Gaussian function (Only for drawing purposes)
//		TF1* resol_gaus = new TF1("resol_gaus", "gaus(0)", Etrue_down, Etrue_up);
//		Double_t norm_factor = 1./(sigma0*Er*TMath::Sqrt(2*TMath::Pi()));
//		resol_gaus->SetParameters(norm_factor, (1-mean0)*Er, sigma0*Er);

	//	TCanvas* resolcanvas = new TCanvas();
	//	resol_gaus->Draw();

		//Double Gaussian function (Only for drawing purposes)
//		TF1* resol_doble_gaus = new TF1("resol_gaus", Double_Gaus, Etrue_down, Etrue_up,6);
//		resol_doble_gaus->SetParameters(const1,(1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er);

		//	TCanvas* resolcanvasdouble = new TCanvas();
		//	resol_doble_gaus->Draw();

		Double_t PDF = 0.;

		//Integral over Etrue
		for (Int_t j=0;j<nbinsEtrue;j++){
			Double_t dEtrue = (Etrue_up-Etrue_down)/nbinsEtrue;
			Double_t Etrue = Etrue_down + (j+0.5)*dEtrue; //Half bin
			Double_t n= 0.;
			if (two_gaus == false)  n = PDFtest(tr,Etrue)*Gaussian(Etrue,(1-mean0)*Er,sigma0*Er);
			else  n = PDFtest(tr,Etrue)*Double_Gaussian(Etrue_down, Etrue_up, Etrue, const1, (1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er);
			PDF += n*dEtrue;
		}

//		resol_gaus->Delete();
//		resol_doble_gaus->Delete();

		if(cout_PDF){
			cout << "PDF for signal " << PDF << endl;
		}

		return PDF;

}

Double_t PDFtestRes_Inter(Double_t tr, Double_t Er)
{
	//Parameters from the double Gaussian fit
	Double_t const1 = Const1Erec->GetBinContent(Const1Erec->FindBin(Er));
	Double_t const2 = Const2Erec->GetBinContent(Const2Erec->FindBin(Er));
	Double_t const0 = Const0Erec->GetBinContent(Const0Erec->FindBin(Er));
	Double_t mean1 = Mean1Erec->GetBinContent(Mean1Erec->FindBin(Er));
	Double_t mean2 = Mean2Erec->GetBinContent(Mean2Erec->FindBin(Er));
	Double_t mean0 = Mean0Erec->GetBinContent(Mean0Erec->FindBin(Er));
	Double_t sigma1 = Sigma1Erec->GetBinContent(Sigma1Erec->FindBin(Er));
	Double_t sigma2 = Sigma2Erec->GetBinContent(Sigma2Erec->FindBin(Er));
	Double_t sigma0 = Sigma0Erec->GetBinContent(Sigma0Erec->FindBin(Er));


	Double_t res = 0.;
	if (two_gaus == false) res = sigma0;
	else res = sigma2;
	Double_t Etrue_down = Er-5*Er*res; //Limits of the integral
	Double_t Etrue_up = Er+5*Er*res;
	if(Etrue_down<0) Etrue_down = 0;

	if(cout_PDF){
		cout << "Const0 Const1 Const2 " << const0 << " " << const1 << " " << const2 << endl;
		cout << "mean0 mean1 mean2 " << mean0 << " " << mean1 << " " << mean2 << endl;
		cout << "sigma0 sigma1 sigma2 " << sigma0 << " " << sigma1 << " " << sigma2 << endl;

		cout << "tr Er " << tr << " " << Er << endl;
		cout << "Resolution " << res << endl;
		cout << "Etrue_down Etrue_up " << Etrue_down << " " << Etrue_up << endl;
		if (two_gaus == false) cout << "I use a single Gaus" << endl;
		if (two_gaus == true) cout << "I use a doble Gaus" << endl;
	}


//	//Gaussian function (Only for drawing purposes)
//	TF1* resol_gaus = new TF1("resol_gaus", "gaus(0)", Etrue_down, Etrue_up);
//	resol_gaus->SetParameters(1./(sigma0*Er*TMath::Sqrt(2*TMath::Pi())), (1-mean0)*Er, sigma0*Er);
//
//	if(cout_PDF && two_gaus == false){
//		TCanvas* resolcanvas = new TCanvas();
//		resol_gaus->Draw();
//	}
//
//	//Double Gaussian function (Only for drawing purposes)
//	TF1* resol_doble_gaus = new TF1("resol_gaus", Double_Gaus, Etrue_down, Etrue_up, 6);
//	resol_doble_gaus->SetParameters(const1,(1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er);
////	Double_t integral = resol_doble_gaus->Integral(Etrue_down, Etrue_up);
////	cout << "double integral: " << integral << endl;
//
//
//	if(cout_PDF && two_gaus == true){
//		TCanvas* resolcanvas2 = new TCanvas();
//		resol_doble_gaus->Draw();
//	}

	//PDF for background -> No need to integrate in Etrue
	Double_t PDF = 0.;

	//Integral over Etrue
	Double_t dEtrue = (Etrue_up-Etrue_down)/nbinsEtrue;
	for (Int_t j=0;j<nbinsEtrue;j++){
		Double_t Etrue = Etrue_down + (j+0.5)*dEtrue; //Half bin
//		cout << "Etrue: " << Etrue << endl;
		Double_t n=0.;
		if (two_gaus == false){
			n = PDFtest_Inter(tr,Etrue)*Gaussian(Etrue,(1-mean0)*Er,sigma0*Er);
//			cout << "Gaus value: " << Gaussian(Etrue,(1-mean0)*Er,sigma0*Er) << endl;
		}
		else{
			n = PDFtest_Inter(tr,Etrue)*Double_Gaussian(Etrue_down, Etrue_up, Etrue, const1,(1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er);
//			cout << "Doble Gaus value: " << Double_Gaussian(Etrue_down, Etrue_up, Etrue, const1,(1-mean1)*Er, sigma1*Er, const2, (1-mean2)*Er, sigma2*Er) << endl;
		}
//		cout << "n: " << n << endl;

		PDF += n*dEtrue;
	}

	return PDF;


}


Double_t Likelihoodtest(){

	cout << "Inside Likelihoodtest" << endl;
	// calculate normalization
	  Double_t Normalization = 0.;
	  const Double_t binsize  = TMath::Log(Emax/Emin)/nbinsE;

	  for (Int_t i=0;i<nbinsE;i++){
		  for (Int_t j=0;j<nbinsT;j++){
			  Double_t dE     = Emin * (TMath::Exp(binsize * (i+1)) - TMath::Exp(binsize*i));
			  Double_t dT     = (tmax-tmin)/nbinsT;
			  Double_t energy = Emin * TMath::Exp(binsize * (i+0.5));
			  Double_t time   = tmin + (j+0.5)*dT;
			  Double_t n=0.;
			  Double_t n_bkg=0.;
			  if(ResolutionInc == false){
				 n=PDFtest(time,energy)*wobble_condition(time);
			  }
			  if(ResolutionInc == true){
				 n= PDFtestRes(time,energy)*wobble_condition(time);
			  }
			  Normalization += n * dE * dT;
		  }
	  }

	  Normalization=TMath::Log(Normalization);

//	  cout << "Signal Normalization for: alpha=" << QGpar << " is: " << Normalization << endl;


	  Double_t NormalizationBkg = 0.;
	  if(include_bkg==1){
		  for (Int_t j=0;j<nbinsT;j++){
				  Double_t dT     = (tmax-tmin)/nbinsT;
				  Double_t time   = tmin + (j+0.5)*dT;
		  //		  cout << "Time " << time << endl;
				  Double_t n_bkg =0.;
				  n_bkg = KDE_function_bkg->Eval(time)*wobble_condition(time);
				  NormalizationBkg += n_bkg * dT;
			  }

//		  cout << "Bkg Normalization for: alpha=" << QGpar << " is: " << NormalizationBkg << endl;

	  }



	  Double_t suma=0.;
	  Double_t temp=0.;
	  for(Int_t i=0; i<numberOfEvents; i++){
		 	  if(ResolutionInc == false){
					if(include_bkg == 1) temp= TMath::Log((PDFtest(t[i],E[i])+PDFeventBkg(t[i],E[i]))/(Normalization+NormalizationBkg));
					if(include_bkg == 0)temp=TMath::Log(PDFtest(t[i],E[i]))- Normalization;
		 	  }
		 	  if(ResolutionInc == true){
		 		  if(include_bkg == 1)temp= TMath::Log((PDFtestRes(t[i],E[i])+PDFeventBkg(t[i],E[i]))/(Normalization+NormalizationBkg));
		 		  if(include_bkg == 0)temp= TMath::Log(PDFtestRes(t[i],E[i]))-Normalization;
		 	  }
		 	 suma+=(-2)*temp;
	  }

	  return suma;
}

Double_t Likelihoodtest_Inter(){

  // calculate normalization
  Double_t Normalization = 0.;
  Double_t NormalizationBkg = 0.;

  const Double_t binsize  = TMath::Log(Emax/Emin)/nbinsE;

  for (Int_t i=0;i<nbinsE;i++){
	  for (Int_t j=0;j<nbinsT;j++){
		  Double_t dE = Emin * (TMath::Exp(binsize * (i+1)) - TMath::Exp(binsize*i));
		  Double_t dT = (tmax-tmin)/nbinsT;
		  Double_t energy = Emin * TMath::Exp(binsize * (i+0.5));
		  Double_t time = tmin + (j+0.5)*dT;
//		  cout << "Energy: " << energy << endl;
//		  cout << "Time " << time << endl;
		  Double_t n=0.;
		  if(ResolutionInc == false){
			 n=PDFtest_Inter(time,energy)*wobble_condition(time);
		  }
		  if(ResolutionInc == true){
			 n= PDFtestRes_Inter(time,energy)*wobble_condition(time);
		  }
//		  cout << "n: " << n << endl;
		  Normalization += n * dE * dT;
	  }
  }

  if(include_bkg==0) Normalization = TMath::Log(Normalization);

  if(include_bkg==1){
	  for (Int_t j=0;j<nbinsT;j++){
			  Double_t dT     = (tmax-tmin)/nbinsT;
			  Double_t time   = tmin + (j+0.5)*dT;
	  //		  cout << "Time " << time << endl;
			  Double_t n_bkg =0.;
			  n_bkg = int_aspline_background_original->Eval(time)*wobble_condition(time)/integral_inter_bkg;
			  NormalizationBkg += n_bkg * dT;
		  }


//	  cout << "Bkg Normalization for: alpha=" << QGpar << " is: " << NormalizationBkg << endl;
	  NormalizationBkg = NormalizationBkg*integral_spectrum_bkg;
//	  cout << "Bkg Normalization corrected = " << NormalizationBkg << endl;

  }

  cout << "Signal Normalization for: alpha=" << QGpar << " is: " << Normalization << endl;


  Double_t suma=0.;
  Double_t temp=0.;
  for(Int_t i=0; i<numberOfEvents; i++){
	 	  if(ResolutionInc == false){
				if(include_bkg == 1) temp= TMath::Log((PDFtest_Inter(t[i],E[i])+PDFeventBkg_Inter(t[i],E[i]))/(Normalization+NormalizationBkg));
				if(include_bkg == 0)temp=TMath::Log(PDFtest_Inter(t[i],E[i]))- Normalization;
	 	  }
	 	  if(ResolutionInc == true){
	 		  if(include_bkg == 1)temp= TMath::Log((PDFtestRes_Inter(t[i],E[i])+PDFeventBkg_Inter(t[i],E[i]))/(Normalization+NormalizationBkg));
	 		  if(include_bkg == 0)temp= TMath::Log(PDFtestRes_Inter(t[i],E[i]))-Normalization;
	 	  }
	 	  if(temp >10000) cout << "nan for " << i << endl;

//		 cout << "For event: " << i << " temp is: " << temp << endl;
	 	 suma+=(-2)*temp;
  }

  return suma;
}

void Normalization_Compute(){

  // calculate time normalization

	integral_inter_signal=0.;
	integral_inter_bkg=0.;
	for (Int_t j=0;j<nbinsT;j++){
		Double_t dT     = (tmax-tmin)/nbinsT;
		Double_t time   = tmin + (j+0.5)*dT;
//		  cout << "Time " << time << endl;
		Double_t n=0.;
		Double_t n_bkg =0.;
		Double_t integ_time_signal = 0.;
		Double_t inter_time_bkg = 0.;
		if(poisson_fluctuations == false){
			integ_time_signal = int_aspline_original->Eval(time)*wobble_condition(time);
			if(include_bkg == 1)inter_time_bkg = int_aspline_background_original->Eval(time)*wobble_condition(time);
		}
		if(poisson_fluctuations == true){
			integ_time_signal = int_aspline->Eval(time)*wobble_condition(time);
			if(include_bkg == 1)inter_time_bkg = int_aspline_background_original->Eval(time)*wobble_condition(time);
		}

		integral_inter_signal += integ_time_signal * dT;
		if(include_bkg == 1)integral_inter_bkg += inter_time_bkg * dT;
	}

	return;

}

Double_t Likelihoodtest_Inter_Energy(Double_t min_Energy, Double_t max_Energy){

  // calculate normalization
  Double_t Normalization = 0.;
  Int_t events = 0;
  const Double_t binsize  = TMath::Log(Emax/Emin)/nbinsE;

  for (Int_t i=0;i<nbinsE;i++){
	  for (Int_t j=0;j<nbinsT;j++){
		  Double_t dE     = Emin * (TMath::Exp(binsize * (i+1)) - TMath::Exp(binsize*i));
		  Double_t dT     = (tmax-tmin)/nbinsT;
		  Double_t energy = Emin * TMath::Exp(binsize * (i+0.5));
		  Double_t time   = tmin + (j+0.5)*dT;
		  Double_t n=0.;
		  if(ResolutionInc == false){
			 n=PDFtest_Inter(time,energy)*wobble_condition(time);
		  }
		  if(ResolutionInc == true){
			 n= PDFtestRes_Inter(time,energy)*wobble_condition(time);
		  }
		  Normalization += n * dE * dT;
	  }
  }

  Normalization = TMath::Log(Normalization);

  Double_t suma=0.;
  Double_t temp=0.;
  for(Int_t i=0; i<numberOfEvents; i++){
	  if(E[i]>min_Energy && E[i]<max_Energy){
		  events++;
		  cout << "i: " << i << endl;
		  cout << "PDF part (Log) "<<TMath::Log(PDFtest_Inter(t[i],E[i])) << endl;
		  cout << "Normalization (Log) "<< Normalization<< endl;
	 	  if(ResolutionInc == false){
	 		 temp= TMath::Log(PDFtest_Inter(t[i],E[i]))-Normalization;
	 	  }
	 	  if(ResolutionInc == true){
	 		 temp= TMath::Log(PDFtestRes_Inter(t[i],E[i]))-Normalization;
	 	  }
	  suma+=(-2)*temp;
	  }
  }
  cout << "events= " << events << endl;

  return suma;
}

void fcn_KDE(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t logL = 0.;
  Double_t temp = 0.;

//  cout << "Inside fcn_KDE" << endl;
  // calculate normalization
  Double_t Normalization = 0.;

  numberOfScapedEvents = 0;
  const Double_t binsize  = TMath::Log(Emax/Emin)/nbinsE;

  for (Int_t i=0;i<nbinsE;i++){
	  for (Int_t j=0;j<nbinsT;j++){
		  Double_t dE     = Emin * (TMath::Exp(binsize * (i+1)) - TMath::Exp(binsize*i));
		  Double_t dT     = (tmax-tmin)/nbinsT;
		  Double_t energy = Emin * TMath::Exp(binsize * (i+0.5));
		  Double_t time   = tmin + (j+0.5)*dT;
		  Double_t n=0.;
		  Double_t n_bkg=0.;
		  if(ResolutionInc == false){
			 n=PDFevent(time,energy, par)*wobble_condition(time);
		  }
		  if(ResolutionInc == true){
			 n= PDFevent_Res(time,energy, par)*wobble_condition(time);
		  }
		  Normalization += n * dE * dT;
	  }
  }

  Normalization = TMath::Log(Normalization);

  Double_t NormalizationBkg = 0.;
    if(include_bkg==1){
  	  for (Int_t j=0;j<nbinsT;j++){
  			  Double_t dT     = (tmax-tmin)/nbinsT;
  			  Double_t time   = tmin + (j+0.5)*dT;
  	  //		  cout << "Time " << time << endl;
  			  Double_t n_bkg =0.;
  			  n_bkg = int_aspline_background_original->Eval(time)*wobble_condition(time);
  			  NormalizationBkg += n_bkg * dT;
  		  }
    }


    for (Int_t i=0; i<numberOfEvents; i++)
    {
  	  if(ResolutionInc == false){
  		  if(include_bkg==1)temp= TMath::Log((PDFevent(t[i],E[i], par)+PDFeventBkg(t[i],E[i]))/(Normalization+NormalizationBkg));
  		  if(include_bkg==0)temp= TMath::Log(PDFevent(t[i],E[i],par))-Normalization;

  	  }
  	  if(ResolutionInc == true){
  		  if(include_bkg==1)temp= TMath::Log((PDFevent_Res(t[i],E[i], par)+PDFeventBkg(t[i],E[i]))/(Normalization+NormalizationBkg));
  		  if(include_bkg==0)temp= TMath::Log(PDFevent_Res(t[i],E[i],par))-Normalization;
  	  }
     logL+=(-2)*temp; //Final probability function (addition of all events).
    }
    f = logL;

  cout <<  par[0]  << " " << f << endl;
}

void fcn_Inter(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

  // calculate normalization
  Double_t Normalization = 0.;
  Double_t NormalizationBkg = 0.;

  const Double_t binsize  = TMath::Log(Emax/Emin)/nbinsE;

  for (Int_t i=0;i<nbinsE;i++){
	  for (Int_t j=0;j<nbinsT;j++){
		  Double_t dE     = Emin * (TMath::Exp(binsize * (i+1)) - TMath::Exp(binsize*i));
		  Double_t dT     = (tmax-tmin)/nbinsT;
		  Double_t energy = Emin * TMath::Exp(binsize * (i+0.5));
		  Double_t time   = tmin + (j+0.5)*dT;
		  Double_t n=0.;
		  if(ResolutionInc == false){
			 n=PDFevent_Inter(time,energy, par)*wobble_condition(time);
		  }
		  if(ResolutionInc == true){
			 n=PDFeventRes_Inter(time,energy, par)*wobble_condition(time);
		  }
		  Normalization += n * dE * dT;
	  }
  }

  if(include_bkg==0) Normalization = TMath::Log(Normalization);

  if(include_bkg==1){
	  for (Int_t j=0;j<nbinsT;j++){
			  Double_t dT     = (tmax-tmin)/nbinsT;
			  Double_t time   = tmin + (j+0.5)*dT;
	  //		  cout << "Time " << time << endl;
			  Double_t n_bkg =0.;
			  n_bkg = int_aspline_background_original->Eval(time)*wobble_condition(time)/integral_inter_bkg;
			  NormalizationBkg += n_bkg * dT;
		  }

//	  cout << "Bkg Normalization for: alpha=" << QGpar << " is: " << NormalizationBkg << endl;
	  NormalizationBkg = NormalizationBkg*integral_spectrum_bkg;
//	  cout << "Bkg Normalization corrected = " << NormalizationBkg << endl;

  }

//  cout << "Signal Normalization for: alpha=" << par[0] << " is: " << Normalization << endl;


  Double_t logL = 0.;
  Double_t temp = 0.;

  for (Int_t i=0; i<numberOfEvents; i++)
  {
	  if(ResolutionInc == false){
		  if(include_bkg==1)temp= TMath::Log((PDFevent_Inter(t[i],E[i], par)+PDFeventBkg_Inter(t[i],E[i]))/(Normalization+NormalizationBkg));
		  if(include_bkg==0)temp= TMath::Log(PDFevent_Inter(t[i],E[i],par))-Normalization;

	  }
	  if(ResolutionInc == true){
		  if(include_bkg==1)temp= TMath::Log((PDFeventRes_Inter(t[i],E[i], par)+PDFeventBkg_Inter(t[i],E[i]))/(Normalization+NormalizationBkg));
		  if(include_bkg==0)temp= TMath::Log(PDFeventRes_Inter(t[i],E[i],par))-Normalization;
	  }
   logL+=(-2)*temp; //Final probability function (addition of all events).
  }
  f = logL;

  cout <<  par[0]  << " " << f << endl;
}

Int_t IndLikelihood(TString fname)
//Int_t IndLikelihood()
{
	//-----------------------------Read time and energy of events----------------------------

	//Txt version (For real data)
	ifstream in;
	in.open(fname);
	Double_t v1,v2;
	Int_t Npoints = 0;
	while(1)
	{
		if (in.eof()) break;
		in >> v1 >> v2;

		//Cuts in the events
		if(v1<0 || v1>4000500)continue;
		if(v2<0 || v2>6000000)continue;

		t[Npoints] = v1;
//		if(t[Npoints]>13000.1) cout << "TOO BIG TIME!  t = " << t[Npoints] << endl;
		E[Npoints] = v2;
		if(E[Npoints]<0.) cout << "NEGATIVE ENERGY EVENT!  E = " << E[Npoints] << endl;
		Npoints++;
	}

	cout << "Number of events " << Npoints << endl;
	numberOfEvents = Npoints;


//	//Tree version (For simulated data)
//	cout << "I will read now" << endl;
//	TString InputfileName = fname;
//	cout << "File name: " << InputfileName << endl;
//
//	 TFile* Eventfile = TFile::Open(InputfileName, "READ");
//	 if (!Eventfile) {
//		 log_file << "ERROR: Eventfile " << Eventfile << " not found..." << endl;
//	 }
//
//	  cout << "I found the file" << endl;
//
//	  //Read the Tree
//	  TTree *tevents = (TTree*)Eventfile->Get("Event_tree.root");
//	  if (!tevents) {
//	  		  log_file << "ERROR: tree file " << tevents << " not found..." << endl;
//	  }
//	  Double_t Es, td;
//	  tevents->SetBranchAddress("E",&Es);
//	  tevents->SetBranchAddress("td",&td);
//
//	  energy_events = new TH1D("energy_events", "energy_events", 100, Emin, Emax);
//	  time_events = new TH1D("time_events", "time_events", 100, tmin, tmax);
//
//	  cout <<"Emin: " << Emin << " Emax: " << Emax << endl;
//
//
//	  Int_t VHE_events = 0;
////	  UInt_t trialEvents= 12994*enhance; //Choose here the number of events to read
//	  //--Events in tree
//	  Long64_t nevents = tevents->GetEntries();
//	  cout << "Number of events in the tree: " << nevents << endl;
//	  UInt_t Npoints = 0;
//	  for(Long64_t i=0;i<nevents;i++){
//		  tevents->GetEntry(i);
//		  E[Npoints] = Es;
//		  if(E[Npoints]<120.)E[Npoints]=120.;
//		  energy_events->Fill(Es);
//		  if(Es>4000)VHE_events++;
//		  t[Npoints] = td;
//		  time_events->Fill(td);
//		  Npoints++;
////		  if(Npoints == trialEvents)break;
//	  }
//	  numberOfEvents = Npoints;
//	  cout << "Number of events " << numberOfEvents << endl;
//	  cout << "Number of VHE events " << VHE_events << endl;
//
//
//
//	  if(plot_results){
//		  TCanvas* events_canvas =  new TCanvas();
//		  events_canvas->Divide(2,1);
//		  events_canvas->cd(1);
//		  events_canvas->cd(1)->SetLogx();
//		  events_canvas->cd(1)->SetLogy();
//		  energy_events->Draw();
//		  events_canvas->cd(2);
//		  time_events->Draw();
//	  }


  //----------------------------------Spectra inputs-------------------------------------------

	  TString pathToReadSpec;
	  if(PIC_job == false) pathToReadSpec = "../Simulation/Mrk421/2014/Fold_outputs/";
	  else pathToReadSpec = "./";
	  const TString fileName = pathToReadSpec + "Output_fold_Zd50_LC120_EPWL.root";
	  TFile* Specfile = TFile::Open(fileName, "READ");
	  if (!Specfile) {
		  log_file << "ERROR: file " << Specfile << " not found..." << endl;
	  }

	  flare_spectrum = dynamic_cast<TF1*>(Specfile->Get("SpectralModel"));
	  if (!flare_spectrum) {
		  log_file << "ERROR: SpectralModel object not found... " << endl;
	  }

//	Double_t bin_value = energy_events->GetBinContent(energy_events->FindBin(E0));
//	flare_spectrum->SetParameter(0,bin_value);

	  cout << "Signal spectrum integral: " << flare_spectrum->Integral(Emin,Emax) << endl;

	  if(plot_results){
		TCanvas* speccanvas = new TCanvas();
		speccanvas->SetLogy();
		speccanvas->SetLogx();
		flare_spectrum->Draw();
	  }


  //-----------------------------------LC inputs (Interpolation/KDE)--------------------------------------

	if(KDE_mode == false){
		TString pathToReadInter;
		if(PIC_job == false) pathToReadInter = "../ReadEvents/Paper_Inter/";
		else pathToReadInter = "./";
		const TString fileName3 = pathToReadInter + "Inter_Paper_36ONevents_width.root";
		TFile* Vectorfile = TFile::Open(fileName3, "READ");
		if (!Vectorfile) {
			log_file << "ERROR: file " << Vectorfile << " not found..." << endl;
		}

		TRandom3 *ran = new TRandom3();
		ran->SetSeed(0);
		//Read the Tree
		TTree *t1 = (TTree*)Vectorfile->Get("t1");
		Double_t x_point, y_point, bin_width;
		t1->SetBranchAddress("x",&x_point);
		t1->SetBranchAddress("y",&y_point);
		t1->SetBranchAddress("width",&bin_width);

		//Fill vectors
		//--First point (Link with the Gaussian), no errors
		Double_t x_value = -90.;
		x_points.push_back(x_value);
		y_points.push_back(gaus_function(x_value));
		y_mod_points.push_back(gaus_function(x_value));

		  //--Interpolation points
		  Long64_t nentries = t1->GetEntries();
		  for(Long64_t i=0;i<nentries;i++){
			  t1->GetEntry(i);
			  x_points.push_back(x_point);
//			  cout << "------y_point: " << y_point << endl;
//			  cout << "Events: " << y_point*bin_width << endl;
			  Double_t new_events = ran->Poisson(y_point*bin_width);
//			  cout << "New Events: " << new_events << endl;
			  Double_t y_mod = new_events/bin_width;
//			  cout << "y_mod: " << y_mod << endl;
			  y_points.push_back(y_point);
			  y_mod_points.push_back(y_mod);
			  bin_widths.push_back(bin_width);
		  }

		  //--Last point (Link with the Gaussian), no error
		  Double_t x_value2 = x_points[nentries]-x_value;
		  x_points.push_back(x_value2);
		  y_points.push_back(gaus_function(x_value2));
		  y_mod_points.push_back(gaus_function(x_value2));



		  log_file << "Number of entries: " << nentries << endl;
		  log_file << "First time point is: " << x_points[1] << " s" << endl;
		  log_file << "First time point value: " << y_points[1] << " s" << endl;
		  log_file << "Last time point is: " << x_points[nentries] << " s" << endl;
		  log_file << "Last time point value: " << y_points[nentries] << " s" << endl;
		  log_file << "New First time point is: " << x_points[0] << " s" << endl;
		  log_file << "New First time point value: " << y_points[0] << " s" << endl;
		  log_file << "New Last time point is: " << x_points[nentries+1] << " s" << endl;
		  log_file << "New Last time point value: " << y_points[nentries+1] << " s" << endl;
		  inter_min = x_points[0];
		  inter_max = x_points[nentries+1];
		  tmin = x_points[1];
		  tmax = x_points[nentries];
		  log_file << "tmin: " << tmin << " tmax: " << tmax << endl;


		  //Do and draw Akima interpolation
		  int_aspline_original = new ROOT::Math::Interpolator(x_points, y_points, ROOT::Math::Interpolation::kAKIMA);
		  int_aspline = new ROOT::Math::Interpolator(x_points, y_mod_points, ROOT::Math::Interpolation::kAKIMA);

		  Int_t number_of_points = 10000;
		  general =  new TGraph(number_of_points);
		  general_original =  new TGraph(number_of_points);
		  general->SetName("Akima interpolation");
		  general->SetTitle("Akima Interpolation");
		  Double_t step = 13000/number_of_points;
		  Double_t point = x_points[0];
		  Int_t iter=0;
		  while(point+iter*step<=x_points[nentries+1]){
			  general->SetPoint(iter, point+iter*step, int_aspline->Eval(point+iter*step));
			  general_original->SetPoint(iter, point+iter*step, int_aspline_original->Eval(point+iter*step));

			  iter++;
		  }

		  if(plot_results){
			  TCanvas* canvasinter =  new TCanvas();
			  canvasinter->Divide(3,1);
			  canvasinter->cd(1);
			  general->SetTitle("Modified Akima Interpolation");
			  general->GetXaxis()->SetTitle("Time(s)");
			  general->GetYaxis()->SetTitle("dN/dt");
			  general->SetLineColor(kBlack);
			  general->SetLineWidth(1);
			  general->Draw("AL");
			  canvasinter->cd(2);
			  general_original->SetTitle("Original Akima Interpolation");
			  general_original->SetLineColor(4);
			  general_original->Draw("AL");
			  canvasinter->cd(3);
			  general->Draw("AL");
			  general_original->Draw("same");
//			  time_events->Draw("same");


			  TCanvas* allcanvas = new TCanvas();
			  general->Draw("AL");
			  general_original->Draw("same");

			  TLegend* const2legend = new TLegend(0.6,0.7,0.9,0.9);
			  const2legend->AddEntry(general_original, "Original interpolation", "l");
			  const2legend->AddEntry(general, "Poisson modified interpolation", "l");
			  const2legend->Draw("");

			  TCanvas* oricanvas = new TCanvas();
			  general_original->SetTitle("");
			  general_original->GetXaxis()->SetTitle("Time(s)");
			  general_original->GetYaxis()->SetTitle("dN/dt");
			  general_original->Draw("AL");
			  general_original->Draw("same");

		  }


		  log_file << "Interpolation is ready" << endl;
	}

	if(KDE_mode == true){

		TString pathToReadKDE;
		if(PIC_job == false) pathToReadKDE = "../ReadEvents/Thesis_KDE/";
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

			cout << "KDE function is ready" << endl;

	}


	if(include_bkg == 1){
		//-----------------------------------Background events--------------------------------------


		//Energy distribution
		TString pathToReadBkg;
		if(PIC_job == false) pathToReadBkg = "../ReadEvents/Paper_Histos/";
		else pathToReadBkg = "./";
		const TString fileNameBkg = pathToReadBkg + "Histo_36OFF_Paper.root";
		TFile* BKgfile = TFile::Open(fileNameBkg, "READ");
		if (!BKgfile) {
			cout << "ERROR: file " << BKgfile << " not found..." << endl;
		}

		energy_dist_bkg_histo = dynamic_cast<TH1D*>(BKgfile->Get("energy_dist"));
		if (!energy_dist_bkg_histo) {
				cout << "ERROR: energy_dist_bkg object not found... " << endl;
		}

		energy_dist_bkg_histo->Rebin();

		Int_t numberEBinsBkg = energy_dist_bkg_histo->GetNbinsX();
		cout << "Number of energy bins in bkg: " << numberEBinsBkg << endl;
		energy_dist_bkg = new TGraph(numberEBinsBkg);
		for(Int_t i=1; i<=numberEBinsBkg; i++ ){
			energy_dist_bkg->SetPoint(i-1, energy_dist_bkg_histo->GetBinCenter(i), energy_dist_bkg_histo->GetBinContent(i));
	//		cout << "For bin " << i << " content: " <<  energy_dist_bkg_histo->GetBinContent(i) << endl;
		}

		Emax_bkg = energy_dist_bkg_histo->GetBinLowEdge(numberEBinsBkg+1);
		cout << "Emax for bkg: " << Emax_bkg << endl;

		if(plot_results){
			TCanvas* canvasEBkg =  new TCanvas();
			canvasEBkg->SetLogx();
			canvasEBkg->SetLogy();
			energy_dist_bkg->SetTitle("Background energy distribution");
			energy_dist_bkg->SetLineColor(kRed);
			energy_dist_bkg->SetLineWidth(2);
			energy_dist_bkg->Draw("");
			energy_dist_bkg_histo->Draw("same");
		}

			cout << "Background Energy is ready" << endl;
			integral_spectrum_bkg = energy_dist_bkg_histo->Integral();


			cout << "Background spectrum integral: " << integral_spectrum_bkg << endl;




		//Time distribution
		if(KDE_mode == false){
				TString pathToReadInter;
				if(PIC_job == false) pathToReadInter = "../ReadEvents/Paper_Inter/";
				else pathToReadInter = "./";
				const TString fileName3 = pathToReadInter + "Inter_Paper_36OFFevents_width.root";
				TFile* Vectorfile = TFile::Open(fileName3, "READ");
				if (!Vectorfile) {
					log_file << "ERROR: file " << Vectorfile << " not found..." << endl;
				}

				TRandom3 *ran = new TRandom3();
				ran->SetSeed(0);
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
	//				  cout << "----y_point: " << y_point << endl;
					  y_point = y_point/3.;
					  x_points_bkg.push_back(x_point);
					  y_points_bkg.push_back(y_point);
					  bin_widths_bkg.push_back(bin_width);
				  }

				  log_file << "Number of entries: " << nentries << endl;
				  log_file << "First time point is: " << x_points_bkg[0] << " s" << endl;
				  log_file << "First time point value: " << y_points_bkg[0] << " s" << endl;
				  log_file << "Last time point is: " << x_points_bkg[nentries-1] << " s" << endl;
				  log_file << "Last time point value: " << y_points_bkg[nentries-1] << " s" << endl;
				  inter_min_bkg = x_points_bkg[0];
				  inter_max_bkg = x_points_bkg[nentries-1];


				  //Do and draw Akima interpolation
				  int_aspline_background_original = new ROOT::Math::Interpolator(x_points_bkg, y_points_bkg, ROOT::Math::Interpolation::kAKIMA);

				  Int_t number_of_points = 10000;
				  general_original_bkg =  new TGraph(number_of_points);
				  general_original_bkg->SetName("Akima interpolation for Background");
				  general_original_bkg->SetTitle("Akima interpolation for Background");
				  Double_t step = 13000/number_of_points;
				  Double_t point = x_points_bkg[0];
				  Int_t iter=0;
				  while(point+iter*step<=x_points_bkg[nentries-1]){
					  general_original_bkg->SetPoint(iter, point+iter*step, int_aspline_background_original->Eval(point+iter*step));
					  iter++;
				  }

				  if(plot_results){
					  TCanvas* allcanvas = new TCanvas();
					  general_original_bkg->Draw("AL");

					  TLegend* const2legend = new TLegend(0.1,0.7,0.4,0.9);
					  const2legend->AddEntry(general_original_bkg, "Original interpolation", "l");
					  const2legend->Draw("");

				  }


				  log_file << "Background Interpolation is ready" << endl;

			}
			if(KDE_mode == true){


				TString pathToReadKDE;
				if(PIC_job == false) pathToReadKDE = "../ReadEvents/Thesis_KDE/";
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

				cout << "KDE Background function is ready" << endl;

			}

	}

	Normalization_Compute();
	cout << "Interpolation Signal: " << integral_inter_signal << endl;
	cout << "Interpolation Background: " << integral_inter_bkg << endl;

  //-----------------------------------IRFs inputs--------------------------------------

	TString pathToReadIRF;
	if(PIC_job == false) pathToReadIRF = "../ReadIRFs/";
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


	//----------Double Gaussian Fit----------
	//Const parameter
	Const1Erec = dynamic_cast<TH1D*>(IRFfile->Get("Const1_histo_Erec"));
	if (!Const1Erec) {
			cout << "ERROR: Const1Erec object not found... " << endl;
	}
	Const2Erec = dynamic_cast<TH1D*>(IRFfile->Get("Const2_histo_Erec"));
	if (!Const2Erec) {
			cout << "ERROR: Const2Erec object not found... " << endl;
	}
	Const0Erec = dynamic_cast<TH1D*>(IRFfile->Get("Const0_histo_Erec"));
	if (!Const0Erec) {
			cout << "ERROR: Const0Erec object not found... " << endl;
	}



	if(plot_results){
		TCanvas* const2canvas = new TCanvas();
		const2canvas->SetLogx();
		Const1Erec->SetStats(0);
		Const1Erec->SetLineColor(kBlue);
		Const1Erec->Draw("");
		Const2Erec->SetLineColor(kRed);
		Const2Erec->Draw("same");
		Const0Erec->SetLineColor(kBlack);
		Const0Erec->Draw("same");

		TLegend* const2legend = new TLegend(0.6,0.7,0.9,0.9);
		const2legend->AddEntry(Const1Erec, "Const paratermer, Gauss 1", "l");
		const2legend->AddEntry(Const2Erec, "Const paratermer, Gauss 2", "l");
		const2legend->AddEntry(Const0Erec, "Const paratermer, 1 Gauss fit", "l");
		const2legend->Draw("");
	}


	//Mean parameter
	Mean1Erec = dynamic_cast<TH1D*>(IRFfile->Get("Mean1_histo_Erec"));
	if (!Mean1Erec) {
			cout << "ERROR: Mean1Erec object not found... " << endl;
	}
	Mean2Erec = dynamic_cast<TH1D*>(IRFfile->Get("Mean2_histo_Erec"));
	if (!Mean2Erec) {
			cout << "ERROR: Mean2Erec object not found... " << endl;
	}
	Mean0Erec = dynamic_cast<TH1D*>(IRFfile->Get("Mean0_histo_Erec"));
	if (!Mean0Erec) {
			cout << "ERROR: Mean0Erec object not found... " << endl;
	}

	if(plot_results){
		TCanvas* mean2canvas = new TCanvas();
		mean2canvas->SetLogx();
		Mean1Erec->SetStats(0);
		Mean1Erec->SetLineColor(kBlue);
		Mean1Erec->Draw("");
		Mean2Erec->SetLineColor(kRed);
		Mean2Erec->Draw("same");
		Mean0Erec->SetLineColor(kBlack);
		Mean0Erec->Draw("same");


		TLegend* main2legend = new TLegend(0.6,0.7,0.9,0.9);
		main2legend->AddEntry(Mean1Erec, "Mean paratermer, Gauss 1", "l");
		main2legend->AddEntry(Mean2Erec, "Mean paratermer, Gauss 2", "l");
		main2legend->AddEntry(Mean0Erec, "Mean paratermer, 1 Gauss fit", "l");
		main2legend->Draw("");
	}

	//Sigma parameter
	Sigma1Erec = dynamic_cast<TH1D*>(IRFfile->Get("Sigma1_histo_Erec"));
	if (!Sigma1Erec) {
			cout << "ERROR: Sigma1Erec object not found... " << endl;
	}
	Sigma2Erec = dynamic_cast<TH1D*>(IRFfile->Get("Sigma2_histo_Erec"));
	if (!Sigma2Erec) {
			cout << "ERROR: Sigma2Erec object not found... " << endl;
	}
	Sigma0Erec = dynamic_cast<TH1D*>(IRFfile->Get("Sigma0_histo_Erec"));
	if (!Sigma0Erec) {
			cout << "ERROR: Sigma0Erec object not found... " << endl;
	}

	if(plot_results){
		TCanvas* sigma2canvas = new TCanvas();
		sigma2canvas->SetLogx();
		Sigma2Erec->SetStats(0);
		Sigma2Erec->SetLineColor(kRed);
		Sigma2Erec->Draw("");
		Sigma1Erec->SetLineColor(kBlue);
		Sigma1Erec->Draw("same");
		Sigma0Erec->SetLineColor(kBlack);
		Sigma0Erec->Draw("same");


		TLegend* sigma2legend = new TLegend(0.6,0.7,0.9,0.9);
		sigma2legend->AddEntry(Sigma1Erec, "Sigma paratermer, Gauss 1", "l");
		sigma2legend->AddEntry(Sigma2Erec, "Sigma paratermer, Gauss 2", "l");
		sigma2legend->AddEntry(Sigma0Erec, "Sigma paratermer, 1 Gauss fit", "l");
		sigma2legend->Draw("");
	}

	if(plot_results){
		//Test for a given energy
		Double_t test_energy = 2000; //In GeV

		Double_t Const = Const0Erec->GetBinContent(Const0Erec->FindBin(test_energy));
		Double_t mean = Mean0Erec->GetBinContent(Mean0Erec->FindBin(test_energy));
		Double_t sigma = Sigma0Erec->GetBinContent(Sigma0Erec->FindBin(test_energy));
		Double_t Erec_down_simple = (mean+1)*test_energy-4*test_energy*sigma; //Limits of the integral
		if(Erec_down_simple<0)Erec_down_simple = 0;
		Double_t Erec_up_simple = (mean+1)*test_energy+4*test_energy*sigma;
		log_file << "Limits simple Gaus" << endl;
		log_file << "UP: " << Erec_down_simple << " DOWN: " << Erec_up_simple << endl;
		TF1* simpleGaus =  new TF1("Simple Gauss", "gaus", Erec_down_simple, Erec_up_simple);
		simpleGaus->SetNpx(10000);

		simpleGaus->SetParameters(Const, (1-mean)*test_energy, sigma*test_energy);

		Double_t const1 = Const1Erec->GetBinContent(Const1Erec->FindBin(test_energy));
		Double_t const2 = Const2Erec->GetBinContent(Const2Erec->FindBin(test_energy));
		Double_t mean1 = Mean1Erec->GetBinContent(Mean1Erec->FindBin(test_energy));
		Double_t mean2 = Mean2Erec->GetBinContent(Mean2Erec->FindBin(test_energy));
		Double_t sigma1 = Sigma1Erec->GetBinContent(Sigma1Erec->FindBin(test_energy));
		Double_t sigma2 = Sigma2Erec->GetBinContent(Sigma2Erec->FindBin(test_energy));
		Double_t Erec_down_doble = (mean1+1)*test_energy-4*test_energy*sigma1; //Limits of the integral
		if(Erec_down_doble<0)Erec_down_doble = 0;
		Double_t Erec_up_doble = (mean1+1)*test_energy+4*test_energy*sigma1;
		log_file << "Limits doble Gaus" << endl;
		log_file << "UP: " << Erec_down_doble << " DOWN: " << Erec_up_doble << endl;
		TF1* dobleGaus =  new TF1("Doble Gauss", Double_Gaus, Erec_down_doble, Erec_up_doble, 6);
		dobleGaus->SetNpx(10000);

		dobleGaus->SetParameters(const1, (1-mean1)*test_energy, sigma1*test_energy, const2, (1-mean2)*test_energy, sigma2*test_energy);

		TCanvas* energytestcanvas =  new TCanvas();
//		energytestcanvas->SetLogx();
		dobleGaus->SetLineColor(kRed);
		dobleGaus->Draw("");
		simpleGaus->SetLineColor(kBlue);
		simpleGaus->Draw("same");

		Double_t Es_simple = simpleGaus->GetRandom();
		Double_t Es_doble = dobleGaus->GetRandom();

		log_file << "Es_simple: " << Es_simple << " Es_doble: " << Es_doble << endl;
	}

	//-----------------------------------EBL table--------------------------------------
	Double_t redshift = 0.031; //Mrk421
	Float_t z[39] = {0.01, 0.02526316, 0.04052632, 0.05578947, 0.07105263, 0.08631579, 0.10157895, 0.11684211, 0.13210526, 0.14736842, 0.16263158, 0.17789474, 0.19315789, 0.20842105, 0.22368421, 0.23894737, 0.25421053, 0.26947368, 0.28473684, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.2, 1.4, 1.6, 1.8, 2.};

	//Read dominguez table
	ifstream f("./tau_dominguez11.dat");
	if (!f.is_open())
	{
		cout << "Cannot open EBL file" << endl;
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


	//-----------------------------------LIKELIHOOD FIT--------------------------------------

	cout.precision(10);

	Double_t Like_res = 0.;
//	if (KDE_mode == false)Like_res = Likelihoodtest_Inter();
//	else if (KDE_mode == true)Like_res = Likelihoodtest();
//	cout << "Initial Likelihood value: " << Like_res << endl;


	TMinuit *myMinuit=new TMinuit(1);
	if (KDE_mode == false) myMinuit->SetFCN(fcn_Inter);
	if (KDE_mode == true) myMinuit->SetFCN(fcn_KDE);


	Double_t arglist[10];
	Int_t ierflg = 0;

	//---Parameter definition---
	Double_t QGdown=-QGpar*50;
	Double_t QGup=QGpar*50;
	myMinuit->mnparm(0, "QGpar", QGpar, 10., QGdown, QGup,ierflg);


	//---Fix parameters---
	//	myMinuit->FixParameter(0);


	myMinuit->SetErrorDef(1);
//	Double_t error;


//	Scan and draw parameters
	myMinuit->Command("SCAn 1 10 -5 5"); //For interpolation
//	myMinuit->Command("SCAn 1 5 0 5"); //For KDE
	TObject* gr1= gMinuit->GetPlot();


	if(plot_results){
		TCanvas* parameters = new TCanvas();
		parameters->SetTitle("QGparameter_MinuitScan");
		gr1->Draw("alp*");
	}



	//Remember to indicate the arguments for every Minuit command!
	arglist[0] = 0;
	myMinuit->mnexcm("SET STRategry", arglist, 1, ierflg);
	cout << "ierflg after SETSTRatregy: " << ierflg << endl;
	arglist[0] = 500;
//	myMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
//	cout << "ierflg after MIGRAD: " << ierflg << endl;
	myMinuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
	cout << "ierflg after SIMPLEX: " << ierflg << endl;
	//Filter the cases that did not converge
	if(ierflg!=0) return 0;


	//Extract results
//	myMinuit->mnexcm("MINOS", arglist ,2,ierflg);
	Double_t ln0, edm, errdef;
	Int_t nvpar, nparx, icstat;
	myMinuit->mnstat(ln0,edm,errdef,nvpar,nparx,icstat);
    TString pname;
    double plbound,pubound;
    int pvari;
	Double_t alpha, erroralpha;
    myMinuit->mnpout(0,pname,alpha,erroralpha,plbound,pubound,pvari);

	// Create a TTree to save results
    Double_t UL,LL;
	TTree tree("t1","Tree with simple variables");
	tree.Branch("alpha",&alpha,"px/D");
	tree.Branch("error",&erroralpha,"py/D");
	tree.Branch("UL",&UL,"pz/D");
	tree.Branch("LL",&LL,"pa/D");
	cout << "Tree created" << endl;

	//Draw and scale correctly the curve: Start in minimum alpha and go to the sides.
	Double_t conf = 4.; //confidence level, in this case 2-sigma two-sided
	TGraph* Like_curve =  new TGraph();
	Double_t right_arm[500] = {0};
	Double_t left_arm[500] = {0};
	Double_t scaling = ln0;
	Double_t delta_L = 0.;

	//Size of the steps to move to the sides
	Double_t step_alpha = 0.;
	if(QGorder == 1) step_alpha = 0.05; //0.01(0.03 quick) for 9Inter //5. for KDE //0.05 for 25Inter //0.2 for 100Inter //0.1 for 64
	if(QGorder == 2) step_alpha = 0.02;//0.003(0.009 quick) for 9Inter //0.5 for KDE //0.01 for 25Inter //0.15 for 100Inter //0.05 for 64
	QGpar = alpha;
	log_file << "alpha is: " << QGpar << endl;
	log_file << "scale is: " << scaling << endl;
	log_file << "step is: " << step_alpha << endl;

	Int_t index = 0;
	Int_t points_left=0;
	Int_t points_right = 0;

	while(delta_L <= conf){
		log_file << "Index: " << index << endl;
		log_file << "QGpar: " << QGpar << endl;
		if(KDE_mode == false)delta_L = Likelihoodtest_Inter()-scaling;
		else if(KDE_mode == true) delta_L = Likelihoodtest()-scaling;
		log_file << "Delta_L: " << delta_L << endl;
		right_arm[index] = delta_L;
		index++;
		UL = QGpar;
		if(index>199) break;
		QGpar+=step_alpha;

	}
	points_right = index-1;
	cout << "Points in the right: " << points_right+1 << endl;
	cout << "UL: " << UL << endl;
	log_file << "Points in the right: " << points_right+1 << endl;
	log_file << "UL: " << UL << endl;


	QGpar=alpha;
	index=0;
	delta_L = 0.;
//	cout << "alpha is: " << QGpar << endl;
	while(delta_L <= conf){
		QGpar-=step_alpha;
		log_file << "QGpar: " << QGpar << endl;
		if(KDE_mode == false) delta_L = Likelihoodtest_Inter()-scaling;
		else if(KDE_mode == true) delta_L = Likelihoodtest()-scaling;
		log_file << "Delta_L: " << delta_L << endl;
		log_file << "Index: " << index << endl;
		LL = QGpar;
		left_arm[index] = delta_L;
		index++;
		if(index>199) break;
	}
	points_left = index-1;
	cout << "Points in the left: " << points_left+1 << endl;
	cout << "LL: " << LL << endl;
	log_file << "Points in the left: " << points_left+1 << endl;
	log_file << "LL: " << LL << endl;
	for(Int_t i=0; i<=points_left+points_right+1; i++){
		if(i<=points_left){
			log_file << "QGpar " << QGpar << endl;
			log_file << "DeltaL " << left_arm[points_left-i] << endl;

			Like_curve->SetPoint(i, QGpar, left_arm[points_left-i]);
			QGpar+=step_alpha;
		}
		else{
//			cout << "QGpar " << QGpar << endl;
//			cout << "DeltaL " << right_arm[i-index] << endl;
			Like_curve->SetPoint(i, QGpar, right_arm[i-index]);
			QGpar+=step_alpha;
		}
	}

	if(plot_results){
		TCanvas* Like_canvas = new TCanvas;
		Like_curve->SetTitle("Likelihood curve");
		Like_curve->GetXaxis()->SetTitle("#alpha");
		Like_curve->GetYaxis()->SetTitle("-2 #Delta L");
		Like_curve->SetLineColor(2);
		Like_curve->SetLineWidth(2);
		Like_curve->Draw("AL");

		Like_canvas->SaveAs("Curve.png");
	}


	tree.Fill();
	//Create de TFile
	TFile* Likelihood_outputfile = new TFile("LikelihoodResult_QuadraticCase_DATA.root","recreate");
	tree.Write();
	Like_curve->Write();

	cout << "Root file created"  << endl;


  return 0;

}


void PDFvsEvsT(Double_t fixed_E, Double_t fixed_t){

	const Int_t numberOfPoints = 2000;
	TGraph* LikevsE = new TGraph(numberOfPoints);
	TGraph* Likevst = new TGraph(numberOfPoints);
	Double_t E_step = (Emax-Emin)/numberOfPoints;
	Double_t t_step = (tmax-tmin)/numberOfPoints;


	for (Int_t i=0; i<numberOfPoints; i++){
		if(ResolutionInc == false){
			LikevsE->SetPoint(i, Emin+i*E_step, (1/100.)*PDFtest_Inter(fixed_t,Emin+i*E_step));
			Likevst->SetPoint(i, tmin+i*t_step, (1/8500.)*PDFtest_Inter(tmin+i*t_step, fixed_E));
		}
		else{
			LikevsE->SetPoint(i, Emin+i*E_step, PDFtestRes_Inter(fixed_t,Emin+i*E_step));
			Likevst->SetPoint(i, tmin+i*t_step, PDFtestRes_Inter(tmin+i*t_step, fixed_E));
		}
	}

	//Draw final plots
	TCanvas* Likevscanvas = new TCanvas();
	Likevscanvas->Divide(2,1);
	Likevscanvas->cd(1);
	Likevscanvas->cd(1)->SetTitle("Fixed time");
	Likevscanvas->cd(1)->SetLogx();
	Likevscanvas->cd(1)->SetLogy();
	LikevsE->SetTitle("PDF vs E");
	LikevsE->Draw("AL");
//	flare_spectrum->Draw("same");
	energy_events->Draw("same");
	Likevscanvas->cd(2);
	Likevscanvas->cd(2)->SetTitle("Fixed energy");
	Likevst->SetTitle("PDF vs t");
	Likevst->Draw("AL");
	time_events->Draw("same");

	//More canvas
	TCanvas* Likevscanvas2 = new TCanvas();
	Likevscanvas2->Divide(2,1);
	Likevscanvas2->cd(1);
	Likevscanvas2->cd(1)->SetTitle("Fixed time");
	Likevscanvas2->cd(1)->SetLogx();
	Likevscanvas2->cd(1)->SetLogy();
	energy_events->Draw("");
	Likevscanvas2->cd(2);
	Likevscanvas2->cd(2)->SetTitle("Fixed energy");
	time_events->Draw("");


}

void DrawLikelihood(Double_t init_alpha, Double_t final_alpha){

	Int_t numberOfPoints = 30;
	Double_t sizeOfSteps = (final_alpha-init_alpha)/numberOfPoints;
	TGraph* Likevsalpha = new TGraph(numberOfPoints);
	for(Int_t i=0; i<numberOfPoints; i++){
		QGpar=init_alpha+i*sizeOfSteps;
		Double_t Like_val = 0.;
		if(KDE_mode == false) Like_val = Likelihoodtest_Inter();
		if(KDE_mode == true) Like_val = Likelihoodtest();
		cout << "Likelihood value: " << Like_val << endl;
		Likevsalpha->SetPoint(i, QGpar, Like_val);
	}


	TCanvas* Likecanvas = new TCanvas();
	Likevsalpha->Draw("APL*");


}

void DrawLikelihood_Energy(Double_t min_Energy, Double_t max_Energy, Double_t init_alpha, Double_t final_alpha){

	Int_t numberOfPoints = 30;
	Double_t sizeOfSteps = (final_alpha-init_alpha)/numberOfPoints;
	TGraph* Likevsalpha = new TGraph(numberOfPoints);
	for(Int_t i=0; i<numberOfPoints; i++){
		QGpar=init_alpha+i*sizeOfSteps;
		Double_t Like_val = Likelihoodtest_Inter_Energy(min_Energy, max_Energy);
		cout << "Likelihood value: " << Like_val << endl;
		Likevsalpha->SetPoint(i, QGpar, Like_val);
	}


	TCanvas* Likecanvas = new TCanvas();
	Likevsalpha->Draw("APL*");


}

void DrawResultLikelihood(){

//Read initial histogram to get the bin distribution
TFile* HistoFile = TFile::Open("/home/lnogues/eclipse-workspace/LIVanalysis/likelihood_code/ReadEvents/FinalEventLists/LC_01-20TeV_ONhistograms_Mrk421_2014.root", "READ");
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


//Get histo with energy and interpolation
TH1D* time_histo_original = new TH1D("time_histo_original","time_histo_original",numberOfBins,arrayBins);
TH1D* time_histo_corr = new TH1D("time_histo_corr","time_histo_corr",numberOfBins,arrayBins);


//Move back the events
Double_t t_corr[20000] = {0};


Double_t QGDelay_result= QGfixedPart*result_alpha/PlanckScale;
cout << "resulting alpha is: " << result_alpha << " the delay is: " << QGDelay_result << endl;
for(Int_t i=0; i<numberOfEvents; i++){
	t_corr[i] = t[i]-(QGDelay_result*E[i]); //Energy in GeV
	//t[index],sim_flare_time->GetBinWidth(1)/sim_flare_time->GetBinWidth(sim_flare_time->FindBin(t[index]))
	time_histo_original->Fill(t[i],time_histo_original->GetBinWidth(1)/time_histo_original->GetBinWidth(time_histo_original->FindBin(t[i])));
	time_histo_corr->Fill(t_corr[i], time_histo_corr->GetBinWidth(1)/time_histo_corr->GetBinWidth(time_histo_corr->FindBin(t_corr[i])));
}


//Create a interpolation with the corrected events (x,y from corrected histogram)
std::vector<double> x_corr;
std::vector<double> y_corr;
for(Int_t i=1; i<=numberOfBins; i++){
	x_corr.push_back(time_histo_corr->GetBinCenter(i));
	y_corr.push_back(time_histo_corr->GetBinContent(i));
//	cout << "x_corr: " << x_corr.at(i-1) << endl;

}
ROOT::Math::Interpolator int_aspline(x_corr, y_corr, ROOT::Math::Interpolation::kAKIMA);
Int_t number_of_points = 4000;
TGraph* new_inter = new TGraph(number_of_points);
 Double_t step = 12500/number_of_points;
 Double_t point = x_corr[0];
  Int_t iter=0;
  while(inter_min+iter*step<x_corr[numberOfBins-2]){
	  new_inter->SetPoint(iter, point+iter*step, int_aspline.Eval(point+iter*step));
//	  cout << "inter_min+iter*step: " << inter_min+iter*step << endl;
	  iter++;
  }



//Draw final result
TCanvas* Resultcanvas = new TCanvas();
//time_histo_original->SetStats(0);
//time_histo_original->SetLineColor(kBlack);
//time_histo_original->SetLineWidth(2);
//time_histo_original->Draw();
//time_histo_corr->SetLineColor(kRed);
//time_histo_corr->Draw("same");
general->SetLineColor(kRed);
new_inter->SetLineColor(1);
general->SetLineWidth(2);
general->Draw("");
new_inter->Draw("same");


TLegend* leg_result = new TLegend(0.1,0.7,0.4,0.9);
//leg_result->AddEntry(time_histo_original, "Simulated time", "l");
//leg_result->AddEntry(time_histo_corr, "Corrected time", "l");
leg_result->AddEntry(general, "Original Interpolation", "l");
leg_result->AddEntry(new_inter, "Resulting Interpolation", "l");
leg_result->Draw();

}

void DrawNormalization(){

	//Draw the LC normalization function vs alpha
	TGraph* norm_graph = new TGraph();
	Int_t index = 0;

	for(Double_t alpha=-20; alpha<=20; alpha+=0.5){
		QGpar = alpha;

		 // calculate normalization
		  Double_t Normalization = 0.;
		  Double_t NormalizationBkg = 0.;

		  const Double_t binsize  = TMath::Log(Emax/Emin)/nbinsE;

		  for (Int_t i=0;i<nbinsE;i++){
			  for (Int_t j=0;j<nbinsT;j++){
				  Double_t dE     = Emin * (TMath::Exp(binsize * (i+1)) - TMath::Exp(binsize*i));
				  Double_t dT     = (tmax-tmin)/nbinsT;
//				  cout << "dE: " << dE << endl;
//				  cout << "dT " << dT << endl;
				  Double_t energy = Emin * TMath::Exp(binsize * (i+0.5));
				  Double_t time   = tmin + (j+0.5)*dT;
//				  cout << "Energy: " << energy << endl;
//				  cout << "Time " << time << endl;
				  Double_t n=0.;
				  Double_t n_bkg=0.;
				  if(ResolutionInc == false){
					 if(KDE_mode == false) n=PDFtest_Inter(time,energy)*wobble_condition(time);
					 if(KDE_mode == true) n=PDFtest(time,energy)*wobble_condition(time);
				  }
				  if(ResolutionInc == true){
					  if(KDE_mode == false)n= PDFtestRes_Inter(time,energy)*wobble_condition(time);
					  if(KDE_mode == true)n= PDFtestRes(time,energy)*wobble_condition(time);
				  }
				  if(include_bkg==1) n_bkg = PDFeventBkg_Inter(time,energy)*wobble_condition(time);
//				  cout << "n: " << n << endl;
//				  cout << "n_bkg: " << n_bkg << endl;
				  Normalization += n * dE * dT;
				  NormalizationBkg += n_bkg * dE * dT;

			  }
		  }

//		  Normalization=TMath::Log(Normalization);
//		  cout << "Normalization: " << Normalization << endl;
		  cout << "NormalizationBkg: " << NormalizationBkg << endl;

		  //Only plot the signal normalization (Bkg normalization does not depend on alpha)
		  norm_graph->SetPoint(index, QGpar, Normalization);
		  cout << "QGpar" << QGpar << endl;
		  cout << "Normalization: " << Normalization << endl;
		  index++;
		}

	TCanvas* norm_canvas = new TCanvas();
	norm_graph->Draw("ACP*");

//	TCanvas* int_canvas = new TCanvas();
//	general->Draw();

}

void DrawPDF(Int_t event_id, Double_t init_alpha, Double_t final_alpha){
	Int_t numberOfPoints = 100;
	cout << "Event -> Energy time: " << E[event_id] << " " << t[event_id] << endl;
	Double_t sizeOfSteps = (final_alpha-init_alpha)/numberOfPoints;
	TGraph* PDFvsalpha = new TGraph(numberOfPoints);
	for(Int_t i=0; i<numberOfPoints; i++){
		QGpar=init_alpha+i*sizeOfSteps;
		Double_t PDF_val = 0.;
		if(ResolutionInc == false) PDF_val= PDFtest_Inter(t[event_id], E[event_id]);
		else PDF_val= PDFtestRes_Inter(t[event_id], E[event_id]);
		PDFvsalpha->SetPoint(i, QGpar, PDF_val);
		cout << "For alpha " << QGpar << " PDF is " << PDF_val<< endl;
	}

	TCanvas* PDFcanvas = new TCanvas();
	PDFvsalpha->Draw("APL*");
}

void EdgesShift(Double_t alpha, Double_t Energy){
	Double_t t_max_shift = tmin;
	Double_t t_min_shift = tmax;
	Double_t delay = GetDelay(Energy, alpha);
	t_max_shift += delay;
	t_min_shift += delay;

	cout << "QGpar = " << alpha << " and energy(GeV) = " << Energy << endl;
	cout << "Inter max= " << tmax << " New max time: " << t_max_shift << " Delta: " << t_max_shift-inter_max << endl;
	cout << "Inter min= " << tmin << " New min time: " << t_min_shift << " Delta: " << t_min_shift-inter_min << endl;

}

void DrawEnergyTimeEvents(){
	TH1D* energy_histo= new TH1D("Energy", "Energy", 30, Emin, Emax);
	TH1D* time_histo= new TH1D("Time", "Time", 30, tmin, tmax);

	for(Int_t i=0; i<numberOfEvents; i++){
		energy_histo->Fill(E[i]);
		time_histo->Fill(t[i]);
	}

	TCanvas* canvas = new TCanvas();
	canvas->Divide(2,1);
	canvas->cd(1);
	energy_histo->Draw();
	canvas->cd(2);
	time_histo->Draw();
}

void DrawLightcurve(Double_t timemin, Double_t timemax){
	Int_t numberOfPoints = 1000;
	Double_t step = (timemax-timemin)/numberOfPoints;
	TGraph* LC = new TGraph (numberOfPoints);
	TGraph* Gaussian =  new TGraph(numberOfPoints);
	Double_t time = 0.;

	for(Int_t i=0; i<numberOfPoints; i++){
		time = timemin + i*step;
		LC->SetPoint(i, time, lc_function(time));
		Gaussian->SetPoint(i, time, gaus_function(time));
	}

	TCanvas* canvas = new TCanvas();
	canvas->Divide(2,1);
	canvas->cd(1);
	LC->Draw();
	canvas->cd(2);
	Gaussian->Draw();
}

void DrawBkgLightcurve(){
	Int_t numberOfPoints = 100;
	Double_t step = (tmax-tmin)/numberOfPoints;
	TGraph* LC = new TGraph (numberOfPoints);
	Double_t time = 0.;

	for(Int_t i=0; i<numberOfPoints; i++){
		time = tmin + i*step;
		LC->SetPoint(i, time, int_aspline_background_original->Eval(time));
	}

	TCanvas* canvas = new TCanvas();
	LC->Draw();
}



