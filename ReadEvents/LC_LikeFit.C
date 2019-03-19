/*
 * LC_LikeFit.C
 *
 *  Created on: Apr 18, 2017
 *      Author: lnogues
 *
 *      Code to use a binned Likelihood to fit the LC described by the data
 */

#include<iostream>
#include<fstream>
#include<istream>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <istream>
#include <TVectorD.h>
#include <TMinuit.h>
#include <TLatex.h>

using namespace std;

Bool_t cut_fit = false;

const Int_t Np = 100000;   // Number of events (hardcoded)
Double_t Ep[Np]; //In GeV
Double_t tp[Np]; //In s

TH1D* timehisto;
TH1D* residualhisto;
TH1D* residualsquarehisto;
TH1D* new_lightcurve;
static TMinuit *minuit = NULL;
TMinuit *myMinuit = NULL;

const Int_t binsnum = 500; //A big value
Int_t numberOfBins = 0; //The real one
Int_t new_numberOfBins = 0; //After cuts
Double_t bin[binsnum+1];
Double_t binwidth=60.; // In seconds
Double_t residual[binsnum]={0};
Double_t residual_squared[binsnum]={0};

//Wobble info
static const  Int_t numberOfWobbles = 14; //Number of wobbles
static const  Int_t numberOfPoints = 2*numberOfWobbles;

static const  Double_t mintimeMJD = 1.77293506944444380e+03; //From wobble info
static const Double_t maxtimeMJD = 1.77308576388889196e+03; //From wobble info

//Wobble limits
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





Double_t mintime = 0;
Double_t maxtime = (maxtimeMJD-mintimeMJD)*24.*60.*60.;
Double_t pvalue = 0;


static Double_t One_gaus(Double_t* x, Double_t*  par){

	//Wobble holes here
	if(x[0]>hole[1] && x[0]<hole[2]) return 0;
	if(x[0]>hole[3] && x[0]<hole[4]) return 0;
	if(x[0]>hole[5] && x[0]<hole[6]) return 0;
	if(x[0]>hole[7] && x[0]<hole[8]) return 0;
	if(x[0]>hole[9] && x[0]<hole[10]) return 0;
	if(x[0]>hole[11] && x[0]<hole[12]) return 0;
	if(x[0]>hole[13] && x[0]<hole[14]) return 0;
	if(x[0]>hole[15] && x[0]<hole[16]) return 0;
	if(x[0]>hole[17] && x[0]<hole[18]) return 0;
	if(x[0]>hole[19] && x[0]<hole[20]) return 0;
	if(x[0]>hole[21] && x[0]<hole[22]) return 0;
	if(x[0]>hole[23] && x[0]<hole[24]) return 0;
	if(x[0]>hole[25] && x[0]<hole[26]) return 0;


	return par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[3],2));
}


static Double_t One_Asymgaus(Double_t* x, Double_t*  par){

	//Wobble holes here
	if(x[0]>hole[1] && x[0]<hole[2]) return 0;
	if(x[0]>hole[3] && x[0]<hole[4]) return 0;
	if(x[0]>hole[5] && x[0]<hole[6]) return 0;
	if(x[0]>hole[7] && x[0]<hole[8]) return 0;
	if(x[0]>hole[9] && x[0]<hole[10]) return 0;
	if(x[0]>hole[11] && x[0]<hole[12]) return 0;
	if(x[0]>hole[13] && x[0]<hole[14]) return 0;
	if(x[0]>hole[15] && x[0]<hole[16]) return 0;
	if(x[0]>hole[17] && x[0]<hole[18]) return 0;
	if(x[0]>hole[19] && x[0]<hole[20]) return 0;
	if(x[0]>hole[21] && x[0]<hole[22]) return 0;
	if(x[0]>hole[23] && x[0]<hole[24]) return 0;
	if(x[0]>hole[25] && x[0]<hole[26]) return 0;


	if(x[0]<par[2]) return par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[3],2));
	else return par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[4],2));
}

static Double_t One_Asymgaus_more(Double_t* x, Double_t*  par){

	//Wobble holes here
	if(x[0]>hole[1] && x[0]<hole[2]) return 0;
	if(x[0]>hole[3] && x[0]<hole[4]) return 0;
	if(x[0]>hole[5] && x[0]<hole[6]) return 0;
	if(x[0]>hole[7] && x[0]<hole[8]) return 0;
	if(x[0]>hole[9] && x[0]<hole[10]) return 0;
	if(x[0]>hole[11] && x[0]<hole[12]) return 0;
	if(x[0]>hole[13] && x[0]<hole[14]) return 0;
	if(x[0]>hole[15] && x[0]<hole[16]) return 0;
	if(x[0]>hole[17] && x[0]<hole[18]) return 0;
	if(x[0]>hole[19] && x[0]<hole[20]) return 0;
	if(x[0]>hole[21] && x[0]<hole[22]) return 0;
	if(x[0]>hole[23] && x[0]<hole[24]) return 0;
	if(x[0]>hole[25] && x[0]<hole[26]) return 0;

	//Parameters
	//par[0] baseline
	//par[1] weigh asym
	//par[2] mean asym
	//par[3] sigma left
	//par[4] sigma right
	//par[5] weigh gaus
	//par[6] mean gaus
	//par[7] sigma gaus
	//par[8] weigh gaus2
	//par[9] mean gaus2
	//par[10] sigma gaus2


	if(x[0]<par[2]){
		Double_t asym = par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[3],2));
		Double_t simgaus = par[5]*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[7],2));
//		Double_t simgaus2 = par[8]*TMath::Exp(-0.5*TMath::Power((x[0]-par[9])/par[10],2));
//		Double_t simgaus3 = par[11]*TMath::Exp(-0.5*TMath::Power((x[0]-par[12])/par[13],2));
//		Double_t simgaus4 = par[14]*TMath::Exp(-0.5*TMath::Power((x[0]-par[15])/par[16],2));
		return asym+simgaus;
//		return asym+simgaus+simgaus2;
//		return asym+simgaus+simgaus2+simgaus3;
//		return asym+simgaus+simgaus2+simgaus3+simgaus4;
	}
	else{
		Double_t asym = par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[4],2));
		Double_t simgaus = par[5]*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[7],2));
//		Double_t simgaus2 = par[8]*TMath::Exp(-0.5*TMath::Power((x[0]-par[9])/par[10],2));
//		Double_t simgaus3 = par[11]*TMath::Exp(-0.5*TMath::Power((x[0]-par[12])/par[13],2));
//		Double_t simgaus4 = par[14]*TMath::Exp(-0.5*TMath::Power((x[0]-par[15])/par[16],2));
		return asym+simgaus;
//		return asym+simgaus+simgaus2;
//		return asym+simgaus+simgaus2+simgaus3;
//		return asym+simgaus+simgaus2+simgaus3+simgaus4;
	}
}

static Double_t One_Asymgaus_more_Save(Double_t* x, Double_t*  par){

	//This is the function that should be saved in the root file!

	if(x[0]<par[2]){
		Double_t asym = par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[3],2));
		Double_t simgaus = par[5]*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[7],2));
		Double_t simgaus2 = par[8]*TMath::Exp(-0.5*TMath::Power((x[0]-par[9])/par[10],2));
		Double_t simgaus3 = par[11]*TMath::Exp(-0.5*TMath::Power((x[0]-par[12])/par[13],2));
		Double_t simgaus4 = par[14]*TMath::Exp(-0.5*TMath::Power((x[0]-par[15])/par[16],2));
//		return asym+simgaus;
//		return asym+simgaus+simgaus2;
//		return asym+simgaus+simgaus2+simgaus3;
		return asym+simgaus+simgaus2+simgaus3+simgaus4;
	}
	else{
		Double_t asym = par[0] + par[1]*TMath::Exp(-0.5*TMath::Power((x[0]-par[2])/par[4],2));
		Double_t simgaus = par[5]*TMath::Exp(-0.5*TMath::Power((x[0]-par[6])/par[7],2));
		Double_t simgaus2 = par[8]*TMath::Exp(-0.5*TMath::Power((x[0]-par[9])/par[10],2));
		Double_t simgaus3 = par[11]*TMath::Exp(-0.5*TMath::Power((x[0]-par[12])/par[13],2));
		Double_t simgaus4 = par[14]*TMath::Exp(-0.5*TMath::Power((x[0]-par[15])/par[16],2));
//		return asym+simgaus;
//		return asym+simgaus+simgaus2;
//		return asym+simgaus+simgaus2+simgaus3;
		return asym+simgaus+simgaus2+simgaus3+simgaus4;
	}
}

Double_t test_statistics(Double_t constant=0, Double_t weight=0, Double_t mean=0, Double_t sigma1=0., Double_t sigma2=0., Double_t weightG = 0., Double_t meanG = 0., Double_t sigmaG = 0., Double_t weightG2 = 0., Double_t meanG2 = 0., Double_t sigmaG2 = 0., Double_t weightG3 = 0., Double_t meanG3 = 0., Double_t sigmaG3 = 0., Double_t weightG4 = 0., Double_t meanG4 = 0., Double_t sigmaG4 = 0.){

	Double_t chiSquare = 0;
	Int_t numberOfCurrentBins = 0;
	Int_t fullBins = 0;

//	TF1* function = new TF1("function", One_gaus, mintime, maxtime, 4);
//	function->SetNpx(1000);
//	function->SetParameters(constant, weight, mean, sigma1);

	TF1* function = new TF1("function", One_Asymgaus, mintime, maxtime, 5);
	function->SetNpx(1000);
	function->SetParameters(constant, weight, mean, sigma1, sigma2);

//	TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 8);
//	function->SetNpx(1000);
//	function->SetParameters(constant, weight, mean, sigma1, sigma2, weightG, meanG, sigmaG);

//	TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 11);
//	function->SetNpx(1000);
//	function->SetParameters(constant, weight, mean, sigma1, sigma2, weightG, meanG, sigmaG, weightG2, meanG2, sigmaG2);

//	  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 14);
//	  function->SetNpx(1000);
//	  Double_t parameters[14] = {constant, weight, mean, sigma1, sigma2, weightG, meanG, sigmaG, weightG2, meanG2, sigmaG2, weightG3, meanG3, sigmaG3};
//	  function->SetParameters(parameters);

//	  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 17);
//	  function->SetNpx(1000);
//	  Double_t parameters[17] = {constant, weight, mean, sigma1, sigma2, weightG, meanG, sigmaG, weightG2, meanG2, sigmaG2, weightG3, meanG3, sigmaG3, weightG4, meanG4, sigmaG4};
//	  function->SetParameters(parameters);

	if(cut_fit == true){
		for(int i = 1; i<=new_numberOfBins; i++){
			numberOfCurrentBins=i;
			Int_t observedValue = new_lightcurve->GetBinContent(i);
			Double_t startBin = new_lightcurve->GetBinLowEdge(i);
			Double_t widthBin = new_lightcurve->GetBinWidth(i);
			Double_t errorBin = new_lightcurve->GetBinError(i);
			if(observedValue!=0){
				fullBins++;
				Double_t expectedValue = (function->Integral(startBin, startBin+widthBin))/widthBin;
				chiSquare+= (TMath::Power(observedValue-expectedValue,2))/expectedValue;
				residual[i-1] = (observedValue-expectedValue)/errorBin;
				residual_squared[i-1] = residual[i-1]*residual[i-1];
		//			cout << "observedValue: " << observedValue << endl;
		//			cout << "expectedValue: " << expectedValue << endl;
		//			cout << "chiSquare: " << chiSquare << endl;
			}
		}
	}

	else{
		for(int i = 1; i<=numberOfBins; i++){
			numberOfCurrentBins=i;
			Int_t observedValue = timehisto->GetBinContent(i);
			Double_t startBin = timehisto->GetBinLowEdge(i);
			Double_t widthBin = timehisto->GetBinWidth(i);
			Double_t errorBin = timehisto->GetBinError(i);
			if(observedValue!=0){
				fullBins++;
				Double_t expectedValue = (function->Integral(startBin, startBin+widthBin))/widthBin;
				chiSquare+= (TMath::Power(observedValue-expectedValue,2))/expectedValue;
				residual[i-1] = (observedValue-expectedValue)/errorBin;
				residual_squared[i-1] = residual[i-1]*residual[i-1];
	//			cout << "observedValue: " << observedValue << endl;
	//			cout << "expectedValue: " << expectedValue << endl;
	//			cout << "chiSquare: " << chiSquare << endl;
			}
		}
	}


	cout << "Number of bins: " << numberOfCurrentBins << endl;
	cout << "Number of non-empty bins: " << fullBins << endl;
	Int_t parnum = 5;
//	Int_t parnum = 8;
//	Int_t parnum = 11;
//	Int_t parnum = 14;
//	Int_t parnum = 17;
	Int_t NDF = fullBins-parnum+1;
	cout << "NDF:" << NDF << endl;

	TF1* chisquare =  new TF1("chisquare", "ROOT::Math::chisquared_cdf(x*[0],[0])", 0, 10);
	chisquare->SetNpx(1000);
	chisquare->SetParameter(0, NDF);
//	TCanvas* chicancas =  new TCanvas();
//	chisquare->Draw();
	cout << "Eval: " << chisquare->Eval(chiSquare/NDF) << endl;
	pvalue = 1-(chisquare->Eval(chiSquare/NDF));

	return chiSquare/NDF;

}

Double_t PDFbin(Int_t numberOfBin, Double_t binValue, Double_t *par)
{

  Double_t Constant = par[0];
  Double_t Weight = par[1];
  Double_t Mean = par[2];
  Double_t Sigma = par[3];
  Double_t Sigma2 = par[4];
  Double_t WeightG = par[5];
  Double_t MeanG = par[6];
  Double_t SigmaG = par[7];
  Double_t WeightG2 = par[8];
  Double_t MeanG2 = par[9];
  Double_t SigmaG2 = par[10];
  Double_t WeightG3 = par[11];
  Double_t MeanG3 = par[12];
  Double_t SigmaG3 = par[13];
  Double_t WeightG4 = par[14];
  Double_t MeanG4 = par[15];
  Double_t SigmaG4 = par[16];


  Double_t PDF=0;

//  TF1* function = new TF1("function", One_gaus, mintime, maxtime, 4);
//  function->SetNpx(1000);
//  function->SetParameters(Constant, Weight, Mean, Sigma);

  TF1* function = new TF1("function", One_Asymgaus, mintime, maxtime, 5);
  function->SetNpx(1000);
  function->SetParameters(Constant, Weight, Mean, Sigma, Sigma2);

//  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 8);
//  function->SetNpx(1000);
//  function->SetParameters(Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG);

//  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 11);
//  function->SetNpx(1000);
//  function->SetParameters(Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG, WeightG2, MeanG2, SigmaG2);

//  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 14);
//   function->SetNpx(1000);
//   Double_t parameters[14] = {Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG, WeightG2, MeanG2, SigmaG2, WeightG3, MeanG3, SigmaG3};
//   function->SetParameters(parameters);

//    TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 17);
//     function->SetNpx(1000);
//     Double_t parameters[17] = {Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG, WeightG2, MeanG2, SigmaG2, WeightG3, MeanG3, SigmaG3,  WeightG4, MeanG4, SigmaG4};
//     function->SetParameters(parameters);

  Double_t startBin;
  Double_t widthBin;
  Double_t funValue;


  if(cut_fit == true){
		startBin = new_lightcurve->GetBinLowEdge(numberOfBin);
		widthBin = new_lightcurve->GetBinWidth(numberOfBin);
		funValue = function->Integral(startBin, startBin+widthBin)/widthBin;
  }

  else{
		startBin = timehisto->GetBinLowEdge(numberOfBin);
		widthBin = timehisto->GetBinWidth(numberOfBin);
		funValue = function->Integral(startBin, startBin+widthBin)/widthBin;
  }

  Double_t exp = TMath::Exp(-funValue);
  Double_t power = TMath::Power(funValue, binValue-(binValue/2));
  PDF = (exp*power*TMath::Power(funValue,(binValue/2)))/TMath::Factorial(binValue);

  return PDF;

}

Double_t PDFbin_test(Int_t numberOfBin)
{

  Int_t binValue = 0;
  if(cut_fit == true) binValue = new_lightcurve->GetBinContent(numberOfBin);
  else binValue = timehisto->GetBinContent(numberOfBin);

  cout << "binValue: " << binValue << endl;

  Double_t Constant = 90.;
  Double_t Weight = 180.;
  Double_t Mean = 4327.45;
  Double_t Sigma = 2214.17;
  Double_t Sigma2 = 1486.49;
  Double_t WeightG = 22.49;
  Double_t MeanG = 1937.54;
  Double_t SigmaG = 220.64;
  Double_t WeightG2 = 17.93;
  Double_t MeanG2 = 4334.91;
  Double_t SigmaG2 = 60.71;
  Double_t WeightG3 = 17.15;
  Double_t MeanG3 = 3122.90;
  Double_t SigmaG3 = 170.33;
  Double_t WeightG4 = 18.11;
  Double_t MeanG4 = 3124.96;
  Double_t SigmaG4 = 176.56;


  Double_t PDF=0.;

//  TF1* function = new TF1("function", One_gaus, mintime, maxtime, 4);
//  function->SetNpx(1000);
//  function->SetParameters(Constant, Weight, Mean, Sigma);


  TF1* function = new TF1("function", One_Asymgaus, mintime, maxtime, 5);
  function->SetNpx(1000);
  function->SetParameters(Constant, Weight, Mean, Sigma, Sigma2);


//  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 8);
//  function->SetNpx(1000);
//  function->SetParameters(Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG);


//  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 11);
//  function->SetNpx(1000);
//  function->SetParameters(Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG, WeightG2, MeanG2, SigmaG2);


//  TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 14);
//  function->SetNpx(1000);
//  Double_t parameters[14] = {Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG, WeightG2, MeanG2, SigmaG2, WeightG3, MeanG3, SigmaG3};
//  function->SetParameters(parameters);

//    TF1* function = new TF1("function", One_Asymgaus_more, mintime, maxtime, 17);
//    function->SetNpx(1000);
//    Double_t parameters[17] = {Constant, Weight, Mean, Sigma, Sigma2, WeightG, MeanG, SigmaG, WeightG2, MeanG2, SigmaG2, WeightG3, MeanG3, SigmaG3, WeightG4, MeanG4, SigmaG4};
//    function->SetParameters(parameters);

  Double_t startBin;
  Double_t widthBin;
  Double_t funValue;


  if(cut_fit == true){
		startBin = new_lightcurve->GetBinLowEdge(numberOfBin);
		widthBin = new_lightcurve->GetBinWidth(numberOfBin);
		funValue = function->Integral(startBin, startBin+widthBin)/widthBin;
  }

  else{
		startBin = timehisto->GetBinLowEdge(numberOfBin);
		widthBin = timehisto->GetBinWidth(numberOfBin);
		funValue = function->Integral(startBin, startBin+widthBin)/widthBin;
  }


  Double_t exp = TMath::Exp(-funValue);
  Double_t power = TMath::Power(funValue, binValue-(binValue/2));
  PDF = (exp*power*TMath::Power(funValue,(binValue/2)))/TMath::Factorial(binValue);

  cout << "binValue: " << binValue << endl;
  cout << "binValue/2: " << binValue/2 << endl;
  cout << "Integral between " << startBin << " and " << startBin+widthBin << endl;
  cout << "funValue: " << funValue << endl;
  cout << "exp: " << exp << endl;
  cout << "power: " << power << endl;
  cout << "PDF: " << PDF << endl;


  return PDF;

}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
 // to avoid warnings only
  gin[0]*=1;
  npar*=1;
  iflag*=1;


  Double_t logL = 0.;
  Double_t pretemp = 0;
  Double_t temp = 0.;
  Int_t bins;


  if(cut_fit == true) bins = new_numberOfBins;
  else bins = numberOfBins;
  for (Int_t i=1; i<=bins; i++)
  {

		  if(cut_fit == true){
			  if(new_lightcurve->GetBinContent(i) > 0){
				  pretemp = PDFbin(i,new_lightcurve->GetBinContent(i), par);
				  temp = TMath::Log(PDFbin(i,new_lightcurve->GetBinContent(i), par));
				  logL+=(-2.)*temp; //Final probability function ( addition of all events).
			  }
		  }
		  else{
			  if(timehisto->GetBinContent(i) > 0){
				  pretemp = PDFbin(i,timehisto->GetBinContent(i), par);
				  temp = TMath::Log(PDFbin(i,timehisto->GetBinContent(i), par));
				  logL+=(-2.)*temp; //Final probability function ( addition of all events).
			  }
		  }

//	  cout << "For bin " << i << " PDF is: "  << pretemp << endl;
//	  cout << "For bin " << i << " LogPDF is: "  << temp << endl;

  }

  f = logL;  //Final function.

  cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " "  << par[4] << " " << f << endl;
//  cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " "  << par[4] << " " << par[5] << " " << par[6] << " "  << par[7] << " " << f << endl;
//  cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " "  << par[4] << " " << par[5] << " " << par[6] << " "  << par[7] << " "<< par[8] << " " << par[9] << " "  << par[10] << f << endl;
//  cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " "  << par[4] << " " << par[5] << " " << par[6] << " "  << par[7] << " "<< par[8] << " " << par[9] << " "  << par[10] << " "<< par[11] << " " << par[12] << " "  << par[13] << f << endl;
//    cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " "  << par[4] << " " << par[5] << " " << par[6] << " "  << par[7] << " "<< par[8] << " " << par[9] << " "  << par[10] << " "<< par[11] << " " << par[12] << " "  << par[13] << " "<< par[14] << " " << par[15] << " "  << par[16]<< f << endl;


}

void fcn_test()
{
  Double_t logL = 0.;
  Double_t pretemp = 0;
  Double_t temp = 0.;
  Int_t bins;
  if(cut_fit == true) bins = new_numberOfBins;
  else bins = numberOfBins;
  for (Int_t i=1; i<=bins; i++)
  {
	  if(cut_fit == true){
		  if(new_lightcurve->GetBinContent(i) > 0){
		  pretemp = PDFbin_test(i);
		  temp = TMath::Log(PDFbin_test(i));
		  cout << "For bin " << i << " PDF is: "  << pretemp << endl;
		  cout << "For bin " << i << " LogPDF is: "  << temp << endl;
		  logL+=(-2.)*temp; //Final probability function ( addition of all events).
		  cout << "Current Likelihood value: " << logL << endl;
		  }
	  }

	  else{
	  		  if(timehisto->GetBinContent(i) > 0){
	  		  pretemp = PDFbin_test(i);
	  		  temp = TMath::Log(PDFbin_test(i));
	  		  cout << "For bin " << i << " PDF is: "  << pretemp << endl;
	  		  cout << "For bin " << i << " LogPDF is: "  << temp << endl;
	  		  logL+=(-2.)*temp; //Final probability function ( addition of all events).
	  		  cout << "Current Likelihood value: " << logL << endl;
	  		  }
	  	  }
  }

  cout << "Likelihood value: " << logL << endl;

}

void LC_LikeFit(){




	//Read Events and histogram. Use its array to create another histogram.
	//--Time histogram--
	TFile* HistoFile = TFile::Open("./Thesis_Histos/Histo_25ON_Kazuma_120.root", "READ");

		if (!HistoFile) {
				cout << "ERROR: file " << HistoFile << " not found..."	<< endl;
			}
		TH1D* lightcurve = dynamic_cast<TH1D*>(HistoFile->Get("ONevents"));
		if (!lightcurve) {
			cout << "ONevents not found" << endl;
		}

	numberOfBins = lightcurve->GetNbinsX();
	const Double_t* arrayBins = lightcurve->GetXaxis()->GetXbins()->GetArray();
	cout << "Number of bins in LC: " << numberOfBins << endl;


	//--List of events--
	ifstream in;
//	in.open("./FinalEventLists/Mrk421_2014flare_ONLikelihoodEvents_01-20TeV.txt");
	in.open("./Thesis_EventList/Event_Kazuma_120.txt");

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


	mintime = lightcurve->GetXaxis()->GetBinLowEdge(1);
	maxtime = lightcurve->GetXaxis()->GetBinUpEdge(numberOfBins);



	//--Cuts for the fit, if necessary

	if(cut_fit == true){
		Double_t Emin = 2600;
		Double_t Emax = 5000;

		cout << "I will do cuts! Emin = " << Emin << " Emax = " << Emax << endl;

		Int_t bin_Emin = lightcurve->FindBin(Emin);
		Int_t bin_Emax = lightcurve->FindBin(Emax);

		cout << "The first bin is " << bin_Emin << " and the last " << bin_Emax << endl;
		cout << "Array should start in " << lightcurve->GetBinLowEdge(bin_Emin) << " and finished in " << lightcurve->GetBinLowEdge(bin_Emax+1) << endl;



		new_numberOfBins = (bin_Emax-bin_Emin)+1;
		Int_t index_new = (bin_Emax-bin_Emin)+1;
		const Int_t index_new2 = index_new+1;
		Double_t new_arrayBins[index_new2];

		cout << "The new histo will have " << new_numberOfBins << " bins and an array of " << index_new << " components." << endl;

		for(Int_t i = 0; i<=index_new; i++){
			Int_t index = (bin_Emin-1)+1*i;
//			cout << "index: " << index << endl;
			Double_t value = *(arrayBins+index);
			new_arrayBins[i]= value;
//			cout << "array index: " << i << endl;
//			cout << "array value: " << new_arrayBins[i] << endl;
		}

		cout << "First array value: " << new_arrayBins[0] << " last " << new_arrayBins[index_new] << endl;


		new_lightcurve =  new TH1D("new_lightcurve", "new_lightcurve", new_numberOfBins, new_arrayBins);
		for(Int_t i=1; i<=new_numberOfBins; i++){
			new_lightcurve->SetBinContent(i, lightcurve->GetBinContent(bin_Emin+(i-1)));
//			cout << "To the bin: " << i << " we put the value: " << lightcurve->GetBinContent(bin_Emin+(i-1)) << endl;
		}

		mintime = new_lightcurve->GetXaxis()->GetBinLowEdge(1);
		maxtime = new_lightcurve->GetXaxis()->GetBinUpEdge(new_numberOfBins);

	}

	cout << "mintime: " << mintime << " maxtime: " << maxtime << endl;



	//Draw function (to deduce the initial value of parameters)

//	TF1* gaus = new TF1("gaus", One_gaus, mintime, maxtime,4);
//	gaus->SetNpx(2000);
//	gaus->SetParameters(24.73, 64.09, 3833.97, 1755.66);

	TF1* gaus = new TF1("gaus", One_Asymgaus, mintime, maxtime,5);
	gaus->SetParameters(0.5, 0.35, 4500., 2000., 1500.);
	gaus->SetNpx(5000);


//	TF1* gaus = new TF1("gaus", One_Asymgaus_more, mintime, maxtime,8);
//	gaus->SetNpx(2000);

//	TF1* gaus = new TF1("gaus", One_Asymgaus_more, mintime, maxtime,11);
//	gaus->SetNpx(2000);

//	TF1* gaus = new TF1("gaus", One_Asymgaus_more, mintime, maxtime,14);
//	gaus->SetNpx(2000);

//	TF1* gaus = new TF1("gaus", One_Asymgaus_more, mintime, maxtime, 17);
//	gaus->SetNpx(2000);

	//Function to save ;)
//	TF1* final_fit = new TF1("gaus", One_Asymgaus_more_Save, mintime, maxtime, 17);
//	final_fit->SetNpx(2000);


	timehisto = new TH1D("Lightcurve", "Lightcurve", numberOfBins, arrayBins);
	timehisto->SetStats(0);
	timehisto->GetXaxis()->SetTitle("Time(s)");
	timehisto->GetYaxis()->SetTitle("# events");
	for(int i = 0; i<Npoints; i++){
		timehisto->Fill(tp[i]);
	}

	residualhisto = new TH1D("Residuals", "Residuals", numberOfBins, arrayBins);
	timehisto->SetStats(0);

	residualsquarehisto = new TH1D("Residuals_Squared", "Residuals_Squared", numberOfBins, arrayBins);
	timehisto->SetStats(0);


	TCanvas* canvas = new TCanvas("canvas", "canvas", 2000,2000);
	TPad* upper =  new TPad("upper", "upper",0.005, 0.3525, 0.995, 0.995);
	TPad* lower =  new TPad("lower", "lower", 0.005, 0.005, 0.995, 0.3475);
	canvas->cd();
	upper->Draw();
	lower->Draw();
	upper->cd();
	timehisto->Draw("ehist");
	gaus->SetLineColor(2);
	gaus->Draw("same");
//	gaus->Draw("");
	if(cut_fit == true){
		new_lightcurve->SetLineColor(kBlack);
		new_lightcurve->SetLineWidth(2);
		new_lightcurve->Draw("same");
	}

	return;

	//------------------------------------Likelihood-------------------------------
	 myMinuit=new TMinuit(5); //(CHANGE HERE)
	 myMinuit->SetFCN(fcn);
	 myMinuit->SetPrintLevel(1);

	 Double_t arglist[10];
	 Int_t ierflg = 0;

	 //---Parameter definition--- (CHANGE HERE)


	 myMinuit->mnparm(0, "Constant", 40., 1., 0., 1000.,ierflg);
	 myMinuit->mnparm(1, "Weight", 70., 1., 0., 10000,ierflg);
	 myMinuit->mnparm(2, "Mean", 4500., 100., 0., maxtime,ierflg);
	 myMinuit->mnparm(3, "Sigma_Left", 2100., 100., 0., 10000.,ierflg);
	 myMinuit->mnparm(4, "Sigma_Right", 1300., 100., 0., 10000.,ierflg);
//	 myMinuit->mnparm(5, "WeightG", 4190.89, 1., 0., 10000.,ierflg);
//	 myMinuit->mnparm(6, "MeanG", 48388.4, 50., 0., maxtime*10.,ierflg);
//	 myMinuit->mnparm(7, "SigmaG", 931.44, 10., 0., 10000.,ierflg);
//	 myMinuit->mnparm(8, "WeightG2", 63.79, 1., 0, 10000.,ierflg);
//	 myMinuit->mnparm(9, "MeanG2", 5066.46, 50., 0, maxtime,ierflg);
//	 myMinuit->mnparm(10, "SigmaG2", 1172.27, 10., 0., 10000.,ierflg);
//	 myMinuit->mnparm(11, "WeightG3", 35.91, 1., 0, 10000.,ierflg);
//	 myMinuit->mnparm(12, "MeanG3", 3129.90, 50., 0, maxtime,ierflg);
//	 myMinuit->mnparm(13, "SigmaG3", 218.13, 10., 0., 10000.,ierflg);
//	 myMinuit->mnparm(14, "WeightG4", 36.22, 1., 0, 10000.,ierflg);
//	 myMinuit->mnparm(15, "MeanG4", 4457.9, 50., 0, maxtime,ierflg);
//	 myMinuit->mnparm(16, "SigmaG4", 362.08, 10., 0., 10000.,ierflg);


	 myMinuit->SetErrorDef(1);
//	 myMinuit->FixParameter(0);
//	 myMinuit->FixParameter(1);
//	 myMinuit->FixParameter(2);
//	 myMinuit->FixParameter(3);
//	 myMinuit->FixParameter(4);
//	 myMinuit->FixParameter(5);
//	 myMinuit->FixParameter(6);
//	 myMinuit->FixParameter(7);
//	 myMinuit->FixParameter(8);
//	 myMinuit->FixParameter(9);
//	 myMinuit->FixParameter(10);
//	 myMinuit->FixParameter(11);
//	 myMinuit->FixParameter(12);
//	 myMinuit->FixParameter(13);
//	 myMinuit->FixParameter(14);
//	 myMinuit->FixParameter(15);
//	 myMinuit->FixParameter(16);

    minuit = myMinuit;

	//Remember to indicate the arguments for every Minuit command! (Check documentation)
	 arglist[0] = 0;
	 minuit->mnexcm("SET STRategy", arglist, 1, ierflg); //the 1 o 2 is the number of arguments used from the arglist
	 cout << "ierflg after SETSTRatregy: " << ierflg << endl;
	 arglist[0] = 3000;
//	 arglist[1] = 0.1;
	 minuit->mnexcm("SIMplex", arglist, 2, ierflg);
	 cout << "ierflg after SIMplex: " << ierflg << endl;

//	 myMinuit->mnexcm("MINIMIZE", arglist ,2,ierflg);
//	 cout << "ierflg after MINIMIZE: " << ierflg << endl;
//
//	 myMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
//	 cout << "ierflg after MIGRAD: " << ierflg << endl;



	 //Recover results: Plot them and test statistics
	 Double_t const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res, weigthG2_res, meanG2_res, sigmaG2_res, weigthG3_res, meanG3_res, sigmaG3_res,  weigthG4_res, meanG4_res, sigmaG4_res;
	 Double_t error_const_res, error_weigth_res, error_mean_res, error_sigma_res, error_sigma2_res, error_weigthG_res, error_meanG_res, error_sigmaG_res, error_weigthG2_res, error_meanG2_res, error_sigmaG2_res,error_weigthG3_res, error_meanG3_res, error_sigmaG3_res, error_weigthG4_res, error_meanG4_res, error_sigmaG4_res;


	 myMinuit->GetParameter(0, const_res, error_const_res);
	 myMinuit->GetParameter(1, weigth_res, error_weigth_res);
	 myMinuit->GetParameter(2, mean_res, error_mean_res);
	 myMinuit->GetParameter(3, sigma_res, error_sigma_res);
	 myMinuit->GetParameter(4, sigma2_res, error_sigma2_res);
//	 myMinuit->GetParameter(5, weigthG_res, error_weigthG_res);
//	 myMinuit->GetParameter(6, meanG_res, error_meanG_res);
//	 myMinuit->GetParameter(7, sigmaG_res, error_sigmaG_res);
//	 myMinuit->GetParameter(8, weigthG2_res, error_weigthG2_res);
//	 myMinuit->GetParameter(9, meanG2_res, error_meanG2_res);
//	 myMinuit->GetParameter(10, sigmaG2_res, error_sigmaG2_res);
//	 myMinuit->GetParameter(11, weigthG3_res, error_weigthG3_res);
//	 myMinuit->GetParameter(12, meanG3_res, error_meanG3_res);
//	 myMinuit->GetParameter(13, sigmaG3_res, error_sigmaG3_res);
//	 myMinuit->GetParameter(14, weigthG4_res, error_weigthG4_res);
//	 myMinuit->GetParameter(15, meanG4_res, error_meanG4_res);
//	 myMinuit->GetParameter(16, sigmaG4_res, error_sigmaG4_res);

	 gaus->SetParameters(const_res, weigth_res, mean_res, sigma_res, sigma2_res);
//	 gaus->SetParameters(const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res);
//	 gaus->SetParameters(const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res, weigthG2_res, meanG2_res, sigmaG2_res);
//	 Double_t parameters[14] = {const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res, weigthG2_res, meanG2_res, sigmaG2_res, weigthG3_res, meanG3_res, sigmaG3_res};
//	 gaus->SetParameters(parameters);
//	 Double_t parameters[17] = {const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res, weigthG2_res, meanG2_res, sigmaG2_res, weigthG3_res, meanG3_res, sigmaG3_res, weigthG4_res, meanG4_res, sigmaG4_res};
//	 gaus->SetParameters(parameters);
//	 final_fit->SetParameters(parameters);
//
	 gaus->Draw("same");

	 Double_t chi = test_statistics(const_res, weigth_res, mean_res, sigma_res, sigma2_res);
//	 Double_t chi = test_statistics(const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res);
//	 Double_t chi = test_statistics(const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res,weigthG2_res, meanG2_res, sigmaG2_res);
//	 Double_t chi = test_statistics(const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res,weigthG2_res, meanG2_res, sigmaG2_res, weigthG3_res, meanG3_res, sigmaG3_res);
//	 Double_t chi = test_statistics(const_res, weigth_res, mean_res, sigma_res, sigma2_res, weigthG_res, meanG_res, sigmaG_res,weigthG2_res, meanG2_res, sigmaG2_res, weigthG3_res, meanG3_res, sigmaG3_res, weigthG4_res, meanG4_res, sigmaG4_res);


	 cout << "Chi2/NDF " << chi << endl;
	 cout << "p-value: " << pvalue << endl;


	 //Here the residuals historgrams
	 for(int i = 1; i<=numberOfBins; i++){
		 residualhisto->SetBinContent(i, residual[i-1]);
		 residualsquarehisto->SetBinContent(i, residual_squared[i-1]);
	 }

	 lower->cd();
	 residualhisto->SetStats(0);
	 residualhisto->GetYaxis()->SetTitle("Simga");
	 residualhisto->GetXaxis()->SetTitle("Time(s)");
	 residualhisto->SetLineColor(4);
	 residualhisto->SetLineWidth(2);
	 residualhisto->Draw("EB");


	 upper->cd();
	 TLatex *txt = new TLatex();
//	 txt->DrawLatex(8000, 350, Form("C_{a} = %.2f #pm %.2f",const_res, error_const_res));
//	 txt->DrawLatex(8000, 330, Form("w_{a} = %.2f #pm %.2f",weigth_res, error_weigth_res));
//	 txt->DrawLatex(8000, 310, Form("#mu_{a} = %.2f #pm %.2f",mean_res,error_mean_res));
//	 txt->DrawLatex(8000, 290, Form("#sigma_{left} = %.2f #pm %.2f",sigma_res, error_sigma_res));
//	 txt->DrawLatex(8000, 270, Form("#sigma_{right} = %.2f #pm %.2f",sigma2_res, error_sigma2_res));

	 txt->DrawLatex(8000, 122, Form("C_{a} = %.2f #pm %.2f",const_res, error_const_res));
	 txt->DrawLatex(8000, 114, Form("w_{a} = %.2f #pm %.2f",weigth_res, error_weigth_res));
	 txt->DrawLatex(8000, 106, Form("#mu_{a} = %.2f #pm %.2f",mean_res,error_mean_res));
	 txt->DrawLatex(8000, 98, Form("#sigma_{left} = %.2f #pm %.2f",sigma_res, error_sigma_res));
	 txt->DrawLatex(8000, 90, Form("#sigma_{right} = %.2f #pm %.2f",sigma2_res, error_sigma2_res));
//	 txt->DrawLatex(8000, 130, Form("w_{s} = %.2f #pm %.2f",weigthG_res, error_weigthG_res));
//	 txt->DrawLatex(8000, 122, Form("#mu_{s} = %.2f #pm %.2f",meanG_res, error_meanG_res));
//	 txt->DrawLatex(8000, 114, Form("#sigma_{s} = %.2f #pm %.2f",sigmaG_res, error_sigmaG_res));
//	 txt->DrawLatex(8000, 106, Form("w2_{s} = %.2f #pm %.2f",weigthG2_res, error_weigthG2_res));
//	 txt->DrawLatex(8000, 98, Form("#mu2_{s} = %.2f #pm %.2f",meanG2_res, error_meanG2_res));
//	 txt->DrawLatex(8000, 90, Form("#sigma2_{s} = %.2f #pm %.2f",sigmaG2_res, error_sigmaG2_res));
//	 txt->DrawLatex(8000, 82, Form("w3_{s} = %.2f #pm %.2f",weigthG3_res, error_weigthG3_res));
//	 txt->DrawLatex(8000, 74, Form("#mu3_{s} = %.2f #pm %.2f",meanG3_res, error_meanG3_res));
//	 txt->DrawLatex(8000, 66, Form("#sigma3_{s} = %.2f #pm %.2f",sigmaG3_res, error_sigmaG3_res));
//	 txt->DrawLatex(8000, 58, Form("w4_{s} = %.2f #pm %.2f",weigthG4_res, error_weigthG4_res));
//	 txt->DrawLatex(8000, 50, Form("#mu4_{s} = %.2f #pm %.2f",meanG4_res, error_meanG4_res));
//	 txt->DrawLatex(8000, 42, Form("#sigma4_{s} = %.2f #pm %.2f",sigmaG4_res, error_sigmaG4_res));
	 txt->DrawLatex(8000, 250, Form("p-value = %.9f",pvalue));
	 txt->DrawLatex(8000, 230, Form("#chi^{2}/NDF = %.2f",chi));
	 txt->SetNDC(kTRUE);
	 txt->SetTextAlign(13);
	 txt->SetTextSize(10);
	 txt->Draw("same");


//	 TCanvas* residualsquearedhisto = new TCanvas();
//	 residualsquarehisto->SetLineColor(2);
//	 residualsquarehisto->SetLineWidth(2);
//	 residualsquarehisto->Draw("EB");




	 //Create a root file to save the final function and print the cov matrix

//	 TFile* fit_outputfile = new TFile("Fit_Kazuma_120_1min.root","recreate");
////	 final_fit->SetName("Fit_function");
////	 final_fit->Write();
//	 canvas->SetName("Fit_canvas");
//	 canvas->Write();
//	 fit_outputfile->Close();


//	 double fmin,fedm;
//	 double errdef;
//	 int    nparv,nparx;
//	 int    fstat;
//
//
//	  myMinuit->mnstat(fmin,fedm,errdef,nparv,nparx,fstat);
//
//	  printf("\n\n");
//	  printf("Results of MINUIT minimisationa\n");
//	  printf("-------------------------------------\n\n");
//
//	  printf(" Minimal function value:              %8.3lf  \n",fmin);
//	  printf(" Estimated difference to true minimum: %11.3le \n",fedm);
//	  printf(" Number of parameters:         %3i     \n",nparv);
//	  printf(" Error definition (Fmin + Delta):      %8.3lf  \n",errdef);
//	   if (fstat==3) {
//	         printf(" Exact covariance matrix.\n");
//	   }
//	   else {
//	        printf(" No/error with covariance matrix.\n");
//	        printf(" Error code: %i3\n",fstat);
//	    }
//	  printf("\n");
//
//	   cout << "Covariance Matrix: " << endl;
//	   myMinuit->mnexcm("SHOw COVariance",arglist,0,ierflg);
//	   cout << "Eigenvalues: " << endl;
//	   myMinuit->mnexcm("SHOw  EIGenvalues",arglist,0,ierflg);
//
//	   cout << "Covariance Matrix only nuisance parameters: " << endl;
//	   myMinuit->FixParameter(0);
//	   myMinuit->mnexcm("SHOw COVariance",arglist,0,ierflg);
//	   cout << "Eigenvalues: " << endl;
//	   myMinuit->mnexcm("SHOw  EIGenvalues",arglist,0,ierflg);

}


