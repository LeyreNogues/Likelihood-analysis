/*
 * LikeManel.C
 *
 *  Created on: Jul 4, 2016
 *      Author: lnogues
 */


/* Parts of the code:
1. Read of events: Cuts and preparation of event list.
2. Collection area. Read and prepare function.
3. Migration matrix: Obtaining a value from two energies and get bias and resolution as function of energy.
4. Spectrum: Intrinsic assumption and MC for parameter fit.
5. Intrinsic time: Compute delay and extract intrinsic time of events.
6. Lightcurve at the source: Fit events and let some parameters free.
7. Normalization
 */

#include <TF1.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TRandom.h>
#include <TStyle.h>
#include<iostream>
using namespace std;


const Int_t numberOfEvents = 1491;

//Options
Int_t QGorder= 1; //Linear or Quadratic
Bool_t ResolutionInc = true; //Include or not resolution

//Variables

Double_t E[numberOfEvents],t[numberOfEvents];
Double_t tmin=0;
Double_t tmax=2750.;
Double_t Emin=150.;
Double_t Emax=12000;
Double_t res=0.22;
Double_t Like_res;

//Spectrum variables
Double_t signalIndex=-2.4;
Double_t bkgIndex=-2.7;
Double_t specIntegralSignal=(1./(signalIndex+1))*(TMath::Power(Emax,signalIndex+1)-TMath::Power(Emin,signalIndex+1));
Double_t specIntegralbkg=(1./(bkgIndex+1))*(TMath::Power(Emax,bkgIndex+1)-TMath::Power(Emin,bkgIndex+1));

//LC inputs
Double_t gausMax=2004.9;
Double_t gausWidth=221.51;
Double_t SBratio=0.369;

//QG inputs
Double_t z=0.034;
Double_t H0=70.;
Double_t PlanckScale=1.2*TMath::Power(10,19);
Double_t QGpar=20.;
//Double_t tau=1.5;
Double_t PcTom=3.085*TMath::Power(10,19);
Double_t QGfixedPart=z*PcTom/H0; //(ETA)

Double_t ResFormula(Double_t* x, Double_t* par){


	//par[0] -> QGpar
	//par[1] -> ratio
	//par[2] -> time of event
	//par[3] -> Erec of event
	//par[4] -> Gaus Maximum
	//par[5] -> Gaus Sigma
	//par[6] -> Res*Erec

	Double_t baseline = par[1]*TMath::Power(x[0],-2.7);

	Double_t flare;
	Double_t QGDelay;
	Double_t final;
	if(QGorder==1){
		QGDelay= QGfixedPart*par[0]/PlanckScale;
		flare = (1-par[1])*50.*TMath::Power(x[0],-2.4)*TMath::Exp(-(TMath::Power((par[2]-QGDelay*x[0]-par[4]),2))/(2*TMath::Power(par[5],2)))/(par[5]*TMath::Sqrt(2.*TMath::Pi()));
	}
	if(QGorder==2){
		QGDelay= QGfixedPart*par[0]/PlanckScale/1000;
		flare = (1-par[1])*50.*TMath::Power(x[0],-2.4)*TMath::Exp(-(TMath::Power((par[2]-QGDelay*x[0]*x[0]-par[4]),2))/(2*TMath::Power(par[5],2)))/(par[5]*TMath::Sqrt(2.*TMath::Pi()));
	}

	Double_t resolution;
	if(x[0]>par[3]-4*par[6] && x[0]<par[3]+4*par[6]){
		resolution =TMath::Exp(-0.5*(TMath::Power((x[0]-par[3]),2)/(TMath::Power(par[6],2))))/(par[6]*TMath::Sqrt(2.*TMath::Pi()));
	}
	else{
		resolution=0.;
	}

	final = (baseline+flare)*resolution;

	return final;
}


void LoopManelLike(char *fname)
{
	//Read time and energy of events

  ifstream in;
  in.open(fname);
  double v1,v2;
  int Npoints = 0;
  while(1)
	{
	  in >> v1 >> v2;
	  if (!in.good()) break;
	  t[Npoints] = v1;
	  E[Npoints] = v2;
	  Npoints++;
	}

//  cout << "Number of events " << Npoints << endl;
//  cout << t[0] << " " << E[0] << endl;
//  cout << t[Npoints-1] << " " << E[Npoints-1] << endl;

  //Likelihood initial value
  //  Like_res = Likelihoodtest();
    if(ResolutionInc == true){
  	  Like_res = Likelihoodtest_Res();
    }

    if(ResolutionInc == false){
    	  Like_res = Likelihoodtest();
      }

    cout << "Res Likelihood computed: " << Like_res << endl;


  TMinuit *gMinuit = new TMinuit(4);
	gMinuit->SetFCN(fcn);
	gMinuit->SetPrintLevel(-1);

	double arglist[10];
	int ierflg = 0;

	//Parameter definition

	//Alpha case
	Double_t QGdown=-QGpar*100;
	Double_t QGup=QGpar*100;
	gMinuit->mnparm(0, "QGpar", QGpar, 5., QGdown, QGup,ierflg);


	//tau case
//	Double_t taudown=-tau*10;
//	Double_t tauGup=tau*10;
//	gMinuit->mnparm(0, "tau", tau, 0.1, taudown, tauGup,ierflg);

	Double_t gausMaxdown=gausMax-3*gausWidth;
	Double_t gausMaxup=gausMax+3*gausWidth;
	gMinuit->mnparm(1, "gausMax", gausMax, 100., gausMaxdown, gausMaxup,ierflg);

	Double_t gausWidthdown=gausWidth/5.;
	Double_t gausWidthup=gausWidth*5.;
	gMinuit->mnparm(2, "gausWidth", gausWidth, 10., gausWidthdown, gausWidthup,ierflg);


	Double_t SBratiodown=0.01;
	Double_t SBratioup=1.;
	gMinuit->mnparm(3, "SBratio", SBratio, 0.5, SBratiodown, SBratioup,ierflg);



	//Fix parameters
//	  gMinuit->FixParameter(0);
//	  gMinuit->FixParameter(1);
//	  gMinuit->FixParameter(2);
//	  gMinuit->FixParameter(3);



	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
//	gMinuit->mnexcm("MINOS", arglist ,2,ierflg);


  //Collect value of parameters
  Double_t QG0, QG, QGErr;
  Double_t Max, MaxErr;
  Double_t With, WithErr;
  Double_t Ratio, RatioErr;

  gMinuit->GetParameter(0, QG0, QGErr);
  gMinuit->GetParameter(1, Max, MaxErr);
  gMinuit->GetParameter(2, With, WithErr);
  gMinuit->GetParameter(3, Ratio, RatioErr);
  cout << "QG = " << QG0 << " Error = " << QGErr << endl;
  cout << "Max = " << Max << " Error = " << MaxErr << endl;
  cout << "With = " << With << " Error = " << WithErr << endl;
  cout << "Ratio = " << Ratio << " Error = " << RatioErr << endl;
  Double_t fmin,edm,errdef,amin, bmin, cmin;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(fmin,edm,errdef,nvpar,nparx,icstat);
  cout << "FCN min value = " << fmin << endl;

  cout << "I fix the QG parameter and start the loop" << endl;

  cout << "***************************************************************" << endl;
  cout << "***************          LOOP RIGHT          ******************" << endl;
  cout << "***************************************************************" << endl;

  TGraph* GraphR = new TGraph();
  Double_t DeltaLklR[100], QGR[100], MaxR[100],WithR[100],RatioR[100];
  //Middle point of the graph
  DeltaLklR[0]= 0;
  QGR[0]= QG0;
  MaxR[0]= Max;
  WithR[0]= With;
  RatioR[0]= Ratio;

  GraphR->SetPoint(0, QGR[0], DeltaLklR[0]);
  Int_t i1 = 1 ;
  bmin = fmin;
  QG = QG0;

  while ( TMath::Abs(bmin-fmin) <= 2.)
  {
	  QG+=1.;
	  gMinuit->mnparm(0, "QG" , QG, 0.1, QGdown, QGup,ierflg);
	  gMinuit->FixParameter(0);
	  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	  gMinuit->mnstat(bmin,edm,errdef,nvpar,nparx,icstat);
	  gMinuit->GetParameter(0, QGR[i1], QGErr);
	  gMinuit->GetParameter(1, MaxR[i1], MaxErr);
	  gMinuit->GetParameter(2, WithR[i1], WithErr);
	  gMinuit->GetParameter(3, RatioR[i1] , RatioErr);
	  DeltaLklR[i1] = TMath::Abs(bmin-fmin);
	  GraphR->SetPoint(i1-1, QGR[i1], DeltaLklR[i1]);

	  cout << "***************************************************************" << endl;
	  cout << "***************           RESULTS            ******************" << endl;
	  cout << "***************************************************************" << endl;

	  cout << "Valor de fcn3 = " << bmin << endl;
	  cout << "DeltaL = " << DeltaLklR[i1] << endl;
	  cout << "QG = " << QGR[i1] << endl;
	  cout << "Max = " << MaxR[i1] << endl;
	  cout << "Width = " << WithR[i1] << endl;
	  cout << "Ratio = " << RatioR[i1] << endl;
	  i1++;

 }
  i1-=1;


cout << "***************************************************************" << endl;
cout << "***************          LOOP LEFT           ******************" << endl;
cout << "***************************************************************" << endl;

TGraph* GraphL = new TGraph();
Double_t DeltaLklL[100], QGL[100], MaxL[100],WithL[100],RatioL[100];
Int_t i2 = 0 ;
cmin = fmin;
QG = QG0;
while ( TMath::Abs(cmin-fmin) <= 2.)
  {
	QG-=1.;
	gMinuit->Release(4);
	gMinuit->mnparm(0, "QG" , QG, 0.1, QGdown, QGup,ierflg);
	gMinuit->FixParameter(0);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnstat(cmin,edm,errdef,nvpar,nparx,icstat);
	gMinuit->GetParameter(0, QGL[i2], QGErr);
	gMinuit->GetParameter(1, MaxL[i2], MaxErr);
	gMinuit->GetParameter(2, WithL[i2], WithErr);
	gMinuit->GetParameter(3, RatioL[i2] , RatioErr);
	DeltaLklL[i2] = TMath::Abs(cmin-fmin);
	GraphL->SetPoint(i2,QGL[i2], DeltaLklL[i2]);

	cout << "***************************************************************" << endl;
	cout << "***************           RESULTS            ******************" << endl;
	cout << "***************************************************************" << endl;

	cout << "Valor de fcn3 = " << cmin << endl;
	cout << "DeltaL = " << DeltaLklL[i2] << endl;
	cout << "QG = " << QGL[i2] << endl;
	cout << "Max = " << MaxL[i2] << endl;
	cout << "Width = " << WithL[i2] << endl;
	cout << "Ratio = " << RatioL[i2] << endl;
	i2++;

  }

i2-=1;

TCanvas* LikeSides = new TCanvas();
LikeSides->Divide(1,2);
LikeSides->cd(1);
GraphR->Draw("*APL");
LikeSides->cd(2);
GraphL->Draw("*APL");


//Fill a vector with the right and left vertor
Double_t DeltaLkl[100],QGcomp[100], Maxcomp[100],Witdhcomp[100],Ratiocomp[100];
TGraph* Graph = new TGraph();
Int_t i3=0;
while(i3<=i2)
{
	 DeltaLkl[i3] = DeltaLklL[i2-i3];
	 QGcomp[i3] = QGL[i2-i3];
	 Maxcomp[i3] = MaxL[i2-i3];
	 Witdhcomp[i3] = WithL[i2-i3];
	 Ratiocomp[i3] = RatioL[i2-i3];
	 i3++;
}

Int_t i4=0;
while(i4<=i1+1)
{
	 DeltaLkl[i3+i4] = DeltaLklR[i4];
	 QGcomp[i3+i4] = QGR[i4];
	 Maxcomp[i3+i4] = MaxR[i4];
	 Witdhcomp[i3+i4] = WithR[i4];
	 Ratiocomp[i3+i4] = RatioR[i4];
	 i4++;
}
Double_t totalPoints=i3+i4-1;


//Create a file with results
FILE *fout;
fout = fopen("results_data_lin_Manel.txt", "wb");
for ( Int_t b=0; b<totalPoints; b++)
 {
   Graph->SetPoint(b, QGcomp[b]*(TMath::Sqrt(1/(8*TMath::Pi()))), DeltaLkl[b]);
   fprintf(fout,"%0.8lf    %0.8lf    %0.8lf    %0.8lf    \n", QGcomp[b], Maxcomp[b], Witdhcomp[b], Ratiocomp[b]);
 }


//Look for the error in the parameter
/*Double_t L_value, p;
Double_t sigmaLik=2;
p=QGcomp[0];
L_value = Graph->Eval(p);
Double_t step=0.002; //Steps for he evaluation of the function
Double_t e=0.05;
Double_t err[2]; //Left and right errors for the parameter
Int_t m=0;
while (L_value<6.)
{
  if( L_value>(sigmaLik-e) && L_value<(sigmaLik+e))
	{
	  cout << "Limit found for parameter  " << p << " at L_value  " << L_value << endl;
	  err[m] = p;
	  m+=1;
	  if(m == 2)break;
	  p = QG0;
	  L_value = Graph->Eval(p);
	  continue;
	}
//      cout << " Parameter " << p << " L_value " << L_value << endl;
  p+=step;
  L_value = Graph->Eval(p);
}*/


//Results
TCanvas* canvas = new TCanvas();
Graph->GetXaxis()->SetTitle("Mp/Mqg");
Graph->GetYaxis()->SetTitle("#Delta L");
Graph->Draw("*APL");

TString rootoutfile("Graph_data_linear");
TFile* outputfile = new TFile(rootoutfile, "recreate");
canvas->Write();
outputfile->Close();

//cout << "Mp/Mqf = " << QG0 << " - " << QG0-err[0] << " + " << err[1]-QG0 << endl;





}

double PDFevent(double tr, double Er, double *par)
{
	//  double tauParameter = par[0];
	  Double_t QGParameter = par[0];
	  Double_t MaxParameter = par[1];
	  Double_t WidthParameter= par[2];
	  Double_t RatioParameter = par[3];

	  Double_t QGDelay = 0;
	if(QGorder==1){ //Linear
		QGDelay= QGfixedPart*QGParameter*Er/PlanckScale;
	}

	if(QGorder==2){ //Quadratic
		QGDelay= QGfixedPart*QGParameter*TMath::Power(Er,2)/PlanckScale/1000;

	}

	  Double_t dt= tr-QGDelay;

	  Double_t PDF= RatioParameter*TMath::Power(Er,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(Er,signalIndex)*TMath::Exp(-(TMath::Power((dt-MaxParameter),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));

	  //Double_t PDF= RatioParameter*TMath::Power(E,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(E,signalIndex)*TMath::Exp(-(TMath::Power((t-tauParameter*E-MaxParameter),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));

	  //Normalization of the PDF

	  PDF=PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);
	  return PDF;

}

double PDFevent_Res(double tr, double Er, double *par)
{
	Double_t QGParameter = par[0];
	Double_t MaxParameter = par[1];
	Double_t WidthParameter= par[2];
	Double_t RatioParameter = par[3];

	Double_t sigma=Er*res;
	TF1* myresol = new TF1("myresol", ResFormula, Er-4*sigma, Er+4*sigma, 7);
	myresol->SetNpx(100000);
	myresol->SetParameters(QGParameter, RatioParameter, tr, Er, MaxParameter, WidthParameter, res*Er);
	myresol->SetLineColor(2);
	Double_t PDF= myresol->Integral(Er-4*sigma, Er+4*sigma);
	//  Double_t PDF= myresol->Integral(Emin/2, Er+4*sigma);
	myresol->Delete();
	PDF=PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);
	return PDF;
}

double PDFtest(double tr, double Er)
{
	Double_t QGParameter = QGpar;
	Double_t MaxParameter = gausMax;
	Double_t WidthParameter= gausWidth;
	Double_t RatioParameter = SBratio;

	Double_t QGDelay = 0;
	if(QGorder==1){ //Linear
		QGDelay= QGfixedPart*QGParameter*Er/PlanckScale;
	}

	if(QGorder==2){ //Quadratic
		QGDelay= QGfixedPart*QGParameter*TMath::Power(Er,2)/PlanckScale/1000;

	}

	Double_t dt= tr-QGDelay;
	Double_t PDF= RatioParameter*TMath::Power(Er,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(Er,signalIndex)*TMath::Exp(-(TMath::Power((dt-MaxParameter),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));


	//Normalization of the PDF
	PDF= PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);
	return PDF;
}

double PDFtestRes(double tr, double Er)
{
	Double_t QGParameter = QGpar;
	Double_t MaxParameter = gausMax;
	Double_t WidthParameter= gausWidth;
	Double_t RatioParameter = SBratio;

	Double_t sigma=Er*res;
	TF1* myresol = new TF1("myresol", ResFormula, Er-4*sigma, Er+4*sigma, 7);
	myresol->SetNpx(100000);
	myresol->SetParameters(QGParameter, RatioParameter, tr, Er, MaxParameter, WidthParameter, res*Er);
	myresol->SetLineColor(2);
	Double_t PDF= myresol->Integral(Er-4*sigma, Er+4*sigma);
	myresol->Delete();
	PDF=PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);
	return PDF;

}

double PDFtesttau(double t, double E)
{

//  double QGParameter = QGpar;
	double tauParameter = 0;
	double MaxParameter = gausMax;
	double WidthParameter= gausWidth;
	double RatioParameter = SBratio;

	Double_t PDF= RatioParameter*TMath::Power(E,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(E,signalIndex)*TMath::Exp(-(TMath::Power((t-(tauParameter*E)-MaxParameter),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));
  //Normalization of the PDF

  PDF= PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);

	  return PDF;


}


double PDFIntegral(){

	double MaxParameter2 = gausMax;//0
	double WidthParameter2= gausWidth;//1
	double RatioParameter2 = SBratio;//2

	TF2 *f2 = new TF2("f","[2]*TMath::Power(y,-2.7)+(1-[2])*50.*TMath::Power(y,-2.4)*TMath::Exp(-(TMath::Power((x-[0]),2))/(2*TMath::Power([1],2)))*1/([1]*TMath::Sqrt(2.*TMath::Pi()))",tmin,tmax,Emin,Emax);
	f2->SetParameters(MaxParameter2,WidthParameter2,RatioParameter2);
	double norm = f2->Integral(tmin,tmax,Emin,Emax);
	double Manelnorm =RatioParameter2*(tmax-tmin)*specIntegralbkg+(1-RatioParameter2)*50.*specIntegralSignal;
	f2->Delete();

	cout << "Integral norm: " << norm << endl;
	cout << "Manel norm: "  << Manelnorm << endl;
	return norm;

}

Double_t Likelihoodtest(){
	Double_t suma=0;
	  for(int i=0; i<numberOfEvents; i++){
		  suma+=(-2)*TMath::Log(PDFtest(t[i], E[i]));
	  }
	  return suma;
}

Double_t Likelihoodtest_Res(){
	Double_t suma2=0;
	  for(int i=0; i<numberOfEvents; i++){
		  suma2+=(-2)*TMath::Log(PDFtestRes(t[i], E[i]));
	  }
	  return suma2;
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{

	Double_t logL = 0.;
	Double_t temp = 0.;
	for (Int_t i=0; i<numberOfEvents; i++)
	{
		if(ResolutionInc == false){
			temp = TMath::Log(PDFevent(t[i],E[i],par));
		}
		if(ResolutionInc == true){
			temp = TMath::Log(PDFevent_Res(t[i],E[i],par));
		}
		logL+=(-2)*temp; //Final probability function ( addition of all events).
	}

	Double_t CHI0 = Like_res;
	f = logL-CHI0;  //Final function.

//	cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << f << endl;


}




