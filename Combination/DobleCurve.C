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

//Variables
const Int_t numberOfEvents_flare1 = 1491;
const Int_t numberOfEvents_flare2 = 1491;
const Int_t numberOfEvents = numberOfEvents_flare1+numberOfEvents_flare2;
Double_t E[numberOfEvents],t[numberOfEvents], flare[numberOfEvents];
Double_t E1[numberOfEvents_flare1],t1[numberOfEvents_flare1];
Double_t E2[numberOfEvents_flare2],t2[numberOfEvents_flare2];
Double_t tmin1=0;
Double_t tmax1=2750.;
Double_t Emin1=150.;
Double_t Emax1=12000;
Double_t tmin2=2750;
Double_t tmax2=5500.;
Double_t Emin2=150.;
Double_t Emax2=12000;

//Spectrum variables
Double_t signalIndex=-2.4;
Double_t bkgIndex=-2.7;
Double_t specIntegralSignal=(1./(signalIndex+1))*(TMath::Power(Emax1,signalIndex+1)-TMath::Power(Emin1,signalIndex+1));
Double_t specIntegralbkg=(1./(bkgIndex+1))*(TMath::Power(Emax1,bkgIndex+1)-TMath::Power(Emin1,bkgIndex+1));

//LC inputs
Double_t gausMax1=2004.9;
Double_t gausMax2=4754.9;
Double_t gausWidth=221.51;
Double_t SBratio=0.369;

//QG inputs
Double_t z=0.034;
Double_t H0=70.;
Double_t PlanckScale=1.2*TMath::Power(10,19);
Double_t QGpar=20.;
Double_t PcTom=3.085*TMath::Power(10,19);
Double_t QGfixedPart=z*PcTom/H0; //(ETA)

void DobleCurve(char *fname)
{
	//Read time and energy of events (2 flares together)

  ifstream in;
  in.open(fname);
  double v1,v2,v3;
  int Npoints = 0;
  while(1)
	{
	  in >> v1 >> v2 >> v3;
	  if (!in.good()) break;
	  flare[Npoints] = v1;
	  t[Npoints] = v2;
	  E[Npoints] = v3;
	  Npoints++;
	}

  cout << "Number of events " << Npoints << endl;
  cout << t[0] << " " << E[0] << endl;
  cout << t[Npoints-1] << " " << E[Npoints-1] << endl;

  //Separar los flares en dos arrays
  Int_t flare1_index=0;
  Int_t flare2_index=0;
  for(int i=0; i<numberOfEvents; i++){
	  if(flare[i]==1){
		  E1[flare1_index]=E[i];
		  t1[flare1_index]=t[i];
		  flare1_index++;
	  }
	  if(flare[i]==2){
		  E2[flare2_index]=E[i];
		  t2[flare2_index]=t[i];
		  flare2_index++;
	  }
  }

  //Representation of the flares
//  TH1D* energyEv1 = new TH1D("Data E dist1", "Data E dist1", 20, Emin1, Emax1);
//  TH1D* timeEv1 = new TH1D("Data t dist1", "Data t dist1", 12, tmin1, tmax1);
//  TH1D* energyEv2 = new TH1D("Data E dist2", "Data E dist2", 20, Emin2, Emax2);
//  TH1D* timeEv2 = new TH1D("Data t dist2", "Data t dist2", 12, tmin2, tmax2);

//  	for ( int a = 0; a < flare1_index; a++)
//  	{
//  		energyEv1->Fill(E1[a]);
//  		timeEv1->Fill(t1[a]);
//  	}
//
//  	for ( int a = 0; a < flare2_index; a++)
//  	  	{
//  	  		energyEv2->Fill(E2[a]);
//  	  		timeEv2->Fill(t2[a]);
//  	  	}

  	/*TCanvas* data1= new TCanvas();
  	data1->Divide(2,1);
  	data1->cd(1);
  	energyEv1->Draw("EP");
  	data1->cd(2);
  	timeEv1->Draw("EP");

  	TCanvas* data2= new TCanvas();
  	data2->Divide(2,1);
  	data2->cd(1);
	energyEv2->Draw("EP");
	data2->cd(2);
	timeEv2->Draw("EP");*/



  TMinuit *gMinuit = new TMinuit(5);
  gMinuit->SetFCN(fcn);

	double arglist[10];
	arglist[0]=0;
	int ierflg = 0;

	//Parameter definition

	Double_t QGdown=-QGpar*100;
	Double_t QGup=QGpar*100;
	gMinuit->mnparm(0, "QGpar", QGpar, 1., QGdown, QGup,ierflg);


	Double_t gausMaxdown1=gausMax1-3*gausWidth;
	Double_t gausMaxup1=gausMax1+3*gausWidth;
	gMinuit->mnparm(1, "gausMax1", gausMax1, 100., gausMaxdown1, gausMaxup1,ierflg);


	Double_t gausMaxdown2=gausMax2-3*gausWidth;
	Double_t gausMaxup2=gausMax2+3*gausWidth;
	gMinuit->mnparm(2, "gausMax2", gausMax2, 100., gausMaxdown2, gausMaxup2,ierflg);


	Double_t gausWidthdown=gausWidth/5.;
	Double_t gausWidthup=gausWidth*5.;
	gMinuit->mnparm(3, "gausWidth", gausWidth, 10., gausWidthdown, gausWidthup,ierflg);


	Double_t SBratiodown=0.01;
	Double_t SBratioup=1.;
	gMinuit->mnparm(4, "SBratio", SBratio, 0.5, SBratiodown, SBratioup,ierflg);



	//Fix parameters
//	  gMinuit->FixParameter(0);
//	  gMinuit->FixParameter(1);
//	  gMinuit->FixParameter(2);
//	  gMinuit->FixParameter(3);
//	  gMinuit->FixParameter(4);


	gMinuit->SetErrorDef(1);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	 //Collect value of parameters
	  Double_t QG0, QG, QGErr;

	  gMinuit->GetParameter(0, QG0, QGErr);

	  cout << "QG = " << QG0 << " Error = " << QGErr << endl;

	  Double_t fmin,edm,errdef,amin, bmin, cmin;
	  Int_t nvpar,nparx,icstat;
	  gMinuit->mnstat(fmin,edm,errdef,nvpar,nparx,icstat);
	  cout << "FCN min value = " << fmin << endl;

	  cout << "I fix the QG parameter and start the loop" << endl;

	  cout << "***************************************************************" << endl;
	  cout << "***************          LOOP RIGHT          ******************" << endl;
	  cout << "***************************************************************" << endl;

	  TGraph* GraphR = new TGraph();
	  Double_t DeltaLklR[100], QGR[100];
	  //Middle point of the graph
	  DeltaLklR[0]= 0;
	  QGR[0]= QG0;

	  GraphR->SetPoint(0, QGR[0], DeltaLklR[0]);
	  Int_t i1 = 1 ;
	  bmin = fmin;
	  QG = QG0;

	  while ( TMath::Abs(bmin-fmin) <= 2.)
	  {
		  QG+=2.;
		  gMinuit->mnparm(0, "QG" , QG, 0.1, QGdown, QGup,ierflg);
		  gMinuit->FixParameter(0);
		  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
		  gMinuit->mnstat(bmin,edm,errdef,nvpar,nparx,icstat);
		  gMinuit->GetParameter(0, QGR[i1], QGErr);
		  DeltaLklR[i1] = TMath::Abs(bmin-fmin);
		  GraphR->SetPoint(i1-1, QGR[i1], DeltaLklR[i1]);

		  cout << "***************************************************************" << endl;
		  cout << "***************           RESULTS            ******************" << endl;
		  cout << "***************************************************************" << endl;

		  cout << "Valor de fcn3 = " << bmin << endl;
		  cout << "DeltaL = " << DeltaLklR[i1] << endl;
		  cout << "QG = " << QGR[i1] << endl;
		  i1++;

	 }
	  i1-=1;


	cout << "***************************************************************" << endl;
	cout << "***************          LOOP LEFT           ******************" << endl;
	cout << "***************************************************************" << endl;

	TGraph* GraphL = new TGraph();
	Double_t DeltaLklL[100], QGL[100];
	Int_t i2 = 0 ;
	cmin = fmin;
	QG = QG0;
	while ( TMath::Abs(cmin-fmin) <= 2.)
	  {
		QG-=2.;
		gMinuit->Release(4);
		gMinuit->mnparm(0, "QG" , QG, 0.1, QGdown, QGup,ierflg);
		gMinuit->FixParameter(0);
		gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
		gMinuit->mnstat(cmin,edm,errdef,nvpar,nparx,icstat);
		gMinuit->GetParameter(0, QGL[i2], QGErr);
		DeltaLklL[i2] = TMath::Abs(cmin-fmin);
		GraphL->SetPoint(i2,QGL[i2], DeltaLklL[i2]);

		cout << "***************************************************************" << endl;
		cout << "***************           RESULTS            ******************" << endl;
		cout << "***************************************************************" << endl;

		cout << "Valor de fcn3 = " << cmin << endl;
		cout << "DeltaL = " << DeltaLklL[i2] << endl;
		cout << "QG = " << QGL[i2] << endl;
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
	Double_t DeltaLkl[300],QGcomp[300];
	TGraph* Graph = new TGraph();
	Int_t i3=0;
	while(i3<=i2)
	{
		 DeltaLkl[i3] = DeltaLklL[i2-i3];
		 QGcomp[i3] = QGL[i2-i3];
		 i3++;
	}

	Int_t i4=0;
	while(i4<=i1)
	{
		 DeltaLkl[i3+i4] = DeltaLklR[i4];
		 QGcomp[i3+i4] = QGR[i4];
		 i4++;
	}
	Double_t totalPoints=i3+i4-1;


	//Create a file with results
//	FILE *fout;
//	fout = fopen("results_alpha50.txt", "wb");
	for ( Int_t b=0; b<totalPoints; b++)
	 {
	   Graph->SetPoint(b, QGcomp[b], DeltaLkl[b]);
//	   fprintf(fout,"%0.8lf    %0.8lf    %0.8lf    %0.8lf    \n", QGcomp[b], Maxcomp[b], Witdhcomp[b], Ratiocomp[b]);
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

	TString rootoutfile("Graph_alpha0");
	TFile* outputfile = new TFile(rootoutfile, "recreate");
	canvas->Write();
	outputfile->Close();

//	cout << "Mp/Mqf = " << QG0 << " - " << QG0-err[0] << " + " << err[1]-QG0 << endl;







}

double PDFevent1(double t1, double E1, double *par)
{

  double QGParameter = par[0];
  double MaxParameter1 = par[1];
  double MaxParameter2 = par[2];
  double WidthParameter= par[3];
  double RatioParameter = par[4];


  Double_t QGDelay1= QGfixedPart*QGParameter*E1/PlanckScale;

  Double_t dt1= t1-QGDelay1;

  Double_t PDF1= RatioParameter*TMath::Power(E1,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(E1,signalIndex)*TMath::Exp(-(TMath::Power((dt1-MaxParameter1),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));


  //Normalization of the PDF

  PDF1=PDF1/(RatioParameter*(tmax1-tmin1)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);

  return PDF1;

}

double PDFevent2(double t2, double E2, double *par) //Hay que ver cómo hacemos para empezarlas a las vez
{

  double QGParameter = par[0];
  double MaxParameter1 = par[1];
  double MaxParameter2 = par[2];
  double WidthParameter= par[3];
  double RatioParameter = par[4];


  Double_t QGDelay2= QGfixedPart*QGParameter*E2/PlanckScale;

  Double_t dt2= t2-QGDelay2;

  Double_t PDF2= RatioParameter*TMath::Power(E2,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(E2,signalIndex)*TMath::Exp(-(TMath::Power((dt2-MaxParameter2),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));


  //Normalization of the PDF
  PDF2=PDF2/(RatioParameter*(tmax2-tmin2)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);

  return PDF2;

}

//double PDFevent_Res(double t, double E, double *par)
//{
//	//Limits of integration
//	Double_t Int_Emin = Emin/2.;
//	Double_t Int_Emax = Emax*2.;
//
////  double tauParameter = par[0];
//  double QGParameter = par[0];
//  double MaxParameter = par[1];
//  double WidthParameter= par[2];
//  double RatioParameter = par[3];
//
//
//	Double_t QGDelay= QGfixedPart*QGParameter/PlanckScale; //Now the energy is inside the function.
//
//	TF1 *f1 = new TF1("f1","([2]*TMath::Power(x,-2.7)+(1-[2])*50.*TMath::Power(x,-2.4)*(TMath::Exp(-(TMath::Power(([6]-[3]*x-[0]),2))/(2*TMath::Power([1],2)))*1/([1]*TMath::Sqrt(2.*TMath::Pi()))))*(TMath::Exp(-(TMath::Power((x-[4]),2))/(2*TMath::Power([5],2)))*(1/([5]*TMath::Sqrt(2.*TMath::Pi()))))",Emin,Emax);
//	f1->SetParameters(MaxParameter,WidthParameter,RatioParameter, QGDelay, E, res*E, t);
//	Double_t PDF= f1->Integral(Int_Emin, Int_Emax);
//	f1->Delete();
//
//  //Normalization of the PDF
//
//  PDF=PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);
//  return PDF;
//
//}
//
double PDFtest1(double t1, double E1)
{

	double QGParameter = QGpar;
    double MaxParameter1 = gausMax1;
    double MaxParameter2 = gausMax2;
    double WidthParameter= gausWidth;
    double RatioParameter = SBratio;


	  Double_t QGDelay1= QGfixedPart*QGParameter*E1/PlanckScale;

	  Double_t dt1= t1-QGDelay1;

	  Double_t PDF1= RatioParameter*TMath::Power(E1,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(E1,signalIndex)*TMath::Exp(-(TMath::Power((dt1-MaxParameter1),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));


	  //Normalization of the PDF

	  PDF1=PDF1/(RatioParameter*(tmax1-tmin1)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);

	  return PDF1;


}

double PDFtest2(double t2, double E2)
{

	double QGParameter = QGpar;
	double MaxParameter1 = gausMax1;
	double MaxParameter2 = gausMax2;
	double WidthParameter= gausWidth;
	double RatioParameter = SBratio;

	Double_t QGDelay2= QGfixedPart*QGParameter*E2/PlanckScale;

  Double_t dt2= t2-QGDelay2;

  Double_t PDF2= RatioParameter*TMath::Power(E2,bkgIndex)+(1-RatioParameter)*50.*TMath::Power(E2,signalIndex)*TMath::Exp(-(TMath::Power((dt2-MaxParameter2),2))/(2*TMath::Power(WidthParameter,2)))*1/(WidthParameter*TMath::Sqrt(2.*TMath::Pi()));


  //Normalization of the PDF
  PDF2=PDF2/(RatioParameter*(tmax2-tmin2)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);

  return PDF2;

}


//double PDFtestRes(double t, double E) //AHORA FUNCIONA BIEN
//{
//	//Limits of integration
//	Double_t Int_Emin = Emin/2.;
//	Double_t Int_Emax = Emax*2.;
//
//	double QGParameter = QGpar;
//	double MaxParameter = gausMax;
//	double WidthParameter= gausWidth;
//	double RatioParameter = SBratio;
//
//	Double_t QGDelay= QGfixedPart*QGParameter/PlanckScale; //Now the energy is inside the function.
//
//	TF1 *f1 = new TF1("f1","([2]*TMath::Power(x,-2.7)+(1-[2])*50.*TMath::Power(x,-2.4)*(TMath::Exp(-(TMath::Power(([6]-[3]*x-[0]),2))/(2*TMath::Power([1],2)))*1/([1]*TMath::Sqrt(2.*TMath::Pi()))))*(TMath::Exp(-(TMath::Power((x-[4]),2))/(2*TMath::Power([5],2)))*(1/([5]*TMath::Sqrt(2.*TMath::Pi()))))",Emin,Emax);
//	f1->SetParameters(MaxParameter,WidthParameter,RatioParameter, QGDelay, E, res*E, t);
//	Double_t PDF= f1->Integral(Int_Emin, Int_Emax);
//	f1->Delete();
//
//
////	Double_t PDF= f1->Eval(E);
//	//Normalization of the PDF
//	PDF=PDF/(RatioParameter*(tmax-tmin)*specIntegralbkg+(1-RatioParameter)*50.*specIntegralSignal);
//
//	return PDF;
//
//}
//
//
//double PDFIntegral(){
//
//	double MaxParameter2 = gausMax;//0
//	double WidthParameter2= gausWidth;//1
//	double RatioParameter2 = SBratio;//2
//
//	TF2 *f2 = new TF2("f","[2]*TMath::Power(y,-2.7)+(1-[2])*50.*TMath::Power(y,-2.4)*TMath::Exp(-(TMath::Power((x-[0]),2))/(2*TMath::Power([1],2)))*1/([1]*TMath::Sqrt(2.*TMath::Pi()))",tmin,tmax,Emin,Emax);
//	f2->SetParameters(MaxParameter2,WidthParameter2,RatioParameter2);
//	double norm = f2->Integral(tmin,tmax,Emin,Emax);
//	double Manelnorm =RatioParameter2*(tmax-tmin)*specIntegralbkg+(1-RatioParameter2)*50.*specIntegralSignal;
//	f2->Delete();
//
//	cout << "Integral norm: " << norm << endl;
//	cout << "Manel norm: "  << Manelnorm << endl;
//	return norm;
//
//}
//
Double_t Likelihoodtest(){
	double logL1, logL2 = 0.;
	double temp1, temp2 = 0.;
	for (int i=0; i<numberOfEvents_flare1; i++)
	{
		temp1 = TMath::Log(PDFtest1(t1[i],E1[i]));
		temp2 = TMath::Log(PDFtest2(t2[i],E2[i]));
		logL1+=(-2)*temp1;
		logL2+=(-2)*temp2;
	}


	  //    cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << logL << endl;

	return logL1+logL2;
}

//Double_t Likelihoodtest_Res(){
//	Double_t suma=0;
//	  for(int i=0; i<numberOfEvents; i++){
//		  suma+=(-2)*TMath::Log(PDFtestRes(t[i], E[i]));
//	  }
//	  return suma;
//}
//
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{

  double logL1=0;
  double logL2=0.;
  double temp1= 0.;
  double temp2 = 0.;

  for (int i=0; i<numberOfEvents_flare1; i++)
  {
    temp1 = TMath::Log(PDFevent1(t1[i],E1[i],par));
    temp2 = TMath::Log(PDFevent2(t2[i],E2[i],par));
    logL1+=(-2)*temp1;
    logL2+=(-2)*temp2;
  }


  f = logL1+logL2;  //Final function.

//  cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] << " " << par[4] <<  " " << f << endl;

}




