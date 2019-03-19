/*
 * twoDLike.C
 *
 *  Created on: Apr 5, 2016
 *      Author: lnogues
 */

#include "TMinuit.h"
using namespace std;

const Int_t Nevents = 2000;
Double_t E[Nevents], t[Nevents];
int Npoints = 0;
Double_t tmin = 0;
Double_t tmax = 2731;
Double_t Emin = 0.15;
Double_t Emax = 11.210;
TGraph* CollAr;

void twoDLike()
{
	char *fname = "./tau0/flaretrue0.txt";
  ifstream in;
  in.open(fname);
  double v1,v2;
  while(1)
    {
      in >> v1 >> v2;
      if (!in.good()) break;
      t[Npoints] = v1;
      E[Npoints] = v2;
      Npoints++;
    }

  cout << Npoints << endl;

  //cout << "t[0]: " << t[0] << " t[955]: " << t[955] << endl;
  //cout << "E[0]: " << E[0] << " E[955]: " << E[955] << endl;

  /*TString myPath = "/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Mrk501/";
  const TString sNewFile   = myPath+"MonoCollAr.root";
  TFile* fNewFile = new TFile(sNewFile);
  CollAr = (TGraph*)fNewFile->Get("Mono");*/

  TMinuit *gMinuit = new TMinuit(4);
  gMinuit->SetFCN(fcn);

  double arglist[10];
  int ierflg = 0;

    //  arglist[0] = 0.5;
    //  gMinuit->mnexcm("SET ERR",arglist, 1, ierflg);

  Double_t vstart[4]={2.,2010,800,3.};            // Parameter starting values
  Double_t vstep[4]={2.,2.,2.,2.};            // Parameter initial search step
  Double_t vllimit[4]={0.,0.,0.,10.};            // Par lower bound
  Double_t vulimit[4]={0.,0.,0.,-10.};      // Par upper bound

  /*norm=1.98;
	mean=2009.;
	sigma=686.9;
	LIV=0.;
   */

  gMinuit->mnparm(0, "Constant", vstart[0], vstep[0], vllimit[0], vulimit[0],ierflg);
  gMinuit->mnparm(1, "tmax"    , vstart[1], vstep[1], vllimit[1], vulimit[1],ierflg);
  gMinuit->mnparm(2, "sigma"   , vstart[2], vstep[2], vllimit[2], vulimit[2],ierflg);
  gMinuit->mnparm(3, "tau"   , vstart[3], vstep[3], vllimit[3], vulimit[3],ierflg);

  gMinuit->FixParameter(0);
  //gMinuit->FixParameter(2);
  gMinuit->FixParameter(1);


  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->SetPrintLevel(1);
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

}

/*Double_t signal2D(Double_t *x, Double_t *par)
{
	Double_t signal = par[0]*14.3*TMath::Power(10.,-6.)*TMath::Power((x[1]/0.3),-2.2)*exp(-0.5*TMath::Power(((x[0]-par[1]-par[3]*x[1])/par[2]),2));
	Double_t area = CollAr->Eval(x[1]);
	return signal*area;
}*/

double function(double x, double y, double *par) //Funcion que normaliza la PDF
{
  double Constant = par[0];
  double tma = par[1];
  double sigma = par[2];
  double tau = par[3];

  TF2* simpler = new TF2("simpler","TMath::Power(y,-2.7)+[0]*TMath::Power(y,-2.4)*exp(-0.5*TMath::Power(((x-[1]-[3]*y)/[2]),2))", tmin, tmax, Emin, Emax);
  simpler->SetParameters(Constant, tma, sigma, tau);
  double norm = simpler->Integral(tmin, tmax, Emin, Emax);

  /*TF2* complex = new TF2("complex", signal2D, tmin, tmax, Emin, Emax, 4);
  complex->SetParameters(Constant, tma, sigma, tau);
  double norm = complex->Integral(tmin, tmax, Emin, Emax);*/

  /*TF2 *f2 = new TF2("f","[0]+[1]*exp(-0.5*((x-[2]-[4]*y)/[3])**2)",0.,100.,0.1,10.);
  f2->SetParameters(Baseline,Constant,tmax,sigma,c1);
  double norm = f2->Integral(0.,100.,0.1,10.);
  f2->Delete();
  double value = Constant*exp(-0.5*((x-tmax-c1*y)/sigma)**2); //(Eval en x,y)*/

  //double value = complex->Eval(x,y);
  double value = simpler->Eval(x,y);
  value /= norm;

  //cout << "value " << value << " norm " << norm << endl;
  //complex->Delete();
  simpler->Delete();

  return value;
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{

  double logl = 0.;
  double temp;
  for (int i=0; i<Npoints; i++)
  {
    temp = TMath::Log(function(t[i],E[i],par)); //Log of the PDF.
    logl += temp; //Final probability function (addition of all events).
  }

  cout <<  par[0] << " " << par[1] << " " << par[2] << " " << par[3] <<  " " << -logl << endl;

  f = -logl;  //Final function.

}


