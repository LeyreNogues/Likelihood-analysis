/*
 * sim.C
 *
 *  Created on: Feb 19, 2016
 *      Author: lnogues
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
#include <TRandom3.h>
#include <TTree.h>
#include <MStatusArray.h>
#include <TStyle.h>
#include <TColor.h>

using namespace std;

//Common variables
const Int_t simnum = 300;
Double_t tmean[simnum], chigaus[simnum], ind[simnum], chipl[simnum], Nlc[simnum], Nsp[simnum];

Double_t gaus(Double_t* x, Double_t* par)
{
	Double_t funci = 1.*exp(-(pow(x[0]-par[0],2)/(2*pow(par[1],2))));
	return funci;
}

Double_t powerLaw(Double_t* x, Double_t* par)
{
  Double_t func = par[0]*TMath::Power(x[0]/par[1],-par[2]);
  return func;
}

Double_t powerLaw2(Double_t* x, Double_t* par)
{
  Double_t func2 = par[0]*TMath::Power(x[0],-par[1]);
  return func2;
}

Double_t line(Double_t* x, Double_t* par)
{
  Double_t constant = par[0];
  return constant;
}

void simulation(Int_t loop)
{

	 //ROOT OPTION
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1111);
	gStyle->SetOptTitle(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetLabelOffset(0.003,"y");
	gStyle->SetLabelOffset(0.01,"x");
	gStyle->SetTitleOffset(1.0,"x");
	gStyle->SetTitleOffset(1.0,"y");
	gStyle->SetTitleXSize(0.05);
	gStyle->SetTitleYSize(0.05);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetTickLength(0.01,"y");
	gStyle->SetStatX(0.99);
	gStyle->SetStatY(0.99);
	gStyle->SetStatH(0.1);
	gStyle->SetStatW(0.1);


	gStyle->SetStatH(0.2);
	gStyle->SetStatW(0.2);
	gStyle->SetStatFontSize(0.03);


//--------------------------------------Lightcurve------------------------------------------


Double_t tmax, sigma, tinit;
tmax = 2005;
sigma = 219;
tinit = 1200;

//TCanvas* lightcurvecanvas = new TCanvas();
TF1 *lightcurve = new TF1("lightcurve", gaus, 1200-tinit, 2731.-tinit, 2);
lightcurve->SetParameters(tmax-tinit, sigma);
lightcurve->SetTitle("Lightcurve");
lightcurve->GetXaxis()->SetTitle("Time(s)");
//lightcurve->Draw();

/*TLatex *txt = new TLatex();

txt->DrawLatex(1250.-tinit, 0.95, "#mu = 805");
txt->DrawLatex(1250.-tinit, 0.85, "#sigma = 219");
txt->SetNDC(kTRUE);
txt->SetTextAlign(13);
txt->SetTextSize(15);

txt->Draw("SAME");*/

//-----------------------------------Generation of time events-------------------------------
gRandom->SetSeed(0);
double t[1490];
for(int i = 0; i<1490; i++)
    {
    t[i]=lightcurve->GetRandom();
    }

TH1D* timehisto = new TH1D("timehisto", "timehisto", 50, 1200.-tinit, 2731.-tinit);
for(int i = 0; i<1490; i++)
{
	timehisto->Fill(t[i]);
}

/*TCanvas* timehistocanvas = new TCanvas();
timehisto->SetTitle("Time generated events");
timehisto->GetXaxis()->SetTitle("Time(s)");
timehisto->GetYaxis()->SetTitle("#events");
timehisto->Fit("gaus");
timehisto->Draw("EP");

TLatex *txt2 = new TLatex();

txt2->DrawText(1250.-tinit, 52., "730+760 = 1490 Events");
txt2->SetNDC(kTRUE);
txt2->SetTextAlign(13);
txt2->SetTextSize(7);

txt2->Draw("SAME");*/

cout << "Time events generated" << endl;

//----------------------------------------------Spectrum------------------------------------

Double_t Emin,Emax;
Emin = 0.15; //First event at 83
Emax = 11.210; //Last event at 11208

Double_t norm, index, E0;
norm = 14.3e-6;
index = 2.2;
E0 = 0.3; //GeV

//TCanvas* spectrumcanvas = new TCanvas();
//spectrumcanvas->SetLogy();
//spectrumcanvas->SetLogx();
//spectrum->Draw();

TF1* spectrum = new TF1("spectrum", powerLaw, Emin, Emax, 3);
spectrum->SetParameters(norm, E0,index);
spectrum->SetTitle("Spectrum");
spectrum->GetXaxis()->SetTitle("Energy(TeV)");
spectrum->GetYaxis()->SetTitle("s^{-1}m^{-2} 0.3TeV^{-1}");

//Used for fits ( for a desired energy range)
TF1* spectrum2 = new TF1("spectrum2", powerLaw2, 0.3, 8., 2);
spectrum2->SetName("spectrum2");
spectrum2->SetParName(0, "norm");
spectrum2->SetParName(1, "index");


/*TLatex *txt3 = new TLatex();

txt3->DrawLatex(1.5, 4.e-5, "f_{0} = 14.3e-6 m^{-2}s^{-1}TeV^{-1}");
txt3->DrawLatex(1.5, 1.7e-5, "E_{0} = 0.3 TeV");
txt3->DrawLatex(1.5, 8.e-6, "#Gamma = 2.2");
txt3->SetNDC(kTRUE);
txt3->SetTextAlign(13);
txt3->SetTextSize(7);

txt3->Draw("SAME");*/


//-----------------------------------Generation of energy events-------------------------------

double XX_step = 0.1;
int n = ceil((log10(Emax)-log10(Emin))/XX_step);
double *XX = new double[n+1];


for(int i=0;i<n+1;i++){
	XX[i]= pow(10,(log10(Emin)+i*XX_step));
}


double E[1490];
for(int i = 0; i<1490; i++)
    {
    E[i]=spectrum->GetRandom(Emin, Emax);
    }

TH1D* energyhisto = new TH1D("energyhisto", "energyhisto", n, XX);
energyhisto->Sumw2();

for(int i = 0; i<1490; i++)
{
	energyhisto->Fill(E[i], 1./E[i]);
}

/*TCanvas* energyhistocanvas = new TCanvas();
energyhistocanvas->SetLogy();
energyhistocanvas->SetLogx();
energyhisto->SetTitle("Energy generated events");
energyhisto->GetXaxis()->SetTitle("Energy(TeV)");
energyhisto->GetYaxis()->SetTitle("#events");
energyhisto->Fit(spectrum, "Q");
energyhisto->SetMarkerColor(1);
energyhisto->SetLineColor(1);
energyhisto->Draw("EP");*/


/*TLatex *txt4 = new TLatex();

txt4->DrawLatex(0.3, 4., "1490 Events (0.15-10 TeV)");
txt4->SetNDC(kTRUE);
txt4->SetTextAlign(13);
txt4->SetTextSize(0.5);

txt4->Draw("SAME");*/

cout << "Energy events generated" << endl;


//----------------------------------- Stereo Collection Area factor-------------------------------

/*TString myPath = "/home/lnogues/workspace_cpp/minuit";
const TString sNewFile   = myPath+"/IRF/CollArMAGIC.root";
TFile* fNewFile = new TFile(sNewFile);
fNewFile->Get("c1");
TH1D *th1 = (TH1D*)c1->GetPrimitive("NewCollAr");*/
/*TCanvas* area = new TCanvas();
area->SetLogx();
area->SetLogy();
th1->Draw();*/


//----------------------------------- Mono Collection Area factor-------------------------------

TGraphErrors *gre = new TGraphErrors(28);
  gre->SetName("Graph_from_mgg1");
  gre->SetTitle("MonoEffectiveCollectionArea");

  Int_t ci;   // for color index setting
  ci = TColor::GetColor("#cc0000");
  gre->SetMarkerColor(ci);
  gre->SetMarkerStyle(8);
  gre->SetMarkerSize(0.7);
  gre->SetPoint(0,5.902727,59.79039);
  gre->SetPointError(0,0.9737387,0);
  gre->SetPoint(1,8.262956,59.79039);
  gre->SetPointError(1,1.353006,0);
  gre->SetPoint(2,11.35279,0.05456552);
  gre->SetPointError(2,1.879995,0.05451783);
  gre->SetPoint(3,15.89225,0.2445037);
  gre->SetPointError(3,2.612245,0.171296);
  gre->SetPoint(4,22.24683,5.255659);
  gre->SetPointError(4,3.629702,1.073794);
  gre->SetPoint(5,30.56576,48.73137);
  gre->SetPointError(5,5.043455,4.17656);
  gre->SetPoint(6,42.7876,344.001);
  gre->SetPointError(6,7.007858,14.2968);
  gre->SetPoint(7,58.7875,1613.114);
  gre->SetPointError(7,9.737387,39.08882);
  gre->SetPoint(8,82.29393,5258.519);
  gre->SetPointError(8,13.53006,89.6552);
  gre->SetPoint(9,115.1995,13050.63);
  gre->SetPointError(9,18.79995,178.9534);
  gre->SetPoint(10,158.2769,25225.41);
  gre->SetPointError(10,26.12245,319.5061);
  gre->SetPoint(11,221.5646,38846.59);
  gre->SetPointError(11,36.29702,494.293);
  gre->SetPoint(12,304.416,53397.6);
  gre->SetPointError(12,50.43455,736.6419);
  gre->SetPoint(13,426.138,67021.43);
  gre->SetPointError(13,70.07858,1041.433);
  gre->SetPoint(14,596.5311,78577.48);
  gre->SetPointError(14,97.37387,1432.469);
  gre->SetPoint(15,819.5966,84121.24);
  gre->SetPointError(15,135.3006,1885.549);
  gre->SetPoint(16,1147.316,86054.76);
  gre->SetPointError(16,187.9995,2476.141);
  gre->SetPoint(17,1576.34,86054.76);
  gre->SetPointError(17,261.2245,3170.326);
  gre->SetPoint(18,2206.647,86054.76);
  gre->SetPointError(18,362.9702,4092.978);
  gre->SetPoint(19,3088.983,96409.73);
  gre->SetPointError(19,504.3455,5401.472);
  gre->SetPoint(20,4244.071,96409.73);
  gre->SetPointError(20,700.7858,7070.007);
  gre->SetPoint(21,5941.081,59.79039);
  gre->SetPointError(21,973.7387,0);
  gre->SetPoint(22,8162.675,59.79039);
  gre->SetPointError(22,1353.006,0);
  gre->SetPoint(23,11426.55,59.79039);
  gre->SetPointError(23,1879.995,0);
  gre->SetPoint(24,15995.51,59.79039);
  gre->SetPointError(24,2612.245,0);
  gre->SetPoint(25,21976.84,59.79039);
  gre->SetPointError(25,3629.702,0);
  gre->SetPoint(26,30764.36,59.79039);
  gre->SetPointError(26,5043.455,0);
  gre->SetPoint(27,42268.32,59.79039);
  gre->SetPointError(27,7007.858,0);

  /*TCanvas *CollArOldcanvas = new TCanvas();
  CollArOldcanvas->SetLogy();
  CollArOldcanvas->SetLogx();*/

  TGraph* CollArOld = new TGraph();
  Double_t x,y;
  for(Int_t a = 0; a<28; a++) //Represent in LogE
  {
	  gre->GetPoint(a, x, y);
	  if(y!= 59.79039) CollArOld->SetPoint(a, x/1000., y);
	  else CollArOld->SetPoint(a, x/1000., 0);
  }
  /*CollArOld->SetTitle("MonoCollectionArea");
  CollArOld->SetName("Mono");
  CollArOld->GetXaxis()->SetTitle("E_{true}(TeV)");
  CollArOld->GetYaxis()->SetTitle("Aeff(m^{2})");
  CollArOld->SetLineColor(2);
  CollArOld->SetLineWidth(2);
  CollArOld->Draw("");*/

  /*TString rootoutfile("MonoCollAr.root");
  TFile* outputfile = new TFile(rootoutfile, "recreate");
  CollArOld->Write();
  outputfile->Close();*/

//Get maximum value
  Double_t fmax = 96409.73;

//cout << "The maximum of the Collection Area is: " << fmax << " m²" << endl;

//-----------------------------------Acceptance factors of the events----------------------------------

Double_t factor[1490];
Double_t total = 0;
for(int i = 0; i<1490; i++)
{
	factor[i] = CollArOld->Eval(E[i])/fmax;
	total += factor[i];
}

Double_t meanfactor = total/1490;
cout << "Mean factor " << meanfactor << endl;

TH1D* factorhisto = new TH1D("factorhisto", "factorhisto", 100, 0., 1.);
for(int i = 0; i<1490; i++)
{
	factorhisto->Fill(factor[i]);
}

TH1D* meanfactorhisto = new TH1D("meanfactorhisto", "meanfactorhisto", 100, 0., 2.);
for(int i = 0; i<1490; i++)
{
	meanfactorhisto->Fill(factor[i]/meanfactor);
}

TH1D* weighthisto = new TH1D("weighthisto", "weighthisto", n, XX);
TH1D* weighthisto2 = new TH1D("weighthisto2", "Energy template",  n, XX);

for(int i = 0; i<1490; i++)
{
	weighthisto->Fill(E[i], factor[i]/E[i]);
	weighthisto2->Fill(E[i], (factor[i]/meanfactor)/E[i]);
}

/*TCanvas* weightfitcanvas = new TCanvas();
weightfitcanvas->SetLogx();
weightfitcanvas->SetLogy();
weighthisto->SetLineColor(2);
weighthisto2->SetLineColor(4);
//weighthisto->Draw("EPsame");
weighthisto2->Fit(spectrum2, "R");
weighthisto2->Draw("EP");*/

/*TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
leg1->AddEntry(energyhisto,"True Events","le");
leg1->AddEntry(weighthisto,"Weighted by f_{A}","le");
leg1->AddEntry(weighthisto2,"Weighted by f'_{A}","le");
leg1->Draw();*/



/*TCanvas* factorcanvas = new TCanvas();
factorcanvas->Divide(1,2);
factorcanvas->cd(1);
factorhisto->SetTitle("Acceptance factor");
factorhisto->GetXaxis()->SetTitle("f_{A}");
factorhisto->GetYaxis()->SetTitle("#events");
factorhisto->Draw("EP");

TLatex *txt4 = new TLatex();
txt4->DrawLatex(0.6, 80., "<f_{A} = 0.52 >");
txt4->SetNDC(kTRUE);
txt4->SetTextAlign(13);
txt4->Draw("SAME");


factorcanvas->cd(2);
meanfactorhisto->SetTitle("Corrected Acceptance factor");
meanfactorhisto->GetXaxis()->SetTitle("f'_{A}");
meanfactorhisto->GetYaxis()->SetTitle("#events");
meanfactorhisto->Draw("EP");*/


cout << "Acceptance factor computed" << endl;



//--------------------------------------Stereo Migration Matrix-------------------------------------


/*TString myPath = "/home/lnogues/workspace_cpp/minuit";
const TString sNewFile   = myPath+"/IRF/MigMatrix.root";
TFile* fNewFile = new TFile(sNewFile);
fNewFile->Get("c1");
TH2D *th2 = (TH2D*)c1->GetPrimitive("migrmatrix_yx__1");
/*TCanvas* area = new TCanvas();
area->SetLogx();
area->SetLogy();
th2->Draw("zcol");*/



//--------------------------------------Mono Migration Matrix-------------------------------------

TFile* fileStatusFluxlc = TFile::Open("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Status_Output_fluxlc.root","READ");
MStatusArray* StatusOutputFluxlc = fileStatusFluxlc->Get("MStatusDisplay");
TCanvas * MigMatrix = (TCanvas *) StatusOutputFluxlc->FindCanvas("Migration Matrix");
TH2D* mm = (TH2D*)MigMatrix->GetPrimitive("mimgplot");
/*TCanvas * c2 = new TCanvas();
c2->SetLogx();
c2->SetLogy();
mm->Draw("colz");*/


//Loop the axis of Etrue
Double_t oldmean[17], oldSigma[17], oldEtrue[17], oldbias[17], oldsigmaE[17];
Int_t oldbin;
Int_t initbin = 5;
TH1D* oldprojX;
TF1* oldmyfunc;

//TCanvas* oldprojXcanvas = new TCanvas();
for(Int_t i = initbin; i<=18; i++)
{
	oldbin = i;
	oldprojX = mm->ProjectionX("oldprojX", oldbin, oldbin, "");
	oldprojX->Fit("gaus", "0Q");

	oldmyfunc = oldprojX->GetFunction("gaus");
	Double_t oldpar1 = oldmyfunc->GetParameter(0);
	Double_t oldpar2 = oldmyfunc->GetParameter(1);
	Double_t oldpar3 = oldmyfunc->GetParameter(2);

	oldEtrue[oldbin-initbin] = mm->GetYaxis()->GetBinCenter(oldbin);
	oldmean[oldbin-initbin] = oldpar2;
	oldSigma[oldbin-initbin] = oldpar3;
	oldbias[oldbin-initbin] = oldmean[oldbin-initbin]/oldEtrue[oldbin-initbin];
	oldsigmaE[oldbin-initbin] = oldSigma[oldbin-initbin]/oldEtrue[oldbin-initbin];

	//cout << " oldEtrue: " << oldEtrue[oldbin-initbin] << " oldMean: " << oldmean[oldbin-initbin] << " oldSigma: " << oldSigma[oldbin-initbin] << endl;
	//cout << "oldErec/oldEtrue: " << oldbias[oldbin-initbin] << endl;
	//cout << "oldsigmaE: " << oldsigmaE[oldbin-initbin] << endl;
}

Int_t count = 0;
Double_t oldadd = 0;
for(Int_t point = 0; point <=oldbin-initbin; point++)
{
	if(oldEtrue[point]>150 &&oldEtrue[point]<11000)
		{
			oldadd+= oldsigmaE[point];
			count++;
		}
}

Double_t res = oldadd/count;
//cout << "Mean resolution: " << res << endl;

cout << "Bias and resolution computed" << endl;


//-------------------------------------Mono Energy bias and resolution-------------------------------

TGraph* shift = new TGraph();
TGraph* resol = new TGraph();

for(Int_t point = 0; point <=oldbin-initbin; point++)
{
	shift->SetPoint(point, oldEtrue[point]/1000., oldbias[point]);
	resol->SetPoint(point, oldEtrue[point]/1000., oldsigmaE[point]);
}

Double_t biasLC, biasFit = 0;
Int_t iLC, iFit = 0;
for(Int_t point = 0; point <=oldbin-initbin; point++)
{
	if(oldEtrue[point]< 150)
	{
		iLC++;
		biasLC += oldbias[point];
	}
	if(oldEtrue[point]>=150 && oldEtrue[point]<= Emax)
	{
		iFit++;
		biasFit += oldbias[point];
	}
}

Double_t finalbias= biasLC/iLC;

//cout << "Mean bias LC: " << biasLC/iLC << " Mean bias Fit: " << biasFit/iFit << endl;


/*TCanvas* shiftcanvas = new TCanvas();
shiftcanvas->SetLogx();
shift->SetTitle("Mono Energy bias");
shift->SetLineColor(2);
shift->SetLineWidth(2);
shift->GetXaxis()->SetTitle("E_{true}(TeV)");
shift->GetYaxis()->SetTitle("E_{rec}/E_{true}");
shift->Draw("");

TCanvas* resolcanvas = new TCanvas();
resolcanvas->SetLogx();
resol->SetTitle("Mono Energy resolution");
resol->SetLineColor(2);
resol->SetLineWidth(2);
resol->GetXaxis()->SetTitle("E_{true}(TeV)");
resol->GetYaxis()->SetTitle("#sigma_{E}/E_{true}");
resol->Draw("");*/

//-------------------------------------Smear of the events-------------------------------


cout << "Starting smearing of events..." << endl;
cout << "Resolution " << res << endl;

Double_t Esmeared[1490];
Double_t limit = 1490;

double XX2_step = 0.1;
int n2 = ceil((log10(Emax)-(log10(Emin)-0.6))/XX2_step);
double *XX2 = new double[n2+1];
for(int i=0;i<n2+1;i++)
{
	XX2[i]= pow(10,((log10(Emin)-0.6)+i*XX2_step));
}
//TCanvas* gauscanvas = new TCanvas();
//gauscanvas->DrawFrame(0.02,0,0.8,1.2);
TF1 *gausian = new TF1("gausian", "gaus", Emin-2., Emax+2.);
//gausian->SetParameters(1, 0.2, res*0.2);
//gausian->Draw("SAME");
//gauscanvas->SaveAs("test.pdf");

for(Int_t a=0; a<limit; a++)
{

	gausian->SetParameters(1, E[a], res*E[a]);
	Esmeared[a] = gausian->GetRandom();
	//cout << "True energy (TeV) " << E[a] << endl;
	//cout << "Reconstructed energy(TeV) " << Esmeared[a] << endl;

}

TH1D* orihisto = new TH1D("orihisto", "orihisto", n2, XX2);
TH1D* recgaushisto = new TH1D("recgaushisto", "recgaushisto", n2, XX2);

for(int i = 0; i<limit; i++)
{
	orihisto->Fill(E[i], 1/E[i]);
	recgaushisto->Fill(Esmeared[i], 1/Esmeared[i]);
}
orihisto->Sumw2();
recgaushisto->Sumw2();

/*TCanvas* reccanvas = new TCanvas();
reccanvas->SetLogy();
reccanvas->SetLogx();
orihisto->SetTitle("True vs Reconstructed events");
orihisto->GetXaxis()->SetTitle("Energy(TeV)");
orihisto->GetYaxis()->SetTitle("#Events");
orihisto->SetLineColor(1);
orihisto->Draw("EP");
recgaushisto->SetLineColor(2);
recgaushisto->Draw("sameEP");

TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
leg1->AddEntry(orihisto,"True Events","le");
leg1->AddEntry(recgaushisto,"Rec Events","le");
leg1->Draw();*/

//cout << "Integralori: " << orihisto->Integral()/orihisto->GetSumOfWeights() << endl;
//cout << "Integralrec: " << recgaushisto->Integral()/orihisto->GetSumOfWeights() << endl;
cout << "Smearing finished" << endl;

//---------------------------------------Final Histogram-----------------------------------
Double_t factor2[1490];
TH1D* final = new TH1D("final", "final", n2, XX2);
for(int i = 0; i<limit; i++)
{
	factor2[i]=factor[i]/meanfactor;
	final->Fill(Esmeared[i], factor2[i]/Esmeared[i]);
}

//TCanvas* finalcanvas = new TCanvas();
//finalcanvas->SetLogx();
//finalcanvas->SetLogy();
final->GetXaxis()->SetTitle("E_{true}(TeV)");
final->GetYaxis()->SetTitle("dN/dE");
final->SetName("PLtemplate");
final->SetMarkerStyle(20);
final->SetMarkerColor(2);
final->SetLineColor(2);
final->Fit(spectrum2, "QR");
//final->Draw("EP");
/*recgaushisto->SetMarkerStyle(20);
recgaushisto->SetMarkerColor(3);
recgaushisto->SetLineColor(3);
recgaushisto->Fit(spectrum2, "R");
recgaushisto->Draw("EP");*/
/*orihisto->SetLineColor(1);
orihisto->SetMarkerStyle(20);
orihisto->SetMarkerColor(1);
orihisto->Fit(spectrum2, "R");
orihisto->Draw("EP");*/
/*weighthisto2->SetLineColor(4);
weighthisto2->SetMarkerStyle(20);
weighthisto2->SetMarkerColor(4);
weighthisto2->Fit(spectrum2, "R");
weighthisto2->Draw("EP");*/

//Obtain parameters of the fit
TF1* fitpl = final->GetFunction("spectrum2");
ind[loop] = fitpl->GetParameter(1);
chipl[loop] = (fitpl->GetChisquare())/(fitpl->GetNDF());

TString rootoutfile(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/tau0/Tree/PLtemplater%i.root", loop));
TFile* outputfile = new TFile(rootoutfile, "recreate");
final->Write();
outputfile->Close();

/*TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
leg1->AddEntry(orihisto,"True Events","lp");
leg1->AddEntry(weighthisto2,"W Events","lp");
leg1->AddEntry(recgaushisto,"S Events","lp");
leg1->AddEntry(final,"W+S Events","lp");
leg1->Draw();*/

cout << "Final Graph drawn" << endl;

//-------------------------------------Flare with EWeight----------------------------------

//Generation of events
/*Double_t Ew[1490];
TH1D* weightsim = new TH1D("weightsim", "weightsim", n, XX);

for(int i = 0; i<limit; i++)
{
	Ew[i] = weighthisto2->GetRandom();
	weightsim->Fill(Ew[i]);
}*/

/*TCanvas* weightflarecanvas = new TCanvas();
weightflarecanvas->SetLogx();
weightflarecanvas->SetLogy();
weighthisto2->SetTitle("Weighted simulated events");
weighthisto2->GetXaxis()->SetTitle("E_{true}(TeV)");
weighthisto2->GetYaxis()->SetTitle("dN/dE");
weighthisto2->SetLineColor(1);
weightsim->SetLineColor(4);
weighthisto2->Draw("");
weightsim->Draw("sameEP");

TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
leg1->AddEntry(weighthisto2,"Histo","l");
leg1->AddEntry(weightsim,"Sim Events","le");
leg1->Draw();*/

//Division of events
/*Int_t LCphotonsw = 0;
Int_t Fitphotonsw = 0;
TH1D* LCtemplatew = new TH1D("LCtemplatew", "LCtemplatew", 50, 1200.-tinit, 2731.-tinit);
TH1D* energysamplew = new TH1D("energysamplew", "energysamplew", n, XX);
TH1D* timesamplew = new TH1D("timesamplew", "timesamplew", 50, 1200.-tinit, 2731.-tinit);

Double_t Esimw[800];
Double_t tsimw[800];
Int_t i2 = 0;*/

/*for(int i = 0; i<limit; i++)
{
	if(Ew[i] >= 0.15 && Ew[i] <= 0.25)
	{
		LCtemplatew->Fill(t[i]);
		LCphotonsw++;
	}
	if(Ew[i] > 0.25 && Ew[i] <= 4.)
	{
		Fitphotonsw++;
		Esimw[i2] = Ew[i];
		tsimw[i2] = t[i];
		timesamplew->Fill(tsimw[i2]);
		energysamplew->Fill(Esimw[i2], factor2[i]/Esimw[i2]);
		i2++;
	}
}*/
/*cout << "Number of photons in LC (Weight Case): " << LCphotonsw << endl;
cout << "Number of photons in Fit (Weight Case): " << Fitphotonsw << endl;
cout << "i2: " << i2 << endl;*/

/*LCtemplatew->Sumw2();
LCtemplatew->Scale(1./LCtemplatew->Integral());*/


//Write a txt file with the evetns
/*FILE *fout;
fout = fopen("weightflare.txt","wb");
fprintf(fout,"%0.3lf \n", LCphotonsw);
fprintf(fout,"%0.3lf \n", Fitphotonsw);
for(int i = 0; i<i2; i++)
{
	fprintf(fout,"%0.3lf    %0.8lf   \n", tsimw[i], Esimw[i]);
}*/


/*TCanvas* LCtemplatecanvasw = new TCanvas();
LCtemplatew->SetTitle("LC template");
LCtemplatew->SetName("LCtemplate");
LCtemplatew->GetXaxis()->SetTitle("Time (s)");
LCtemplatew->GetYaxis()->SetTitle("dN/dt");
LCtemplatew->Fit("gaus", "Q");
LCtemplatew->Draw("EP");*/

/*TString rootoutfile("LCtemplatew.root");
TFile* outputfile = new TFile(rootoutfile, "recreate");
LCtemplatew->Write();
outputfile->Close();*/

/*TCanvas* Fitcanvasw = new TCanvas();
Fitcanvasw->Divide(1,2);
Fitcanvasw->cd(1);
energysamplew->SetTitle("Energy events for the fit");
energysamplew->GetXaxis()->SetTitle("E_{true}(TeV)");
energysamplew->GetYaxis()->SetTitle("dN/dE");
Fitcanvasw->SetLogx();
Fitcanvasw->SetLogy();
energysamplew->Draw("EP");
Fitcanvasw->cd(2);
timesamplew->SetTitle("Time events for the fit");
timesamplew->GetXaxis()->SetTitle("Time (s)");
timesamplew->GetYaxis()->SetTitle("dN/dt");
timesamplew->Draw("EP");*/




//-------------------------------------Flare with Erec-----------------------------------

//Generation of events
Double_t Er[1490];
TH1D* recsim = new TH1D("recsim", "recsim", n2, XX2);

for(int i = 0; i<limit; i++)
{
	Er[i] = final->GetRandom();
	recsim->Fill(Er[i]);
}

/*TCanvas* recflarecanvas = new TCanvas();
recflarecanvas->SetLogx();
recflarecanvas->SetLogy();
final->SetTitle("Rec simulated events");
final->GetXaxis()->SetTitle("E_{true}(TeV)");
final->GetYaxis()->SetTitle("dN/dE");
final->SetLineColor(1);
recsim->SetLineColor(4);
final->Draw("");
recsim->Draw("sameEP");

TLegend* leg1 = new TLegend(0.1,0.7,0.4,0.9);
leg1->AddEntry(weighthisto2,"Histo","l");
leg1->AddEntry(weightsim,"Sim Events","le");
leg1->Draw();*/

//Division of events
Int_t LCphotonsr = 0;
Int_t Fitphotonsr = 0;
TH1D* LCtemplater = new TH1D("LCtemplater", "LCtemplater", 50, 1200.-tinit, 2731.-tinit);
TH1D* energysampler = new TH1D("energysampler", "energysampler", n, XX);
TH1D* timesampler = new TH1D("timesampler", "timesampler", 50, 1200.-tinit, 2731.-tinit);

Double_t Esimr[760];
Double_t tsimr[760];
Int_t i3 = 0;

for(int i = 0; i<limit; i++)
{
	if(Er[i] >= 0.15 && Er[i] <= 0.3)
	{
		LCtemplater->Fill(t[i]);
		LCphotonsr++;
	}
	if(Er[i] > 0.3 && Er[i] <= 8.)
	{
		Esimr[i3] = Er[i];
		tsimr[i3] = t[i];
		Fitphotonsr++;
		timesampler->Fill(tsimr[i3]);
		energysampler->Fill(Esimr[i3], factor2[i]/Esimr[i3]);
		i3++;
	}
}
cout << "Number of photons in LC (Rec Case): " << LCphotonsr << endl;
cout << "Number of photons in Fit (Rec Case): " << Fitphotonsr << endl;
cout << "i3: " << i3 << endl;

LCtemplater->Sumw2();
LCtemplater->Scale(1./LCtemplater->Integral());

Nlc[loop] = LCphotonsr;
Nsp[loop] = Fitphotonsr;

//Write a txt file with the evetns
FILE *fout;
fout = fopen(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/tau0/Tree/recflare%i.txt", loop),"wb");
fprintf(fout,"%0.3lf \n", LCphotonsr);
fprintf(fout,"%0.3lf \n", Fitphotonsr);
for(int i = 0; i<i3; i++)
{
	fprintf(fout,"%0.3lf    %0.8lf   \n", tsimr[i], Esimr[i]);
}

//TCanvas* LCtemplatecanvasr = new TCanvas();
LCtemplater->SetTitle("LC template");
LCtemplater->SetName("LCtemplate");
LCtemplater->GetXaxis()->SetTitle("Time (s)");
LCtemplater->GetYaxis()->SetTitle("dN/dt");
LCtemplater->Fit("gaus", "Q");
//LCtemplater->Draw("EP");

//Obtain parameters of the fit
TF1* fitlc = LCtemplater->GetFunction("gaus");
tmean[loop] = fitlc->GetParameter(1);
chigaus[loop] = (fitlc->GetChisquare())/(fitlc->GetNDF());

TString rootoutfile(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/tau0/Tree/LCtemplater%i.root", loop));
TFile* outputfile = new TFile(rootoutfile, "recreate");
LCtemplater->Write();
outputfile->Close();

/*TCanvas* Fitcanvasr = new TCanvas();
Fitcanvasr->Divide(1,2);
Fitcanvasr->cd(1);
energysampler->SetTitle("Energy events for the fit");
energysampler->GetXaxis()->SetTitle("E_{true}(TeV)");
energysampler->GetYaxis()->SetTitle("dN/dE");
Fitcanvasr->SetLogx();
Fitcanvasr->SetLogy();
energysampler->Draw("EP");
Fitcanvasr->cd(2);
timesampler->SetTitle("Time events for the fit");
timesampler->GetXaxis()->SetTitle("Time (s)");
timesampler->GetYaxis()->SetTitle("dN/dt");
timesampler->Draw("EP");*/

cout << "Resimulation and saving of the Likelihood inputs done" << endl;

delete timehisto;
delete energyhisto;
delete factorhisto;
delete meanfactorhisto;
delete weighthisto;
delete weighthisto2;
delete recsim;
delete energysampler;
delete timesampler;


}

void sim()
{
	Int_t iter;
	///home/lnogues/workspace_cpp/LIVanalysis/Simulation/tau0/Tree/
	TString rootoutfile("Simtau0.root");
	TFile* outputfile = new TFile(rootoutfile, "recreate");
	TTree* tree = new TTree("tree", "tree");
	Double_t t, chig, in, chip, nl, nf;
	tree->Branch("tmean", &t, "t/D");;
	TBranch *branchcg = tree->Branch("chigaus", &chig, "chig/D");
	TBranch *branchind = tree->Branch("ind", &in, "in/D");
	TBranch *branchcp = tree->Branch("chipl", &chip, "chip/D");
	TBranch *branchnlc = tree->Branch("Nlc", &nl, "nl/D");
	TBranch *branchnsp = tree->Branch("Nsp", &nf, "nf/D");


	for(Int_t i = 0; i< 10; i++)
	{
		iter = i;
		cout << "----------------------------- SIMULATION " << iter+1 << " --------------------------" << endl;
		simulation(iter);
		cout << "tmean " <<tmean[iter] << endl;
		cout << "chigaus " << chigaus[iter] << endl;
		cout << "ind " << ind[iter] << endl;
		cout << "chipl " << chipl[iter] << endl;
		cout << "Nlc " << Nlc[iter] << endl;
		cout << "Nsp " << Nsp[iter] << endl;
	}

	Double_t meantmean, meanchigaus, meanind, meanchipl, meanNlc, meanNsp;
	Int_t sumtmean, sumchigaus, sumind, sumchipl, sumNlc, sumNsp = 0;

	for(Int_t i = 0; i<=iter; i++)
	{
		sumtmean+=tmean[i];
		sumchigaus+=chigaus[i];
		sumind+=ind[i];
		sumchipl+=chipl[i];
		sumNlc+=Nlc[i];
		sumNsp+=Nsp[i];
	}

	meantmean=sumtmean/iter;
	meanchigaus=sumchigaus/iter;
	meanind=sumind/iter;
	meanchipl=sumchipl/iter;
	meanNlc=sumNlc/iter;
	meanNsp=sumNsp/iter;

	TH1D* tmeanhisto = new TH1D("tmeanhisto", "tmeanhisto", 50, 760, 860);
	TH1D* chigaushisto = new TH1D("chigaushisto", "chigaushisto", 50, 0, meanchigaus+2*meanchigaus);
	TH1D* indhisto = new TH1D("indhisto", "indhisto", 50, 1., 3.);
	TH1D* chiplhisto = new TH1D("chiplhisto", "chiplhisto", 50, 0, meanchipl+2*meanchipl);
	TH1D* Nlchisto = new TH1D("Nlchisto", "Nlchisto", 50, 500, 900);
	TH1D* Nsphisto = new TH1D("Nsphisto", "Nsphisto", 50, 200, 800);
	TGraph* lcchi = new TGraph();
	TGraph* plchi = new TGraph();

	for(Int_t i = 0; i<=iter; i++)
	{
		tmeanhisto->Fill(tmean[i]);
		chigaushisto->Fill(chigaus[i]);
		indhisto->Fill(ind[i]);
		chiplhisto->Fill(chipl[i]);
		Nlchisto->Fill(Nlc[i]);
		Nsphisto->Fill(Nsp[i]);
		lcchi->SetPoint(i, tmean[i], chigaus[i]);
		plchi->SetPoint(i, ind[i], chipl[i]);

		t = tmean[i];
		tree->Fill();

	}

	tree->Write();
	outputfile->Close();



	gStyle->SetOptStat(111111);

	/*TCanvas* big = new TCanvas("tau0sim", "tau0sim", 1000, 1000);
	big->Divide(2,4);
	big->cd(1);
	tmeanhisto->SetTitle("LC Tmax");
	tmeanhisto->Draw("EP");
	big->cd(2);
	chigaushisto->SetTitle("Chi2/NDF (LC)");
	chigaushisto->Draw("EP");
	big->cd(3);
	indhisto->SetTitle("PL index");
	indhisto->Draw("EP");
	big->cd(4);
	chiplhisto->SetTitle("Chi2/NDF (PL)");
	chiplhisto->Draw("EP");
	big->cd(5);
	lcchi->SetTitle("LC Tmax vs Chi2/NDF (LC)");
	lcchi->Draw("AP*");
	big->cd(6);
	plchi->SetTitle("PL index vs Chi2/NDF (PL)");
	plchi->Draw("AP*");
	big->cd(7);
	Nlchisto->SetTitle("Events LC");
	Nlchisto->Draw("EP");
	big->cd(8);
	Nsphisto->SetTitle("Events Fit");
	Nsphisto->Draw("EP");
	big->SaveAs("bigcanvas.pdf");*/



}



