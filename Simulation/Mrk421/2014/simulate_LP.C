#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>

#include "MBinning.h"

const Double_t ln10 = 2.302585093;
const Double_t E0 = 364.;
const Int_t nbins = 20;
//
//
// if X is a random variable that is distributed like f(X),
// then a variable transformation y = g(X) yields:
// h(y) = f (g^-1(y)) * d(g^-1(y))/dy.
//
// In our case: 
// f(E) = (E/E_0) ^(a -b * (b * log10(E/E_0)))  = (E/E_0) ^(a - b^2/ln(10) * ln(E/E_0))  
// y = ln(E/E_0)
// g^-1(y) = e^y * E_0   ,   d(g^-1(y))/dy = E_0 * e^y
// 
// h(y) = (e^y)^(a - b^2/ln(10)*y) * E_0 * e^y = e^(y * (a - b^2/ln(10) * y)) * E_0 * e^y = E_0 * e^(y * (a+1. -b^2/ln(10) * y)) = E_0 * e^(y * (a + 1. -b^2/ln(10) * y))
//
TF1 *get_function(TString file)
{
  TFile *f = new TFile(file);
  TF1 *g = (TF1*)f->Get("SpectralModel");

  return g;
}

Double_t get_E(Double_t y)
{
  return exp(y)*E0;
}

void simulate_LP(TString file="Output_fold_LP.root", const Double_t emin=150., const Double_t emax=30000., const ULong_t N=1E6)
{

  TF1 *fold = get_function(file);
  TString logform = fold->GetTitle();
  cout << "Fold output: " << logform << endl;
  cout << "Fold parameters: " << fold->GetParameter(0) << ", " << fold->GetParameter(1) << ", " << fold->GetParameter(2) << endl;

  // check also the SED
  // TString sed(logform);
  // sed.Append(Form("*x/%.0f*x/%.0f",1.,1.));
  logform.ReplaceAll("[0]*","");  
  TString exponent(Form("([1]-([2]*([2]*log10(x/%.0f))))",E0));
  logform = Form("[0]+%s*log(x/%.1f)",exponent.Data(),E0);
  // sed = Form("log(%s)",sed.Data());
  // TF1 *fsed = new TF1("fsed",sed,emin,emax);
  // fsed->SetLineColor(kPink+1);
  // fsed->SetParameters(fold->GetParameter(0)*N,fold->GetParameter(1),fold->GetParameter(2));
  
  Double_t ymin = log(emin/E0);
  Double_t ymax = log(emax/E0);
  
  TString lognewform = (Form("%.8f+(x*(%f-%f*(%f*x)))",log(E0),fold->GetParameter(1),fold->GetParameter(2),fold->GetParameter(2)/ln10));
  cout << "New logarithmic form: ln(P(y))=" << lognewform << endl;
  TF1 *flognew = new TF1("flognew",lognewform,ymin,ymax);

  TCanvas *c = new TCanvas("c","c",1200,400);
  c->Divide(3,1);
  c->cd(1);
  fold->SetRange(emin,emax);
  fold->GetHistogram()->SetTitle("Original function from fold;E (GeV);flux");
  fold->Draw();
  gPad->Modified();
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->Update();
  
  TRandom3 *ran = new TRandom3();
  ran->SetSeed(0);

  // histogram for tests
  MBinning bins;
  bins.SetEdgesLog(nbins,emin,emax);
  
  TH1F *hist    = new TH1F("hist","simulated spectrum (linear scale);E (GeV);dN/dE",nbins,bins.GetEdges());
  hist->SetStats(0);
  TH1F *loghist = new TH1F("loghist","simulated spectrum (log scale);E (GeV);ln(dN/dE)",nbins,bins.GetEdges());
  loghist->SetStats(0);  
  //  TH1F *hsed    = new TH1F("sed","simulated SED",nbins,bins.GetEdges());


  // initialize simulation:
  Double_t deltay = flognew->Eval(ymin)-flognew->Eval(ymax);
  Double_t pmin   = flognew->Eval(ymax);

  //  cout << "pmin: " << pmin << " pmax: " << flognew->Eval(ymin) << endl;
	     
  ULong_t cnt=0;							  
  while (cnt<N)
    {
      Double_t y = ran->Uniform(ymin,ymax);
      Double_t p = ran->Uniform(deltay);
      
      if (flognew->Eval(y) < (log(p)-pmin))
	continue;

      hist->Fill(get_E(y));
      cnt++;
    }
  
  // Fill the test histogram loghist
  for (Int_t b=1;b<=hist->GetNbinsX();b++)
    {
      Double_t content    = hist->GetBinContent(b);
      Double_t logcontent = (content > 0. ? log(content) : -0.69);
      loghist->SetBinContent(b,logcontent);
      loghist->SetBinError(b,content > 1. ? 1./content*hist->GetBinError(b) : 1.41);
      // Double_t E = hist->GetXaxis()->GetBinCenter(b);
      // hsed->SetBinContent(b,content > 0. ? log(content * E * E) : log(0.5*E*E));
      // hsed->SetBinError  (b,loghist->GetBinError(b));
      // cout << "H" << loghist->GetBinContent(b) << " +- " << loghist->GetBinError(b) << endl;
      // cout << hsed->GetBinContent(b) << " +- " << hsed->GetBinError(b) << endl;
      // cout << fsed->Eval(E) << endl;
    }
  cout << endl;

  c->cd(2);

  // loghist->Draw();
  // return;

  cout << "Fitting function: " << logform << endl;

  TLegend *leg = new TLegend(0.15,0.2,0.4,0.45,"","brNDC");

  TF1 *flogfit = new TF1("flogfit",logform,emin,emax);
  flogfit->SetParLimits(1,-5.,-1.);
  flogfit->SetParLimits(2,0.,2.);
  // flogfit->FixParameter(1,fold->GetParameter(1));
  // flogfit->FixParameter(2,fold->GetParameter(2));
  flogfit->SetParameter(1,fold->GetParameter(1));
  flogfit->SetParameter(2,fold->GetParameter(2));
  flogfit->SetNpx(1000);
  flogfit->SetParameter(0,log(hist->GetEntries())/1.3);
  flogfit->SetLineColor(kOrange+2);  
  loghist->SetLineColor(kGreen+2);
  loghist->Fit(flogfit,"");

  leg->AddEntry(flogfit,"fitted","L");

  
  TF1 *flogfitold = (TF1*)flogfit->Clone("flogfitold");
  flogfitold->SetParameter(1,fold->GetParameter(1));
  flogfitold->SetParameter(2,fold->GetParameter(2));
  flogfitold->SetParameter(0,flogfit->GetParameter(0));  
  flogfitold->SetLineColor(kBlue);
  flogfitold->SetNpx(1000);
  flogfitold->Draw("same");

  leg->AddEntry(flogfitold,"orignal","L");
  leg->Draw();
  
  // hsed->SetLineColor(kCyan+1);
  // hsed->Draw("same");  
  // fsed->Draw("same");

  gPad->SetLogx();
  
  c->cd(3);
  
  TF1 *ffit = (TF1*)fold->Clone("ffit");
  ffit->SetParameter(1,flogfit->GetParameter(1));
  ffit->SetParameter(2,flogfit->GetParameter(2));
  ffit->SetParameter(0,exp(flogfit->GetParameter(0)));
  ffit->SetNpx(1000);

  // ffit->SetParameter(1,fold->GetParameter(1));
  // ffit->SetParameter(2,fold->GetParameter(2));
  fold->SetParameter(0,ffit->GetParameter(0));
  fold->SetLineColor(kBlue);
  //  hist->Fit(ffit);
  hist->Draw();
  ffit->SetLineColor(kOrange+2);
  ffit->Draw("same");
  fold->Draw("same");

  leg = new TLegend(0.15,0.2,0.4,0.45,"","brNDC");
  leg->AddEntry(ffit,"fitted"  ,"L");
  leg->AddEntry(fold,"original","L");
  leg->Draw();
  
  cout << "Test results: " << endl;
  cout << Form("p1:  %.5f (old)    %.5f+-%.5f (new)",fold->GetParameter(1),flogfit->GetParameter(1),flogfit->GetParError(1)) << endl;
  cout << Form("p2:  %.5f (old)    %.5f+-%.5f (new)",fold->GetParameter(2),flogfit->GetParameter(2),flogfit->GetParError(2)) << endl;

  gPad->SetLogy();
  gPad->SetLogx();

} 
