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
#include<TKDE.h>
#include<TTree.h>
#include <vector>
#include <iostream>
#include <numeric>

using namespace std;

const Int_t Np = 12520;   // Number of events (hardcoded)
std::vector<double> E;
std::vector<double> t;
Long64_t ncomponents = 0;
Int_t QGorder = 2; //Linear or Quadratic


void PairView(TString fname)
{
	cout.precision(10);
	Long64_t MAX = Np*(Np-1)/2;
	cout << "Number of pairs: " << MAX << endl;
	std::vector<double> ratio(MAX);
	ratio.clear();

	ifstream in;
	in.open(fname);
	Double_t v1,v2;
	while(1)
	{
		in >> v1 >> v2;
		if (in.eof()) break;
		t.push_back(v1);
		E.push_back(v2);
	}

  std::cout << "Number of events: " << E.size() << std::endl;
  //Remove last event
  E.erase(E.begin()+E.size()-1);
  std::cout << "New E size: " << E.size() << std::endl;


  Int_t Npoints = E.size();

  std::cout << "t[0]: " << t[0] << std::endl;
  std::cout << "t[1]: " << t[1] << std::endl;
  std::cout << "t[2]: " << t[2] << std::endl;
  std::cout << "t[Npoints-4]: " << t[Npoints-4] << std::endl;
  std::cout << "t[Npoints-3]: " << t[Npoints-3] << std::endl;
  std::cout << "t[Npoints-2]: " << t[E.size()-2] << std::endl;
  std::cout << "t[Npoints-1]: " << t[E.size()-1] << std::endl;
  std::cout << "E[0]: " << E[0] << std::endl;
  std::cout << "E[1]: " << E[1] << std::endl;
  std::cout << "E[2]: " << E[2] << std::endl;
  std::cout << "E[Npoints-4]: " << E[Npoints-4] << std::endl;
  std::cout << "E[Npoints-3]: " << E[Npoints-3] << std::endl;
  std::cout << "E[Npoints-2]: " << E.at(Npoints-2) << std::endl;
  std::cout << "E[Npoints-1]: " << E.at(Npoints-1) << std::endl;


  //Compute the ratio (deltaT/deltaE) and fill the array ratio[] with the results
  cout << "QGorder: " << QGorder << endl;
  for ( Int_t i=0; i<Npoints; i++)//Npoints
    {
    for( Int_t j=0; j<Npoints; j++)
      {
    	if(i>j)
    	{
    		if(E.at(i)==E.at(j)) continue;
    		if(QGorder == 1)ratio.push_back((t.at(i)-t.at(j))/(E.at(i)-E.at(j)));
    		if(QGorder == 2)ratio.push_back((t.at(i)-t.at(j))/((E.at(i)*E.at(i))-(E.at(j)*E.at(j))));
    		if(ratio[ncomponents]>10000000000){
    		std::cout << "ratio " << ratio[ncomponents] << std::endl;
    		std::cout << "E[i]: " << E.at(i) << " E[j]: " << E.at(j) << std::endl;
    		std::cout << "t[i]: " << t.at(i) << " t[j]: " << t.at(j) << std::endl;
    		}
    		ncomponents++;
    	 }
      }
    }

  std::cout << "Numero total de componentes: " << ncomponents << std::endl;
  std::cout << "Size of ratio: " << ratio.size() << std::endl;

  Double_t average = accumulate( ratio.begin(), ratio.end(), 0.0/ ratio.size());
  cout << "The average of ratio is: " << average << endl;

  //Get maximum and miminum value in the ratio[] array

  Double_t max = 0., min = 100.;
  for ( Int_t a = 0; a < ratio.size(); a++)
    {
      if(ratio.at(a)<min)
      {
    	  min=ratio.at(a);
      }
	  if(ratio.at(a)>max)
	  {
		  max=ratio.at(a);
	  }
    }

  cout << "The minimum is " << min << " and the maximum is " << max << endl;
  cout << "And the middle point is: " << min + (max-min)/2. << endl;

//  return;

  Int_t numOfBins = 1000;

  //Fill and histogram with the results
  Double_t up_limit = 0., down_limit = 0.;
  if(QGorder == 1){
	  up_limit = 0. + TMath::Abs(average)*0.00000005;
	  down_limit = 0. - TMath::Abs(average)*0.00000005;
  }
  if(QGorder == 2){
	  up_limit = 0 + TMath::Abs(average)*0.00000005;
	  down_limit = 0 - TMath::Abs(average)*0.00000005;
  }
  cout << "down_limits: " << down_limit << " up_limit: " << up_limit << endl;
  TH1D* histo = new TH1D("histo", "PairView", numOfBins, down_limit, up_limit);

  for (Long64_t a = 0; a < ncomponents; a++){
      histo->Fill(ratio.at(a));
    }

  TCanvas* histocanvas =  new TCanvas();
  histo->SetXTitle("l_{i,j}");
  histo->SetYTitle("Number of event pairs");
  histo->SetStats(0);
  histo->Draw();

  //Look for the maximum in the histogram
  Double_t value = 0.;
  Double_t mid = 0.;

  for (Int_t bin=1; bin<=numOfBins; bin++)
    {
      Double_t ex = histo->GetBinContent(bin);
      if(ex>value)
        {
          value = ex;
          mid = histo->GetBinCenter(bin);
        }
      else continue;
    }
  std::cout << "The maximum of the histogram is: " << value << " located at tau value: " << mid << std::endl;

	// Create a TTree to save results
  Double_t tau;
  TTree tree("t1","Tree with simple variables");
  tree.Branch("tau",&tau,"px/D");
  cout << "Tree created" << endl;
  tau = mid;

  tree.Fill();
  //Create de TFile
//  TFile* Likelihood_outputfile = new TFile("PairView_KazumaRandom_Quadratic_DATA.root","recreate");
//  histo->Write();
//  tree.Write();

  ratio.clear();
  ratio.resize(0);
  cout << "Ration size: " << ratio.size() << endl;
  t.clear();
  t.resize(0);
  E.clear();
  E.resize(0);
  ncomponents=0;


  cout << "Root file created"  << endl;

}
