//TF1* const hk = new TF1();
TF1* hk;
Double_t Scaling(Double_t* x, Double_t* par)
{
  return  hk->Eval(x[0])/par[0];
}

const Int_t Np = 1000;   // Number of events (hardcoded)                                                                                                                                                                                    
Double_t xp[Np];
Double_t yp[Np];
Int_t comp = 10000000000;

//Double_t ratio[10000000000];
Int_t ind = 0;

void KernelPair(char *fname)
{
  std::vector<double> ratio(comp);
  //Read file and fill the arrays xp and yp
  ifstream in;
  in.open(fname);
  Double_t v1,v2;
  Int_t Npoints = 0;
  while(1)
    {
      in >> v1 >> v2;
      if (!in.good()) break;
      xp[Npoints] = v1;//time
      yp[Npoints] = v2;//energy
      Npoints++;
    }

  //Compute the ratio (deltaT/deltaE) and fill the array ratio[] with the results
  for ( Int_t i=0; i<Np; i++)
    {
    for( Int_t j=0; j<Np; j++)
      {
	if(i>j)
	  {
	    ratio[ind]=((xp[i]-xp[j])/(yp[i]-yp[j]));
	    ind++;
	  }
      }
    }
  cout << "Numero total de componentes: " << ind << endl;
  
  //Get maximum and miminum value in the ratio[] array

  Double_t max, min;
  max = min = ratio[0];
  for ( Int_t a = 1; a < ind; a++)
    {
      if(ratio[a]<min)
	{
	  min=ratio[a];
	}
      if(ratio[a]>max)
	{
	  max=ratio[a];
	}
    }
  cout << "The minimum is " << min << " and the maximum is " << max << endl;

  //Fill and histogram with the results
  TH1D* histo = new TH1D("histo", "PairView", 200, -4., 4.);
  for ( Int_t a = 0; a < ind; a++)
    {
      histo->Fill(ratio[a]);
      histo->SetXTitle("l_{i,j}");
    }

  TCanvas* canvas = new TCanvas();
  canvas->SetTitle("TKDE test");
  histo->Scale(1./histo->Integral(), "width");
  histo->SetStats(false);
  // histo->SetTitle("l_{i,j} distribution");
  histo->Draw("");

  //Create the TKDE Class (Kernel Density Estimation)
  
  cout << "Creating TKDE class" << endl; 
  double rho = 1.0; //default value                                                                                                  
  TKDE * kde = new TKDE(ind, &ratio[0], -4., 4., "", rho);
  cout << "TKDE finised" << endl;
  //kde->Draw("ConfidenceInterval@0.95 Same");                                                                                       
  //kde->Draw("SAME");

  hk = kde->GetFunction(1000);                                                                                        
  hk->SetLineColor(kRed);                                                                                                        
  //hk->Draw("SAME");

  //Normalize the TKDE function
  Double_t norm = hk->Integral(-4, 4.);
  cout << "Norm is: " << norm << endl;
 
  TF1* const hknorm = new TF1("khnorm", Scaling, -4., 4., 1); 
  hknorm->SetParameter(0, norm);
  hknorm->Draw("SAME");

  cout << "The maximum is located at " << hk->GetMaximumX() << endl;

  /*TF1 * fl = kde->GetLowerFunction(0.684);
  TF1 * fu = kde->GetUpperFunction(0.684);
  fl->Draw("SAME");
  fl->SetLineColor(kBlue-5);
  fu->SetLineColor(kBlue-5);
  fu->Draw("SAME");*/


  TLegend * legend = new TLegend(0.6,0.7,0.9,0.95);
  //legend->AddEntry(kde->GetDrawnFunction(),"TKDE");
  legend->AddEntry(hknorm,"TKDE");                                                                  
  //legend->AddEntry(fl,"TKDE - #sigma");                                                                  
  //legend->AddEntry(fu,"TKDE + #sigma");
  legend->Draw();
  // return;

  gPad->Update(); 
}
