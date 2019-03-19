const Int_t Np = 1491;   // Number of events (hardcoded)
const Int_t numbin = 500;
Double_t t[Np];
Double_t E[Np];
const Int_t numOfIte = 10000000;
Double_t ratio[numOfIte];
Int_t ind;

void Pairloop()
{

  //Loop to repeat for many flares
  const int numOfFlares = 100;
  Double_t results[numOfFlares];
  Int_t Npoints;
  
for(Int_t l=0; l<numOfFlares; l++)
{
  char *fname;
  ifstream in;
  Double_t v1,v2;
//   if(l==1)
//  {
  //Read file and fill the arrays t and E
  std::cout << "antes " << t[0] << " " << E[0] << std::endl;
  memset(t, 0, Np);
  memset(E, 0, Np);
  std::cout << "despues " << t[0] << " " << E[0] << std::endl;
  Npoints = 0;
  ind = 0;
  fname = Form("/home/lnogues/workspace_cpp/LIVanalysis/likelihood_code/Simulation/Mrk501/Manel_flares/alpha/with_Resol/Manel_flare_alpha-20_%i.txt", l); //Name of the simulations here.
  in.open(fname);
  std::cout << "Flare: " << fname << std::endl;

  while(1)
    {
	  in >> v1 >> v2;
	  if (!in.good()) break;
	  if(v1>1000){
		  t[Npoints] = v1;
		  E[Npoints] = v2;
		  Npoints++;
	  }
    }

  std::cout << "Npoints: " << Npoints-1 << std::endl;

  //Compute the ratio (deltaT/deltaE) and fill the array ratio[] with the results
  for ( Int_t i=0; i<Npoints; i++)
    {
    for( Int_t j=0; j<Npoints; j++)
      {
    	if(i>j)
    	{
    		ratio[ind]=((t[i]-t[j])/(E[i]-E[j]));
    		// std::cout << "ratio " << ratio[ind] << std::endl;
//    		if(ind ==0) std::cout << "ratio " << ratio[ind] << std::endl;
    		ind++;
    	}
      }
    }
  std::cout << "Numero total de componentes: " << ind << std::endl;

  //Compute mean ratio to define histogram
  Double_t sum_ind=0;
  Double_t mean_ind=0;
  for(Int_t i=0; i<ind; i++){
	  sum_ind+=ratio[i];
  }
  mean_ind = sum_ind/(ind-1);
  std::cout << "Sum ind value: " << sum_ind << std::endl;
  std::cout << "Mean ind value: " << mean_ind << std::endl;

  TH1D* histo = new TH1D("histo", "PairView", numbin, 0-TMath::Abs(mean_ind)*2, 0+TMath::Abs(mean_ind)*2);
  for ( Int_t a = 0; a < ind; a++)
    {
      histo->Fill(ratio[a]);
      histo->SetXTitle("l_{i,j}");
    }
//  TCanvas* hi = new TCanvas();
//  histo->Draw();

  //Look for the maximum in the histogram
  Double_t value = 0.;
  Double_t mid;

  for (Int_t bin=1; bin<=numbin; bin++)
    {
      Double_t ex = histo->GetBinContent(bin);
      if(ex>value)
      {
    	  value = ex;
    	  mid = histo->GetBinCenter(bin);
      }
      else continue;
    } 
  std::cout << "The maximum of the histogram is: " << value << " located at QG scale of: " << mid << std::endl;
  results[l] = mid;
  std::cout << "Result: " << results[l] << std::endl;
  
  histo->Delete();
  value = 0.;
  mid = 0.;
//  } //This should be removed for tests
}

//Computation of the final results and the error
// std::cout << "results[0]: " << results[0] << std::endl;
 Double_t mean=0;
 Double_t sum=0;
 Double_t dev = 0;
 Double_t error = 0.;

 TH1D* final = new TH1D("100 sim", "100 sim", 400, -6., 6.);

 for(int i=0; i<numOfFlares; i++){

	 sum+= results[i];
	 final->Fill(results[i]);
 }
 std::cout << "sum: " << sum << std::endl;
 mean = sum/numOfFlares;
 for(int i=0; i<numOfFlares; i++){
	 dev+= pow((results[i]-mean),2);
 }
 error = sqrt(dev/numOfFlares);
 std::cout << "Mean: " << mean << " Error: " << error << std::endl;
 gStyle->SetOptFit(1111);
 TCanvas* finalcanvas =  new TCanvas();
 final->Draw();
 final->Fit("gaus");
 final->GetFunction("gaus")->SetNpx(500);
 finalcanvas->Update();

}
