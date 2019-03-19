/*
 * readevents.C
 *
 *  Created on: Feb 17, 2016
 *      Author: Leyre Nogués
 *
 *  This code allows to read a Event List (the one created by readEventsfromFlute.C)
 *  and allow to do cut, tranform MJD to second and create histograms.
 *  Finally it makes another final Event List for the Likelihood Code.
 *
 */

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
#include<iostream>
#include<fstream>
#include<istream>
#include<TVectorD.h>

using namespace std;


const Int_t Np = 100000;   // Number of events (hardcoded)
Double_t Ep[Np]; //In GeV
Double_t tp[Np]; //In s
Double_t MJD[Np]; //In MJD


void generalRead(char *fname)
{

//---------------------------------------Read the events--------------------------------------

	cout.precision(15);
	cout << "Reading Event List in: " << fname << endl;

  ifstream in;
  in.open(fname);
  Double_t v1,v2;
  int Npoints = 0;
  while(1)
    {
      in >> v1 >> v2;
      if (!in.good()) break;
      MJD[Npoints] = v1;
      Ep[Npoints] = v2;
      Npoints++;
    }

  cout << "Npoints(Number of events): " << Npoints << endl;


  //-----------------------------Maximum and minimum time and energy---------------------------

  //TIME
  Double_t maxt, mint;
  maxt = mint = MJD[0];
  cout << "MJD[0]: " << MJD[0] << endl;
  for ( int a = 1; a < Npoints; a++)
     {
       if(MJD[a]<mint)
 	{
           mint=MJD[a];
 	}
       if(MJD[a]>maxt)
         {
           maxt=MJD[a];
         }
     }
   cout << "The minimum time is " << mint << " and the maximum time is " << maxt << " MJD"<< endl;

   //ENERGY
     Double_t maxE, minE;
     maxE = minE = Ep[0];
     for ( int a = 1; a < Npoints; a++)
        {
          if(Ep[a]<minE)
    	{
              minE=Ep[a];
    	}
          if(Ep[a]>maxE)
            {
              maxE=Ep[a];
            }
        }

     cout << "The minimum energy is " << minE << " and the maximum energy is " << maxE << " GeV"<< endl;


      //-----------------------------Histograms for time and energy---------------------------

     //Cuts in Energy
     Double_t cutInE_Low = 120.;
     Double_t cutInE_High = 100000.;
     cout << "Cuts in energy will be: " << cutInE_Low << " GeV - " << cutInE_High << " GeV"<< endl;

     Double_t XX_step = 0.1;
     Int_t n = 0;
     if(cutInE_High<maxE) n = ceil((log10(cutInE_High)-log10(cutInE_Low))/XX_step);
     else n = ceil((log10(maxE)-log10(cutInE_Low))/XX_step);
     Double_t *XX = new double[n+1];
     for(int i=0;i<n+1;i++)
     {
    	 XX[i]= pow(10,(log10(cutInE_Low)+i*XX_step));
     }

     cout << "Log bins ready" << endl;

     Double_t time_bins = 200;
     TH1D* time = new TH1D("time", "time", time_bins, mint, maxt);
     TH1D* energy = new TH1D("energy", "energy", n, XX);
     energy->Sumw2();
     Int_t selected_events = 0;

     cout << "Histograms ready" << endl;

     for ( int a = 0; a < Npoints; a++)
     {
    	 if(Ep[a]>=cutInE_Low && Ep[a]<=cutInE_High)
    	 {

    		 time->Fill(MJD[a]);
    		 energy->Fill(Ep[a],1/Ep[a]);
    		 selected_events++;
    	 }

     }

     cout<< "Selected events: " << selected_events<<endl;

     TCanvas* dist = new TCanvas();
     dist->SetTitle("Energy and time distributions");
     dist->Divide(1,2);
     dist->cd(1);
     time->GetXaxis()->SetTitle("Time(Mjd)");
     time->Draw("EP");
     dist->cd(2);
     dist->cd(2)->SetLogx();
     dist->cd(2)->SetLogy();
     energy->GetXaxis()->SetTitle("Energy(GeV)");
     energy->Draw("EP");




      //-----------------------------From MJD to seconds---------------------------

     Double_t mintime = 3497.887615740743058; //From wobble info
     Double_t maxtime = 3498.057141203702486; //From wobble info
     Double_t maxt_sec = 0.;
     for (int a = 1; a < Npoints; a++)
     {
    	 tp[a]=(MJD[a]-MJD[0])*24.*60.*60.;
    	 if(maxt_sec<tp[a])maxt_sec=tp[a];
     }

     cout << "The minimum time is " << tp[0] << " and the maximum time is " << maxt_sec << " s"<< endl;


     //Create time histogram with wobbles and analysis cuts information.
     const Int_t binsnum = 500; //A big value
     Int_t numberOfBins; //The real one
     Double_t bin[binsnum+1];
     Double_t binwidth=120.; // In seconds(60)
     const Int_t numberOfWobbles = 3; //Number of wobbles
     const Int_t numberOfPoints = 2*numberOfWobbles;
     Double_t hole[numberOfPoints];

 	for(Int_t i=0;i<=binsnum;i++) bin[i]=-1.; //Initialize
 	hole[0]=0.; //Beggining of first wobble
	hole[1] = (3497.901354166664532 - mintime)*24.*60.*60.;//End of W1
	hole[2] = (3497.901782407410792 - mintime)*24.*60.*60.; //Wobble2
	hole[3] = (3497.915520833332266 - mintime)*24.*60.*60.; //End of W2
	hole[4] = (3497.91594907407125 - mintime)*24.*60.*60.; //Wobble3
	hole[5] = (3497.92969907407678 - mintime)*24.*60.*60.;
	hole[6] = (3497.930127314815763 - mintime)*24.*60.*60.; //Wobble4
	hole[7] = (3497.943865740737238 - mintime)*24.*60.*60.;
	hole[8] = (3497.944305555553001 - mintime)*24.*60.*60.; //Wobble 5
	hole[9] = (3497.958032407404971 - mintime)*24.*60.*60.;
	hole[10] = (3497.958437499997672 - mintime)*24.*60.*60.; //Wobble6
	hole[11] = (3497.972210648149485 - mintime)*24.*60.*60.;
	hole[12] = (3497.972673611111531 - mintime)*24.*60.*60.; //Wobble7
	hole[13] = (3497.986400462963502 - mintime)*24.*60.*60.;
	hole[14] = (3497.986805555556202 - mintime)*24.*60.*60.; //Wobble8
	hole[15] = (3498.000543981484952 - mintime)*24.*60.*60.;
	hole[16] = (3498.00098379629344 - mintime)*24.*60.*60.; //Wobble9
	hole[17] = (3498.014687499999127 - mintime)*24.*60.*60.;
	hole[18] = (3498.015081018515048 - mintime)*24.*60.*60.; //Wobble10
	hole[19] = (3498.028842592590081 - mintime)*24.*60.*60.;
	hole[20] = (3498.02944444444438 - mintime)*24.*60.*60.; //Wobble11
	hole[21] = (3498.042986111111531 - mintime)*24.*60.*60.;
	hole[22] = (3498.043414351850515 - mintime)*24.*60.*60.; //Wobble12
	hole[23] = (maxtime - mintime)*24.*60.*60.; //End of last wobble (not last event)

	cout << "hole[1]: " << hole[1] << endl;
	cout << "hole[23]: " << hole[23] << endl;




 	//Divide every interval and fill the bins.
 	const Int_t numberOfIntervals = numberOfWobbles + (numberOfWobbles-1);
 	Int_t I_bins[numberOfIntervals];
 	Int_t currentBinIndex=0;
 	Int_t lastBinIndex = 0;


 	for(Int_t j = 0; j <numberOfIntervals; j+=2){
 		I_bins[j] = TMath::Nint((hole[j+1]-hole[j])/binwidth);
 		cout << "number of bins in I" << j << ": " << I_bins[j] << endl;

 		if(I_bins[j]<1){
 //			cout << "no more bins" << endl;
 			bin[lastBinIndex] = maxtime;
 			break;
 		}

 			for (Int_t i=0; i <I_bins[j]; i++){
 				bin[currentBinIndex] = hole[j];
 				bin[i+currentBinIndex]= bin[currentBinIndex]+((hole[j+1]-hole[j])/I_bins[j])*i;
// 				cout << "bin[" << i+currentBinIndex << "]: " << bin[i+currentBinIndex] << endl;
 				if(i==I_bins[j]-1){
 					bin[i+currentBinIndex+1] = hole[j+1];
// 					cout << "bin[" << i+currentBinIndex+1 << "]: " << bin[i+currentBinIndex+1] << endl;
 				}
 			}
 			currentBinIndex+= I_bins[j]+1;
// 			cout << "currentBinIndex: " << currentBinIndex << endl;
 			lastBinIndex = currentBinIndex-1;

 	}

 	cout << "lastBinIndex: " << lastBinIndex << endl;
 	numberOfBins = lastBinIndex-1;


 	TH1D* time_insec = new TH1D("Lightcurve", "Lightcurve", numberOfBins, bin);
 	time_insec->SetStats(0);


	  for ( int a = 0; a < Npoints; a++)
	  {
		  if(Ep[a]>=cutInE_Low && Ep[a]<=cutInE_High)
		  {
			  time_insec->Fill(tp[a]);
		  }
	  }

	  TCanvas* selectedtimecanvas = new TCanvas();
	  time_insec->GetXaxis()->SetTitle("Time(s)");
	  time_insec->GetYaxis()->SetTitle("Selected events");
	  time_insec->Draw();

	  TCanvas* selectedenergycanvas = new TCanvas();
	  selectedenergycanvas->SetLogx();
	  selectedenergycanvas->SetLogy();
	  energy->SetTitle("Spectrum");
	  energy->SetStats(0);
	  energy->GetXaxis()->SetTitle("Energy(GeV)");
	  energy->GetYaxis()->SetTitle("Selected events");
	  energy->Draw("EP");

      //-----------------------------Create event list for Likelihood---------------------------

	  Int_t selected = 0;
	  FILE *fout;
	  fout = fopen("Events_GRB_ON_Koji_seconds.txt","wb");
	  for(int i = 0; i<Npoints; i++)
	  {
		  if(Ep[i]>=cutInE_Low && Ep[i]<=cutInE_High){
			  fprintf(fout,"%f    %f   \n", tp[i], Ep[i]);
			  selected++;
		  }
	  }
	  cout << "Selected events: " << selected << endl;



      //-----------------------------Create root files---------------------------


//	  	  TFile* ON_outputfile = new TFile("Histo_25OFF_Kazuma_120.root","recreate");
//	  	  time_insec->SetName("LC_insecs");
//	  	  time_insec->Write();
//	  	  energy->SetName("Energy_dist");
//	  	  energy->Write();
//	  	  ON_outputfile->Close();


}
