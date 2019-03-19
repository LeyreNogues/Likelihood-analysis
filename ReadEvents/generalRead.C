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

     Double_t mintime = 1.77293506944444380e+03; //From wobble info
     Double_t maxtime = 1.77308576388889196e+03; //From wobble info
     Double_t maxt_sec = 0.;
     for (int a = 1; a < Npoints; a++)
     {
    	 tp[a]=(MJD[a]-mintime)*24.*60.*60.;
//    	 tp[a]=(MJD[a]-1772.93512107)*24.*60.*60.; //For OFF case

    	 if(maxt_sec<tp[a])maxt_sec=tp[a];
//    	 if(tp[a]>12093.9073919988 && tp[a]< 12133.68376126309158){
//    		 cout << tp[a] << " is a in hole" << endl;
//    	 }
     }

     cout << "The minimum time is " << tp[0] << " and the maximum time is " << maxt_sec << " s"<< endl;


     //Create time histogram with wobbles and analysis cuts information.
     const Int_t binsnum = 500; //A big value
     Int_t numberOfBins; //The real one
     Double_t bin[binsnum+1];
     Double_t binwidth=60.; // In seconds(60)
     const Int_t numberOfWobbles = 14; //Number of wobbles
     const Int_t numberOfPoints = 2*numberOfWobbles;
     Double_t hole[numberOfPoints];

 	for(Int_t i=0;i<=binsnum;i++) bin[i]=-1.; //Initialize
 	hole[0]=0.; //Beggining of first wobble
 	hole[1] = (1.77294540509259241e+03 - mintime)*24.*60.*60.;
 	hole[2] = (1.77294579861110833e+03 - mintime)*24.*60.*60.; //Wobble2
 	hole[3] = (1.77295609953703388e+03 - mintime)*24.*60.*60.;
 	hole[4] = (1.77295817129629722e+03 - mintime)*24.*60.*60.; //Wobble3
 	hole[5] = (1.77296849537036906e+03 - mintime)*24.*60.*60.;
 	hole[6] = (1.77296908564814657e+03 - mintime)*24.*60.*60.; //Wobble4
 	hole[7] = (1.77297917824074102e+03 - mintime)*24.*60.*60.;
 	hole[8] = (1.77297956018518744e+03 - mintime)*24.*60.*60.; //Wobble 5
 	hole[9] = (1.77298983796295943e+03 - mintime)*24.*60.*60.;
 	hole[10] = (1.77299023148148262e+03 - mintime)*24.*60.*60.; //Wobble6
 	hole[11] = (1.77300049768518511e+03 - mintime)*24.*60.*60.;
 	hole[12] = (1.77300089120370103e+03 - mintime)*24.*60.*60.; //Wobble7
 	hole[13] = (1.77301114583333401e+03 - mintime)*24.*60.*60.;
 	hole[14] = (1.77301156250000349e+03 - mintime)*24.*60.*60.; //Wobble8
 	hole[15] = (1.77302184027777548e+03 - mintime)*24.*60.*60.;
 	hole[16] = (1.77302224537036818e+03 - mintime)*24.*60.*60.; //Wobble9
 	hole[17] = (1.77303248842592438e+03 - mintime)*24.*60.*60.;
 	hole[18] = (1.77303289351851708e+03 - mintime)*24.*60.*60.; //Wobble10
 	hole[19] = (1.77304315972221957e+03 - mintime)*24.*60.*60.;
 	hole[20] = (1.77304371527778130e+03 - mintime)*24.*60.*60.; //Wobble11
 	hole[21] = (1.77305380787036847e+03 - mintime)*24.*60.*60.;
 	hole[22] = (1.77305427083333052e+03 - mintime)*24.*60.*60.; //Wobble12
 	hole[23] = (1.77306445601851738e+03 - mintime)*24.*60.*60.;
 	hole[24] = (1.77306501157407183e+03 - mintime)*24.*60.*60.; //Wobble 13
 	hole[25] = (1.77307512731481256e+03 - mintime)*24.*60.*60.;
 	hole[26] = (1.77307554398148204e+03 - mintime)*24.*60.*60.; //Wobble14
 	hole[27] = (maxtime - mintime)*24.*60.*60.; //End of last wobble (not last event)

 	cout << "hole[1]: " << hole[1] << endl;



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


	  //-------------------------Very high energy event list (for camera display)---------------------------

	  	  Int_t VHE_selected = 0;
	  	  FILE *fout2;
	  	  fout2 = fopen("Mrk421_Paper_ONVHEevents4_Zd50.txt","wb");
	  	  for(int i = 0; i<Npoints; i++)
	  	  {
	  		  if(Ep[i]>4000.){
	  			  fprintf(fout2,"%f    %f   \n", tp[i], Ep[i]);
	  			VHE_selected++;
	  		  }
	  	  }
	  	  cout << "Selected VHE events: " << VHE_selected << endl;

      //-----------------------------Create event list for Likelihood---------------------------

	  Int_t selected = 0;
	  FILE *fout;
	  fout = fopen("Event_Paper_Zd50_ON.txt","wb");
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
