/*
 * BinningTest.C
 *
 *  Created on: Feb 17, 2016
 *      Author: Leyre Nogués
 *
 *  This code allows to read a Event List (the one created by readEventsfromFlute.C)
 *  and allow to do cut, tranform MJD to second and create histograms.
 *  Finally it makes another final Event List for the Likelihood Code.
 *  UPDATE: This code also produces a histogram with fixed number of events per bin and
 *  its correspoding histogram showing dN/dt.
 *  It also saves the energy distribution of events (useful for the OFF region)
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
#include<TLegend.h>

using namespace std;


const Int_t Np = 100000;   // Number of events (hardcoded)
Double_t Ep[Np]; //In GeV
Double_t tp[Np]; //In s
Double_t MJD[Np]; //In MJD


void BinningTest(char *fname)
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
     Double_t cutInE_High = 200000;
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
//     TH1D* energy = new TH1D("energy", "energy", n, XX);
     TH1D* energy_original = new TH1D("energy_original", "energy_original", 2000, minE, maxE);
     TH1D* energy = new TH1D("energy", "energy", 2000, minE, maxE);


     energy->Sumw2();
     Int_t selected_events = 0;

     cout << "Histograms ready" << endl;

     for ( int a = 0; a < Npoints; a++)
     {
    	 if(Ep[a]>=cutInE_Low && Ep[a]<=cutInE_High)
    	 {

    		 time->Fill(MJD[a]);
//    		 energy->Fill(Ep[a],1/Ep[a]);
    		 energy_original->Fill(Ep[a]);
    		 selected_events++;
    	 }

     }

     for(Int_t i=1; i<=1000; i++){
    	 energy->SetBinContent(i, energy_original->GetBinContent(i)/3.);
     }

     cout<< "Selected events: " << selected_events<<endl;
     cout << "Integral of energy: " << energy->Integral();

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
     const Int_t binsnum = 2000; //A big value
     Int_t numberOfBins; //The real one
     Double_t bin[binsnum+1] = {0};
     Double_t binwidth=60.; // In seconds(60)
     const Int_t numberOfWobbles = 14; //Number of wobbles
     const Int_t numberOfPoints = 2*numberOfWobbles;
     Double_t hole[numberOfPoints] = {0};

 	for(Int_t i=0;i<=binsnum;i++) bin[i]=-1.; //Initialize
 	hole[0]=0.; //Beggining of first wobble
 	hole[1] = (1.77294540509259241e+03 - mintime)*24.*60.*60.;//End of W1
 	hole[2] = (1.77294579861110833e+03 - mintime)*24.*60.*60.; //Wobble2
 	hole[3] = (1.77295609953703388e+03 - mintime)*24.*60.*60.; //End of W2
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
 	cout << "hole[27]: " << hole[27] << endl;


 	//Divide every interval and fill the bins.
 	const Int_t numberOfIntervals = numberOfWobbles + (numberOfWobbles-1);
 	Int_t I_bins[numberOfIntervals] = {0};
 	Int_t currentBinIndex=0;
 	Int_t lastBinIndex = 0;
 	cout << "numberOfIntervals: " << numberOfIntervals << endl;


 	for(Int_t j = 0; j <numberOfIntervals; j+=2){
 		I_bins[j] = TMath::Nint((hole[j+1]-hole[j])/binwidth);
// 		cout << "number of bins in I" << j << ": " << I_bins[j] << endl;

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
 			lastBinIndex = currentBinIndex;

 	}

// 	cout << "lastBinIndex: " << lastBinIndex << endl;
 	numberOfBins = lastBinIndex-1;
 	cout << "Number of bins in LightCurve: " << numberOfBins << endl;


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
	  time_insec->GetYaxis()->SetTitle("Number of events");
	  time_insec->Draw();

	  TCanvas* selectedenergycanvas = new TCanvas();
	  selectedenergycanvas->SetLogx();
	  selectedenergycanvas->SetLogy();
	  energy->SetStats(0);
	  energy->SetTitle("Spectrum");
	  energy->GetXaxis()->SetTitle("Energy(GeV)");
	  energy->GetYaxis()->SetTitle("Number of events");
	  energy->Draw("EP");

	  //-----------------------------Order events in time---------------------------
	  cout << "Now I order the events in increasing time..."  << endl;
	  Double_t temp_t;
	  Double_t temp_E;
	  for(Int_t i=0;i<Npoints;i++)
	  {
//		  cout << i << endl;
		  for(Int_t j=0;j<Npoints-i;j++)
		  {
			  if(tp[j]>tp[j+1])
			  {
				  temp_t = tp[j];
				  tp[j] = tp[j+1];
				  tp[j+1]=temp_t;
				  temp_E = Ep[j];
				  Ep[j] = Ep[j+1];
				  Ep[j+1]=temp_E;
			  }
		  }
	  }

	  cout << "First event: " << tp[0] << " last event: " << tp[Npoints-1] << endl;

	  //-----------------------------Cut Events in Energy---------------------------
	  Double_t E[14500];
	  Double_t t[14500];
	  Int_t insideEvents = 0;


	  for(Int_t i=0;i<Npoints;i++){
		  if(Ep[i]>=cutInE_Low && Ep[i]<=cutInE_High){
			  E[insideEvents] = Ep[i];
			  t[insideEvents] = tp[i];
			  insideEvents++;

		  }
	  }

	  cout << "Inside: " << insideEvents << endl;


	  //-----------------------------Time bins for fixed - ON events---------------------------

	  for(Int_t i = 0; i<numberOfIntervals; i++){
		  cout << "I_bins: " << I_bins[i] << " ";
	  }
	  cout << endl;
	  cout << "numberOfIntervals: " << numberOfIntervals << endl;

	  Int_t counter = 0;
	  Int_t event_index = 0;
	  Int_t ON_events = 36;  //HERE YOU INDICATE THE NUMBER OF EVENTS PER BIN
	  cout << "ON_events: " << ON_events<< endl;
	  Int_t temp_ONevents = ON_events;
	  std::vector<double> miniBins;
	  Int_t bin_counter = 0;

	  //First bin
	  miniBins.push_back(hole[0]);

	  //Loop over wobbles
	  Int_t numberOfIntraBins = 0;
	  Int_t rest= 0;
	  for(Int_t a = 0; a<numberOfIntervals; a=a+2){//numberOfIntervals
		  ON_events = temp_ONevents;
		  Int_t eventsInWobble = 0;
		  cout << "---------------Wobble-----------------: " << a << endl;
		  cout << "Number of bins in the wobble: " << I_bins[a] << endl;
		  for (Int_t i=bin_counter+1; i <=bin_counter+I_bins[a]; i++){
//			  cout << "Content in bin: " << i << " is " << time_insec->GetBinContent(i) << endl;
			  eventsInWobble+=time_insec->GetBinContent(i);

		  }
		  cout << "The number of events in the the wobble is: " << eventsInWobble << endl;
		  numberOfIntraBins = TMath::Nint(eventsInWobble/ON_events);
		  rest = eventsInWobble%ON_events;
		  cout << "Cociente: " << numberOfIntraBins << endl;
		  cout << "Resto: " << rest << endl;
		  if(numberOfIntraBins<=rest){
			  cout << "Resto>Cociente" << endl;
			if(numberOfIntraBins == 0){
				  ON_events = eventsInWobble;
				  numberOfIntraBins = TMath::Nint(eventsInWobble/ON_events);
				  rest = eventsInWobble%ON_events;
			  }
			else{
			  ON_events = TMath::Nint(eventsInWobble/numberOfIntraBins);
			  numberOfIntraBins = TMath::Nint(eventsInWobble/ON_events);
			  rest = eventsInWobble%ON_events;
			}

			  cout << "New ON: " << ON_events << endl;
			  cout << "Cociente: " << numberOfIntraBins << endl;
			  cout << "Resto: " << rest << endl;

		  }
		  eventsInWobble+=time_insec->GetBinContent(bin_counter+I_bins[a]+1);
		  cout << "Events adding one more (empty) bin: " << eventsInWobble << endl;
//		  cout << "Wobble gap beggining: " << time_insec->GetBinLowEdge(bin_counter+I_bins[a]+1) << endl;
//		  cout << "Wobble gap end: " << time_insec->GetBinLowEdge(bin_counter+I_bins[a]+2) << endl;
		  bin_counter+=I_bins[a]+1;

		  //Loop over events inside the wobble
		  Int_t points_9 = 0;
		  Int_t points_10 = 0;
//		  cout << "event_index before loop: " << event_index << endl;
		  for(Int_t j=event_index;j<event_index+eventsInWobble;j++){
//				  cout << "counter: " << counter << endl;
				  if(points_10<rest){
					  if(counter==ON_events+1){
//					  	cout << "Point for array of 10 events" << endl;
//					  	cout << "event index: " << event_index << endl;
//					  	cout << "t event: " << t[j] << endl;
						  miniBins.push_back(t[j]);
						  points_10++;
						  counter = 0;
					  }
				  }
				  else{
					  if(counter==ON_events){
//					  	cout << "Point for array of 9 events" << endl;
//					  	cout << "event index: " << event_index << endl;
						  if(points_10+points_9 ==numberOfIntraBins-1){
							  miniBins.push_back(hole[a+1]);
							  counter=0;
//							  cout << "event index before break: " << event_index << endl;

							  break;
						  }
						  miniBins.push_back(t[j]);
//						  cout << "t event: " << t[j] << endl;
						  points_9++;
						  counter = 0;

					  }
				  }
				  counter++;
				  event_index++; //Inside the energy range
//				  cout << "Counter: " << counter << endl;
//				  cout << "event_index: " << event_index << endl;
			  }

		  	  if(a<numberOfIntervals-2)miniBins.push_back(hole[a+2]);
		  }

//	  for(Int_t i = 0; i<miniBins.size(); i++){
//		  cout << miniBins.at(i) << " ";
//	  }


	  cout << "Last event index: " << event_index << endl;
	  cout << "Size of miniBins: " << miniBins.size() << endl;
	  cout << "Number of miniBins: " << miniBins.size()-1 << endl;
	  Double_t* miniBins_pointer = miniBins.data();
	  cout << "Last component: " << miniBins.at(miniBins.size()-1) << endl;
	  cout << "Pre-Last component: " << miniBins.at(miniBins.size()-2) << endl;



	  //Draw histogram to see if it is correct.
	 TH1D* time_ONevents = new TH1D("time_ONevents", "time_ONevents", miniBins.size()-1, miniBins_pointer);
	 time_ONevents->SetStats(0);
	 TH1D* time_ONvariation = new TH1D("time_ONvariation", "time_ONvariation", miniBins.size()-1, miniBins_pointer);

	 Double_t min_binwidth = 1000.;
	 for(Int_t i=0; i<insideEvents;i++){
	  time_ONevents->Fill(t[i]);
	 }
	 for(Int_t i = 1; i<miniBins.size(); i++){
		 time_ONvariation->SetBinContent(i, time_ONevents->GetBinContent(i)/time_ONevents->GetBinWidth(i));
		 if(time_ONevents->GetBinWidth(i)<min_binwidth)min_binwidth=time_ONevents->GetBinWidth(i);
	 }

	 cout << "Miminum binwidth: " << min_binwidth << endl;

	 TCanvas* patata = new TCanvas();
	 patata->Divide(2,1);
	 patata->cd(1);
	 time_ONevents->SetTitle("");
	 time_ONevents->GetXaxis()->SetTitle("Time(s)");
	 time_ONevents->GetYaxis()->SetTitle("Events");
	 time_ONevents->Draw();
	 patata->cd(2);
	 time_ONvariation->SetStats(0);
	 time_ONvariation->SetTitle("");
	 time_ONvariation->GetXaxis()->SetTitle("Time(s)");
	 time_ONvariation->GetYaxis()->SetTitle("dN/dt");
	 time_ONvariation->Draw();


   //-------------------------------Check smoothness of the histogram---------------------------------

//	Int_t number_bins = 0;
//	Double_t difference =  0.;
//	for(Int_t i = 1; i<miniBins.size(); i++){
//		if(time_ONvariation->GetBinContent(i+1)!=0 && time_ONvariation->GetBinContent(i)!=0){
//			difference += TMath::Abs(time_ONvariation->GetBinContent(i+1) - time_ONvariation->GetBinContent(i));
//			number_bins++;
//		}
//		else cout << "Empty bin!" << endl;
//	}
//	cout << "Diference: " << difference << endl;
//	cout << "Bins: " << number_bins << endl;
//	Double_t Se = difference/number_bins;
//	cout << "Smooth estimator: " << Se << endl;


  //-----------------------------Compute rate of time-fixed histogram---------------------------

//	TH1D* time_insec_rate = new TH1D("Lightcurve_rate", "Lightcurve_rate", numberOfBins, bin);
//	for(Int_t i=1; i<=numberOfBins; i++){
////		cout << "Index" << i << endl;
////		cout << "BIn: " << bin[i-1] << endl;
////
////		cout << "Point: " << time_insec->GetBinContent(i)/time_insec->GetBinWidth(i) << endl;
//		time_insec_rate->SetBinContent(i, time_insec->GetBinContent(i)/time_insec->GetBinWidth(i));
//	}
//
//	 TCanvas* cebolla = new TCanvas();
//	 time_ONvariation->GetXaxis()->SetTitle("Time(s)");
//	 time_ONvariation->GetYaxis()->SetTitle("dN/dt");
//	 time_ONvariation->SetLineColor(4);
//	 time_ONvariation->Draw();
//	 time_insec_rate->SetLineColor(2);
//	 time_insec_rate->Draw("same");
//
//	TLegend* const2legend = new TLegend(0.1,0.7,0.4,0.9);
//	const2legend->AddEntry(time_ONvariation, "Fixed events", "l");
//	const2legend->AddEntry(time_insec_rate, "Fixed time", "l");
//	const2legend->Draw("");

  //-----------------------------Compute CDFs of both histograms---------------------------
//	TGraph* CDFs_variation =  new TGraph(miniBins.size());
//	TGraph* CDFs_time =  new TGraph(numberOfBins+1);
//	TH1D* time_insec_rate_copy = (TH1D*)time_insec_rate->Clone("Lightcurve_rate_copy");
//	time_insec_rate_copy->Scale(1./time_insec_rate_copy->Integral("width"));
//	TH1D* time_ONvariation_copy = (TH1D*)time_ONvariation->Clone("time_ONvariation_copy");
//	time_ONvariation_copy->Scale(1./time_ONvariation_copy->Integral("width"));
//
//	 TCanvas* ajo = new TCanvas();
//	 time_ONvariation_copy->GetXaxis()->SetTitle("Time(s)");
//	 time_ONvariation_copy->GetYaxis()->SetTitle("dN/dt");
//	 time_ONvariation_copy->SetLineColor(4);
//	 time_ONvariation_copy->Draw();
//	 time_insec_rate_copy->SetLineColor(2);
//	 time_insec_rate_copy->Draw("same");
//
//	 cout << "Integral variation: " << time_ONvariation_copy->Integral("width") << endl;
//	 cout << "Integral rate: " << time_insec_rate_copy->Integral("width") << endl;
//
//	 Double_t diff_variation = 0.;
//	 CDFs_variation->SetPoint(0,0.,0.);
//	 for(Int_t i=1; i<=miniBins.size()-1;i++){
////		 cout << "Bin center: " << time_ONvariation_copy->GetBinCenter(i) << endl;
//		 diff_variation+=time_ONvariation_copy->GetBinContent(i)*time_ONvariation_copy->GetBinWidth(i);
//		 CDFs_variation->SetPoint(i,time_ONvariation_copy->GetBinCenter(i),diff_variation);
//	 }
//
//	 Double_t diff_rate = 0.;
//	 CDFs_time->SetPoint(0,0.,0.);
//	 for(Int_t i=1; i<=numberOfBins;i++){
////		 cout << "Bin center: " << time_insec_rate_copy->GetBinCenter(i) << endl;
//		 diff_rate+=time_insec_rate_copy->GetBinContent(i)*time_insec_rate_copy->GetBinWidth(i);
//		 CDFs_time->SetPoint(i,time_insec_rate_copy->GetBinCenter(i),diff_rate);
//	 }
//
//	 TCanvas* berenjena = new TCanvas();
//	 CDFs_variation->GetXaxis()->SetTitle("Time(s)");
//	 CDFs_variation->GetYaxis()->SetTitle("CDF");
//	 CDFs_variation->SetMarkerColor(4);
//	 CDFs_variation->SetMarkerStyle(33);
//	 CDFs_variation->Draw("APL");
//	 CDFs_time->SetMarkerColor(2);
//	 CDFs_time->SetMarkerStyle(34);
//	 CDFs_time->Draw("PLsame");
//
//	TLegend* const3legend = new TLegend(0.1,0.7,0.4,0.9);
//	const3legend->AddEntry(CDFs_variation, "CDF Fixed events", "p");
//	const3legend->AddEntry(CDFs_time, "CDF Fixed time", "p");
//	const3legend->Draw("");
//
//	 Double_t max_diff = 0.;
//	 Int_t number_of_points = 1000;
//	 Double_t diff_step = miniBins.at(miniBins.size()-1)/number_of_points;
//	 cout << "diff_step: " << diff_step << endl;
//	 Double_t diff_CDF = 0.;
//	 Double_t current_point = 0.;
//	 for(Int_t i = 0; i<number_of_points; i++){
//		 Double_t diff_CDF = TMath::Abs(CDFs_variation->Eval(current_point)-CDFs_time->Eval(current_point));
//		 if(diff_CDF>max_diff)max_diff = diff_CDF;
//		 current_point+=diff_step;
//	 }
//
//	 cout << "Max diff: " << max_diff << endl;

//      //-----------------------------Create event list for Likelihood---------------------------
//
//	  Int_t selected = 0;
//	  FILE *fout;
//	  fout = fopen("Mrk421_2014flare_ONEvents_Kazuma_Leyre_noHECut.txt","wb");
//	  for(int i = 0; i<Npoints; i++)
//	  {
//		  if(Ep[i]>=cutInE_Low && Ep[i]<=cutInE_High){
//			  fprintf(fout,"%f    %f   \n", tp[i], Ep[i]);
//			  selected++;
//		  }
//	  }
//	  cout << "Selected events: " << selected << endl;
//
//

//      //-----------------------------Create root files---------------------------
//
//
	  	  TFile* ON_outputfile = new TFile("Histo_36ON_Paper.root","recreate");
//	  	  time_ONevents->SetName("ONevents");
//	  	  time_ONevents->Write();
//	  	  time_ONvariation->SetName("ONvariation");
	  	  time_ONvariation->Write("ONvariation");
	  	  energy->Write("energy_dist");
	  	  ON_outputfile->Close();


}
