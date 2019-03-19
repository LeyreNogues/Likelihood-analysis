/*
 * Smoothness_test.C
 *
 *  Created on: Jul 24, 2018
 *      Author: lnogues
 */

#include <Rtypes.h>
#include <stddef.h>
#include <TArrayD.h>
#include <TAttLine.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <TRandom3.h>
#include<TSystem.h>
#include<TLegend.h>



using namespace std;

const Int_t Np = 100000;   // Number of events (hardcoded)
vector<double> E, t, E_sel, t_sel;
Int_t numberOfEvents = 1474; //Kazuma, 1 OFF
Double_t Se = 0;
Double_t binwidth=60.; // In seconds(60)
Int_t ON_events = 9;
Double_t average_width;


void Smoothness_test(char *fname){

	E.clear();
	t.clear();
	E_sel.clear();
	t_sel.clear();

	//Read original OFF data
	cout.precision(15);
//	cout << "Reading Event List in: " << fname << endl;

	ifstream in;
	in.open(fname);
	Double_t v1,v2;
	int Npoints = 0;
	while(1)
	{
		in >> v1 >> v2;
		if(v2<120.) continue;
		if (!in.good()) break;
		t.push_back(v1);
		E.push_back(v2);
		Npoints++;
	}

//	cout << "Npoints(Number of events): " << Npoints << endl;
//	cout << "Length of the array E: " << E.size() << endl;
//	cout << "Length of the array t: " << t.size() << endl;

	//Select sub-sample
//	cout << "Going to selection..." << endl;
		TRandom3 *ran = new TRandom3();
		ran->SetSeed(0);
		for(Int_t i=0; i<numberOfEvents; i++){
			Int_t number = E.size();
	//		cout << "Legth: " << number << endl;
			Int_t component = ran->Uniform(0,number);
	//		cout << "component: " << component << endl;
	//		cout << "E[component]: " << E[component] << endl;

			E_sel.push_back(E.at(component));
			t_sel.push_back(t.at(component));
//			E.erase(E.begin()+component);
//			t.erase(t.begin()+component);

		  }

//		cout << "Length of E: " << E.size() << endl;
//		cout << "Length of E_sel: " << E_sel.size() << endl;
//		cout << "Length of t: " << t.size() << endl;
//		cout << "Length of t_sel: " << t_sel.size() << endl;


		//Create the time histogram
		Double_t mintime = 1.77293506944444380e+03; //From wobble info
		Double_t maxtime = 1.77308576388889196e+03; //From wobble info
		const Int_t binsnum = 3000; //A big value
		Int_t numberOfBins; //The real one
		Double_t bin[binsnum+1] = {0};
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

//		cout << "hole[1]: " << hole[1] << endl;
//		cout << "hole[27]: " << hole[27] << endl;


		//Divide every interval and fill the bins.
		const Int_t numberOfIntervals = numberOfWobbles + (numberOfWobbles-1);
		Int_t I_bins[numberOfIntervals] = {0};
		Int_t currentBinIndex=0;
		Int_t lastBinIndex = 0;
//		cout << "numberOfIntervals: " << numberOfIntervals << endl;


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

//		 	cout << "lastBinIndex: " << lastBinIndex << endl;
		 	numberOfBins = lastBinIndex-1;
//		 	cout << "Number of bins in LightCurve: " << numberOfBins << endl;


		 	TH1D* time_insec = new TH1D("Lightcurve", "Lightcurve", numberOfBins, bin);
		 	time_insec->SetStats(0);

		 	Double_t initial_time = time_insec->GetBinLowEdge(1);
		 	Double_t final_time = time_insec->GetBinLowEdge(numberOfBins+1);
		 	Double_t total_time = final_time-initial_time;
//		 	cout << "Initial: " << initial_time << " final: " << final_time << endl;
//		 	cout << "Total time: " << total_time << endl;


			  for ( int a = 0; a < numberOfEvents; a++)
			  {
					  time_insec->Fill(t_sel[a]);
			  }

//			  TCanvas* selectedtimecanvas = new TCanvas();
//			  time_insec->GetXaxis()->SetTitle("Time(s)");
//			  time_insec->GetYaxis()->SetTitle("Selected events");
//			  time_insec->Draw("");

			  	  Double_t ave_width = 0.;
			  	  Int_t number_of_widths= 0;
			 	//Compute effective time.
			 	for(Int_t i=1; i<=numberOfBins; i++){
			 		if(time_insec->GetBinContent(i)==0){
//			 			cout << "hole" << endl;
			 			Double_t width = time_insec->GetBinWidth(i);
			 			total_time-= width;
			 		}
			 		else{
			 			ave_width+=time_insec->GetBinWidth(i);
			 			number_of_widths++;
			 		}
			 	}

//			 	cout << "Effective time: " << total_time << endl;
//			 	cout << "Average_width: " << ave_width/number_of_widths << endl;
			 	average_width = ave_width/number_of_widths;


			  //Order events in increasing time
//			  cout << "Now I order the events in increasing time..."  << endl;
			  Double_t temp_t;
			  Double_t temp_E;
			  for(Int_t i=0;i<numberOfEvents;i++)
			  {
		//		  cout << i << endl;
				  for(Int_t j=0;j<numberOfEvents-i;j++)
				  {
					  if(t_sel[j]>t_sel[j+1])
					  {
						  temp_t = t_sel[j];
						  t_sel[j] = t_sel[j+1];
						  t_sel[j+1]=temp_t;
						  temp_E = E_sel[j];
						  E_sel[j] = E_sel[j+1];
						  E_sel[j+1]=temp_E;
					  }
				  }
			  }

//			  cout << "First event: " << t_sel[0] << " last event: " << t_sel[numberOfEvents-1] << endl;

			  //-----------------------------Time bins for fixed - ON events---------------------------

//			  	  for(Int_t i = 0; i<numberOfIntervals; i++){
////			  		  cout << "I_bins: " << I_bins[i] << " ";
//			  	  }
//			  	  cout << endl;
////			  	  cout << "numberOfIntervals: " << numberOfIntervals << endl;
//
//			  	  Int_t counter = 0;
//			  	  Int_t event_index = 0;
//			  	  cout << "ON_events: " << ON_events<< endl;
//			  	  Int_t temp_ONevents = ON_events;
//			  	  std::vector<double> miniBins;
//			  	  Int_t bin_counter = 0;
//
//			  	  //First bin
//			  	  miniBins.push_back(hole[0]);
//
//			  	  //Loop over wobbles
//			  	  Int_t numberOfIntraBins = 0;
//			  	  Int_t rest= 0;
//			  	  for(Int_t a = 0; a<numberOfIntervals; a=a+2){//numberOfIntervals
//			  		  ON_events = temp_ONevents;
//			  		  Int_t eventsInWobble = 0;
////			  		  cout << "---------------Wobble-----------------: " << a << endl;
////			  		  cout << "Number of bins in the wobble: " << I_bins[a] << endl;
//			  		  for (Int_t i=bin_counter+1; i <=bin_counter+I_bins[a]; i++){
//			  //			  cout << "Content in bin: " << i << " is " << time_insec->GetBinContent(i) << endl;
//			  			  eventsInWobble+=time_insec->GetBinContent(i);
//
//			  		  }
////			  		  cout << "The number of events in the the wobble is: " << eventsInWobble << endl;
//			  		  numberOfIntraBins = TMath::Nint(eventsInWobble/ON_events);
//			  		  rest = eventsInWobble%ON_events;
//			  		  cout << "Cociente: " << numberOfIntraBins << endl;
//			  		  cout << "Resto: " << rest << endl;
//			  		  if(numberOfIntraBins<=rest){
//			  			  cout << "Resto>Cociente" << endl;
//			  			if(numberOfIntraBins == 0){
//							  ON_events = eventsInWobble;
//							  numberOfIntraBins = TMath::Nint(eventsInWobble/ON_events);
//							  rest = eventsInWobble%ON_events;
//						  }
//			  			else{
//			  			  ON_events = TMath::Nint(eventsInWobble/numberOfIntraBins);
//			  			  numberOfIntraBins = TMath::Nint(eventsInWobble/ON_events);
//			  			  rest = eventsInWobble%ON_events;
//			  			}
//
//			  			  cout << "New ON: " << ON_events << endl;
//			  			  cout << "Cociente: " << numberOfIntraBins << endl;
//			  			  cout << "Resto: " << rest << endl;
//
//			  		  }
//			  		  eventsInWobble+=time_insec->GetBinContent(bin_counter+I_bins[a]+1);
////			  		  cout << "Events adding one more (empty) bin: " << eventsInWobble << endl;
//			  //		  cout << "Wobble gap beggining: " << time_insec->GetBinLowEdge(bin_counter+I_bins[a]+1) << endl;
//			  //		  cout << "Wobble gap end: " << time_insec->GetBinLowEdge(bin_counter+I_bins[a]+2) << endl;
//			  		  bin_counter+=I_bins[a]+1;
//
//			  		  //Loop over events inside the wobble
//			  		  Int_t points_9 = 0;
//			  		  Int_t points_10 = 0;
//			  //		  cout << "event_index before loop: " << event_index << endl;
//			  		  for(Int_t j=event_index;j<event_index+eventsInWobble;j++){
//			  //				  cout << "counter: " << counter << endl;
//			  				  if(points_10<rest){
//			  					  if(counter==ON_events+1){
//			  //					  	cout << "Point for array of 10 events" << endl;
//			  //					  	cout << "event index: " << event_index << endl;
//			  //					  	cout << "t event: " << t[j] << endl;
//			  						  miniBins.push_back(t_sel[j]);
//			  						  points_10++;
//			  						  counter = 0;
//			  					  }
//			  				  }
//			  				  else{
//			  					  if(counter==ON_events){
//			  //					  	cout << "Point for array of 9 events" << endl;
//			  //					  	cout << "event index: " << event_index << endl;
//			  						  if(points_10+points_9 ==numberOfIntraBins-1){
//			  							  miniBins.push_back(hole[a+1]);
//			  							  counter=0;
//			  //							  cout << "event index before break: " << event_index << endl;
//
//			  							  break;
//			  						  }
//			  						  miniBins.push_back(t_sel[j]);
//			  //						  cout << "t event: " << t[j] << endl;
//			  						  points_9++;
//			  						  counter = 0;
//
//			  					  }
//			  				  }
//			  				  counter++;
//			  				  event_index++; //Inside the energy range
//			  //				  cout << "Counter: " << counter << endl;
//			  //				  cout << "event_index: " << event_index << endl;
//			  			  }
//
//			  		  	  if(a<numberOfIntervals-2)miniBins.push_back(hole[a+2]);
//			  		  }
//
//			  //	  for(Int_t i = 0; i<miniBins.size(); i++){
//			  //		  cout << miniBins.at(i) << " ";
//			  //	  }
//
//
////			  	  cout << "Last event index: " << event_index << endl;
////			  	  cout << "Size of miniBins: " << miniBins.size() << endl;
////			  	  cout << "Number of miniBins: " << miniBins.size()-1 << endl;
//			  	  Double_t* miniBins_pointer = miniBins.data();
////			  	  cout << "Last component: " << miniBins.at(miniBins.size()-1) << endl;
////			  	  cout << "Pre-Last component: " << miniBins.at(miniBins.size()-2) << endl;
//
//
//
//			  //Draw histogram to see if it is correct.
//			 TH1D* time_ONevents = new TH1D("time_ONevents", "time_ONevents", miniBins.size()-1, miniBins_pointer);
//			 time_ONevents->SetStats(0);
//			 TH1D* time_ONvariation = new TH1D("time_ONvariation", "time_ONvariation", miniBins.size()-1, miniBins_pointer);
//
//			 for(Int_t i=0; i<numberOfEvents;i++){
//				 time_ONevents->Fill(t_sel[i]);
//			 }
//			 for(Int_t i = 1; i<miniBins.size(); i++){
//				 time_ONvariation->SetBinContent(i, time_ONevents->GetBinContent(i)/time_ONevents->GetBinWidth(i));
//			 }

//			 TCanvas* patata = new TCanvas();
//			 patata->Divide(2,1);
//			 patata->cd(1);
//			 time_ONevents->Draw();
//			 patata->cd(2);
//			 time_ONvariation->Draw();


			TH1D* time_insec_rate = new TH1D("Lightcurve_rate", "Lightcurve_rate", numberOfBins, bin);
			for(Int_t i=1; i<=numberOfBins; i++){
		//		cout << "Index" << i << endl;
		//		cout << "BIn: " << bin[i-1] << endl;
		//
		//		cout << "Point: " << time_insec->GetBinContent(i)/time_insec->GetBinWidth(i) << endl;
				time_insec_rate->SetBinContent(i, time_insec->GetBinContent(i)/time_insec->GetBinWidth(i));
			}


		   //-------------------------------Check smoothness of the histogram---------------------------------

//			Int_t number_bins = 0;
//			Double_t difference =  0.;
//			for(Int_t i = 1; i<miniBins.size(); i++){
//				if(time_ONvariation->GetBinContent(i+1)!=0 && time_ONvariation->GetBinContent(i)!=0){
//					difference += TMath::Abs(time_ONvariation->GetBinContent(i+1) - time_ONvariation->GetBinContent(i));
//					number_bins++;
//				}
////				else cout << "Empty bin!" << endl;
//			}

			Int_t number_bins = 0;
			Double_t difference =  0.;
			for(Int_t i = 1; i<numberOfBins; i++){
				if(time_insec_rate->GetBinContent(i+1)!=0 && time_insec_rate->GetBinContent(i)!=0){
					difference += TMath::Abs(time_insec_rate->GetBinContent(i+1) - time_insec_rate->GetBinContent(i));
//					difference += time_insec_rate->GetBinContent(i+1) - time_insec_rate->GetBinContent(i);

					number_bins++;
				}
//				else cout << "Empty bin!" << endl;
			}

//			cout << "Diference: " << difference << endl;
//			cout << "Bins: " << number_bins << endl;
			Se = difference/number_bins;
//			cout << "Smooth estimator: " << Se << endl;

//			time_ONvariation->Delete();
//			time_ONevents->Delete();
			time_insec->Delete();
			time_insec_rate->Delete();


			E.clear();
			E.resize(0);
			t.clear();
			t.resize(0);
			E_sel.clear();
			E_sel.resize(0);
			t_sel.clear();
			t_sel.resize(0);
}

void Loop_samples(Int_t number_of_samples){

	gSystem->Load("libMathMore");
	//Arrays to save results
	const Int_t iterations = 10;
	Double_t Sigmas[iterations] = {0.};
	Double_t Se_mean[iterations] = {0.};
	Double_t Se_RMS[iterations] = {0.};
	Double_t sigma_RMS[iterations] = {0.};
	Double_t width[iterations] = {0.};
	Double_t ave_width[iterations] = {0.};



	//Flare properties
	Int_t numberONevents = 12520;
	Int_t numberOFFevents = 4422;
	Double_t effective_time = 12180.3;//12180.3;// 13020.;//12376.9999990007;//12180.3; //seconds
	Double_t ratio_ON = numberONevents/effective_time;
	Double_t ratio_3OFF =(numberOFFevents/effective_time);
	Double_t ratio_1OFF = ratio_3OFF/3.;

	cout << "ratio_ON: " << ratio_ON << endl;
	cout << "ratio_3OFF: " << ratio_3OFF << endl;
	cout << "ratio_1OFF: " << ratio_1OFF << endl;


	for(Int_t sigma=3; sigma<=12; sigma++){ //8

		Sigmas[sigma-3] = sigma;
		ON_events = sigma*sigma;
		binwidth=ON_events/ratio_ON;

//		cout << "ON_events: " << ON_events << endl;
//		cout << "binwidth: " << binwidth << endl;

		for(Int_t i=0; i<number_of_samples; i++){

		Smoothness_test("./FinalEventLists/All_Events/Mrk421_2014flare_OFFEvents_Kazuma_Leyre_noHECut.txt");
//		cout << "--------------------------------Se: " << Se << endl;
		Se_mean[sigma-3]+=Se;
		Se_RMS[sigma-3]+=Se*Se;

		}

		Se_mean[sigma-3] = Se_mean[sigma-3]/number_of_samples;
		Se_RMS[sigma-3] = 2*TMath::Sqrt((Se_RMS[sigma-3]-(Se_mean[sigma-3]*Se_mean[sigma-3]*number_of_samples))/(number_of_samples-1.));
		ave_width[sigma-3] = average_width;
		cout << "ON_events: " << ON_events << endl;
		cout << "binwidth: " << binwidth << endl;
		cout << "ave_binwidth: " << ave_width[sigma-3] << endl;


		cout << "-------------Result -> Mean_Se: " << Se_mean[sigma-3] << " RMS: " << Se_RMS[sigma-3] << "------------" << endl;
	}

//	TGraphErrors* sigmavsSe = new TGraphErrors(iterations, Sigmas, Se_mean, sigma_RMS, Se_RMS);
//	sigmavsSe->SetTitle("Se - fixed-time-rate histogram");
//	sigmavsSe->GetXaxis()->SetTitle("sigma");
//	sigmavsSe->GetYaxis()->SetTitle("Se");
//	sigmavsSe->SetMarkerStyle(20);
//	sigmavsSe->SetTitle("");
//
//	TCanvas* rosa = new TCanvas();
//	sigmavsSe->Draw("AP");

    TF1* bessel_0= new TF1("J_0", "ROOT::Math::cyl_bessel_i(0,2*x)*TMath::Exp(-2*x)", 0., 100.);
    TF1* bessel_1= new TF1("J_1", "ROOT::Math::cyl_bessel_i(1,2*x)*TMath::Exp(-2*x)", 0., 100.);
    TF1* sumbessel =  new TF1("SUM" , "ROOT::Math::cyl_bessel_i(0,2*x)*TMath::Exp(-2*x)+ROOT::Math::cyl_bessel_i(1,2*x)*TMath::Exp(-2*x)", 0. ,100.);


    TCanvas* canvas22 = new TCanvas();
	TGraph* graph = new TGraph(10);
	for(Int_t i = 0; i< iterations; i++){
		Double_t mu=Sigmas[i]*Sigmas[i]/2.83/3.;
		width[i]=Sigmas[i]*Sigmas[i]/ratio_ON;
		Double_t y = (2*mu/width[i])*sumbessel->Eval(mu);
		graph->SetPoint(i, Sigmas[i], y);
	}
	graph->SetLineColor(2);
	graph->Draw("AL");

	Int_t N=100;
	TGraph* graph_ON2 = new TGraph(N);
	for(Int_t i = 0; i< N; i++){
		Double_t mu=i*i;
		Double_t y = (2*ratio_ON)*TMath::Sqrt(1./TMath::Pi()/mu);
		graph_ON2->SetPoint(i, i, y);
	}
	graph_ON2->SetLineColor(9);
	graph_ON2->SetLineWidth(2);
	graph_ON2->Draw("L");

	TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
	legend->AddEntry(graph,"Background estimation","l");
	legend->AddEntry(graph_ON2,"Signal estimation","l");
	legend->Draw();


//	TGraph* graph2 = new TGraph(iterations);
//	for(Int_t i = 0; i< iterations; i++){
//		Double_t mu=Sigmas[i]*Sigmas[i];
//		width[i]=ave_width[i];
//		Double_t y = (2.*mu/width[i])*sumbessel->Eval(mu);
//		graph2->SetPoint(i, Sigmas[i], y);
//	}
//	graph2->SetLineColor(4);
//	graph2->Draw("L");



//    TCanvas* besselcanvas =  new TCanvas();
//    bessel_0->Draw();
//    bessel_1->SetLineColor(1);
//    bessel_1->Draw("same");
//    sumbessel->SetLineColor(4);
//    sumbessel->Draw("same");

}

