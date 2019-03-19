/*
 * readEventsfromFlute.C
 *
 *  Created on: May 10, 2017
 *      Author: Leyre Nogués
 *
 *  This code is to extract from Flute files the list of ON and OFF events
 *  with a cut in theta, the one indicated in the rc file.
 *  The cut in hadroness is already done by Flute.
 *  The code allows to do a cut in Energy, if desired.
 *  It is valid for BkgMode0 (Wobble partner) and BkgMode1 (Simultaneous).
 */

#include <iostream>
#include <fstream>

#include <TApplication.h>
#include <TArrayD.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCollection.h>
#include <TFile.h>
#include <TGClient.h>
#include <TGraph2D.h>
#include <TGraphAsymmErrors.h>
#include <TGTab.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TLatex.h>
#include <TNtupleD.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TVectorD.h>
#include "TLegend.h"

#include "MAnalysisProblems.h"
#include "MAnalysisWizard.h"
#include "MArgs.h"
#include "MAstro.h"
#include "MAverageEnergy.h"
#include "MBinning.h"
#include "MCalcExcess.h"
#include "MContinue.h"
#include "MFillH.h"
#include "MGeomCamMagic.h"
#include "MGraphicsWizard.h"
#include "MEnergyEst.h"
#include "MEnv.h"
#include "MEvtLoop.h"
#include "MF.h"
#include "MFHadAlpha.h"
#include "MH3.h"
#include "MHadAlphaCut.h"
#include "MHadronness.h"
#include "MHAlphaEnergyTheta.h"
#include "MHEffectiveOnTime.h"
#include "MHExcessEnergyTheta.h"
#include "MHFlux.h"
#include "MHillas.h"
#include "MHMcCollectionArea.h"
#include "MHMcEnergyMigration.h"
#include "MLog.h"
#include "MLogManip.h"
#include "MMcCollectionAreaCalc.h"
#include "MParList.h"
#include "MPointingPos.h"
#include "MRawRunHeader.h"
#include "MReadMarsFile.h"
#include "MReadTree.h"
#include "MSrcPosCalc.h"
#include "MStatusDisplay.h"
#include "MStereoPar.h"
#include "MTaskList.h"
#include "MTheta2vsEest.h"
#include "MStatusArray.h"


using namespace std;

void readEventsfromFlute(TString sOutputFluteName){

	cout.precision(12);

	TFile* fileOutputFlute = TFile::Open(sOutputFluteName, "READ");
		if (!fileOutputFlute) {
			cout << "ERROR: file " << sOutputFluteName << " not found..."	<< endl;
		}

	//-------------------------------Create text files for events----------------------------------


	cout << "Obtaining energy and time from events on ON an OFF regions" << endl;

	//Get the value of the alpha cut
	MHadAlphaCut* mHadAlphaCut =
			dynamic_cast<MHadAlphaCut*>(fileOutputFlute->Get("HadTheta2Cuts"));
	if (!mHadAlphaCut) {
		cout << "ERROR: MHadAlphaCut object not found... " << endl;
		cout << "HadTheta2Cuts not found" << endl;
	}


	MTheta2vsEest* mTheta2vsEst =
			dynamic_cast<MTheta2vsEest*>(fileOutputFlute->Get("MTheta2vsEest")); //This name depends on MARS version!
	if (!mTheta2vsEst) {
		cout << "ERROR: MTheta2vsEest object not found... " << endl;
		cout << "MTheta2vsEest not found" << endl;
	}


	Double_t tau_Bkg1 = 1.;

	// Read background mode
	//		0 = off from wobble partner
	// 		1 = simultaneous background
	Int_t bckgMode = mTheta2vsEst->GetBckgMode();
	if(bckgMode == 0) cout << "Bkg mode is OFF from wobble partner" << endl;
	else if(bckgMode == 1) cout  << "Bkg mode is simultaneous background" << endl;
	else if(bckgMode != 0 && bckgMode!=1) cout << "Bkg mode is unknown" << endl;


	//Flute create NumberfOfPointingxNumberfOfPointing NTuples with indices i,j.
	//The i,i NTuples correspond to ON regions and i,j for OFF regions.
	//The OFF region NTuples are only filled if BkgMode = 0.


	//Number of pointings computation and NTuple reading
	const Int_t numPointings = (Int_t) TMath::Sqrt((Double_t) mTheta2vsEst->GetHist().GetSize());
	cout << "Number of Pointings: " << numPointings << endl;
	TTree* ntupleTree[(Int_t)TMath::Power(numPointings*1.,2)];

	Int_t index=0;
	for(Int_t i=0;i<numPointings;i++)
	{
		for(Int_t j=0;j<numPointings;j++)
		{
			index = numPointings*i+j;
			ntupleTree[index] = (TTree*) mTheta2vsEst->GetNtuple(i, j);
		}
	}

	//Select maximum and minimum time of the events
	Double_t mintime_bkg0 = 10000000;
	Double_t maxtime_bkg0 = 0;

	//Array to save the events (One array per pointing)
	TArrayD EventEnergy_ON[numPointings]={0};
	TArrayD EventTime_ON[numPointings]={0};
	TArrayD EventEnergy_OFF[numPointings]={0};
	TArrayD EventTime_OFF[numPointings]={0};
	Int_t numEvt_ON[numPointings] = {0};
	Int_t numEvt_OFF[numPointings] = {0};

	Int_t selected_events[numPointings] = {0};

	if (bckgMode == 0){

		//Get and normalize theta plots(to compute tau)
		//This should be done independently for every wobble position
		Double_t tau[numPointings] = {0.};
		Double_t binWidth=0.005;
		Int_t numBinsTheta2 = (TMath::Sqrt(mTheta2vsEst->GetNormRangeTheta2Max())-TMath::Sqrt(mTheta2vsEst->GetNormRangeTheta2Min()))/binWidth;
		TH1D* theta2NormOn[numPointings] = {0};
		TH1D* theta2NormOff[numPointings] = {0};
		Double_t numEventsNormOn[numPointings]= {0};
		Double_t numEventsNormOff[numPointings]= {0};

		//Reading the Ntuplas (ON and OFF). Cut in hadronness is already done.
		Int_t ntuplaInteger=0;
		for(Int_t xpoint = 0; xpoint < numPointings; xpoint++){

			//Security checks (Check if time is ordered inside wobbles, careful! Wobbles can be not ordered in time themselves)

			theta2NormOn[xpoint] = new TH1D(Form("theta2NormOn%i",xpoint),"",numBinsTheta2,mTheta2vsEst->GetNormRangeTheta2Min(),mTheta2vsEst->GetNormRangeTheta2Max());
			theta2NormOff[xpoint] = new TH1D(Form("theta2NormOff%i",xpoint),"",numBinsTheta2,mTheta2vsEst->GetNormRangeTheta2Min(),mTheta2vsEst->GetNormRangeTheta2Max());

			for (Int_t ypoint = 0; ypoint < numPointings; ypoint++){

				Double_t currenttime = 0.;
				Int_t ntuplaInteger=(Int_t) numPointings*xpoint+ypoint;

				Double_t zenith, azimuth, eest, time, dirx, diry, srcx, srcy;

				ntupleTree[ntuplaInteger]->SetBranchAddress("zenith", &zenith);
				ntupleTree[ntuplaInteger]->SetBranchAddress("azimuth", &azimuth);
				ntupleTree[ntuplaInteger]->SetBranchAddress("eest", &eest);
				ntupleTree[ntuplaInteger]->SetBranchAddress("time", &time);
				ntupleTree[ntuplaInteger]->SetBranchAddress("dirx", &dirx);
				ntupleTree[ntuplaInteger]->SetBranchAddress("diry", &diry);
				ntupleTree[ntuplaInteger]->SetBranchAddress("srcx", &srcx);
				ntupleTree[ntuplaInteger]->SetBranchAddress("srcy", &srcy);

				Int_t nEventsNtuppla=ntupleTree[ntuplaInteger]->GetEntries();

				cout << "Reading Ntuppla (" << xpoint << "," << ypoint <<"). Total number of events = " << nEventsNtuppla << " (before theta cut);" << endl;


				//Loop inside the events in the NTupla and do cut in theta. Also compute tau.
				Int_t counter=0;
				for (Int_t i = 0; i < (Int_t)nEventsNtuppla; i++){
					ntupleTree[ntuplaInteger]->GetEntry(i);

//					cout << "Event: " << i << endl;

					//Time check
					if(currenttime < time) currenttime = time;
					else if(currenttime > time){
						cout << "WOBBLE EVENTS NOT ORDERED IN TIME" << endl;
						cout << "Im am in NTuple: " << xpoint << "," << ypoint << " at the event: " << i << endl;
						cout << "Last time was: " << currenttime << " and the time now is " << time << endl;
						return;
					}

					//Cuts in theta (thetaMin2Cut has to be < 0, except for Ring MC)
					Double_t thetaMin2Cut = mHadAlphaCut->GetAlphaMinCut(eest,zenith);
					Double_t theta2Cut = mHadAlphaCut->GetAlphaCut(eest, zenith);
//					cout << "thetaMinCut: " << thetaMin2Cut << " thetaCut: " << theta2Cut << endl;
					if(thetaMin2Cut > 0){
						cout << "PROBLEM! THETHA MIN CUT IS BIGGER THAN 0" << endl;
						return;
					}

					//Compute the theta value of the events (given the source position)
					//The source position is already adapted for ON and OFF regions in BkgMode0
					Double_t theta2 = mTheta2vsEst->GetTheta2(dirx, diry, srcx, srcy);
//					cout << "Theta value of the event: " << theta2 << endl;
					if(theta2<0.) continue; //Wrong event, do not take it into account

					//Do the cut in theta (ON and OFF regions)
					if (theta2 < theta2Cut){
						if(time < mintime_bkg0) mintime_bkg0 = time;
						if(time > maxtime_bkg0) maxtime_bkg0 = time;
						if(xpoint==ypoint)
						{
//							cout << "SELECTED as ON" << endl;
							EventEnergy_ON[xpoint].Set(numEvt_ON[xpoint] + 1);
							EventTime_ON[xpoint].Set(numEvt_ON[xpoint] + 1);
							EventEnergy_ON[xpoint].SetAt(eest, numEvt_ON[xpoint]);
							EventTime_ON[xpoint].SetAt(time, numEvt_ON[xpoint]);
							numEvt_ON[xpoint]++;
							selected_events[xpoint]++;
						}
						else
						{
//							cout << "SELECTED as OFF" << endl;
							EventEnergy_OFF[xpoint].Set(numEvt_OFF[xpoint] + 1);
							EventTime_OFF[xpoint].Set(numEvt_OFF[xpoint] + 1);
							EventEnergy_OFF[xpoint].SetAt(eest, numEvt_OFF[xpoint]);
							EventTime_OFF[xpoint].SetAt(time, numEvt_OFF[xpoint]);
							numEvt_OFF[xpoint]++;
							selected_events[xpoint]++;
						}
					}
					else if(theta2>mTheta2vsEst->GetNormRangeTheta2Min() && theta2<mTheta2vsEst->GetNormRangeTheta2Max())
					{// ON/OFF normalization
//						cout << "Normalization theta region is: " << mTheta2vsEst->GetNormRangeTheta2Min() << " - " << mTheta2vsEst->GetNormRangeTheta2Max() << endl;
						if(xpoint==ypoint)	theta2NormOn[xpoint]->Fill(theta2);
						else				theta2NormOff[xpoint]->Fill(theta2);
					}

					if(xpoint == 0 && ypoint == 0 && i == 0){ //Only one time
						theta2NormOn[xpoint]->Sumw2();
						theta2NormOff[xpoint]->Sumw2();
					}
					counter++;

				}//Events loop

//				cout << "Events in this Ntupla: " << counter << endl;
				numEventsNormOn[xpoint]=theta2NormOn[xpoint]->GetEntries();
				numEventsNormOff[xpoint]=theta2NormOff[xpoint]->GetEntries();
				tau[xpoint] = numEventsNormOff[xpoint]/numEventsNormOn[xpoint];

				cout << "tau for Ntuple " << xpoint << " " << ypoint << " is: "  << tau[xpoint] << endl;

			}//ypoint loop

			cout << "Selected events (in the wobble): " << selected_events[xpoint] << endl;
			cout << "As ON: " << numEvt_ON[xpoint] << endl;
			cout << "As OFF: " << numEvt_OFF[xpoint] << endl;

		}//xpoint loop

	}//Bkg 0 condition

	if(bckgMode == 1){

		cout << "I am in BkgMode 1" << endl;

		Int_t simultaneousBkgPositions = mTheta2vsEst->GetNumSimultaneousBgPositions();
		cout << "Number of Simultaneous positions: " << simultaneousBkgPositions << endl;

		//For simultaneous background mode, tau is common for every wobble.
		Double_t tau_value;
		tau_value=1./(double)(1.*simultaneousBkgPositions);
		tau_Bkg1 = tau_value;
		cout << "tau is: " << tau_Bkg1 << endl;

		//Reading the Ntuplas (only ON). Cut in hadronness is already done.
		Int_t ntuplaInteger=0;
		for(Int_t xpoint = 0; xpoint < numPointings; xpoint++){

			//Security checks (Check if time is ordered inside wobbles)
			Double_t currenttime = 0.;

			ntuplaInteger = (Int_t)numPointings*xpoint+xpoint;
//			cout << "ntuplaInteger: " << ntuplaInteger << endl;

			Double_t zenith, azimuth, eest, time, dirx, diry, srcx, srcy;

			ntupleTree[ntuplaInteger]->SetBranchAddress("zenith", &zenith);
			ntupleTree[ntuplaInteger]->SetBranchAddress("azimuth", &azimuth);
			ntupleTree[ntuplaInteger]->SetBranchAddress("eest", &eest);
			ntupleTree[ntuplaInteger]->SetBranchAddress("time", &time);
			ntupleTree[ntuplaInteger]->SetBranchAddress("dirx", &dirx);
			ntupleTree[ntuplaInteger]->SetBranchAddress("diry", &diry);
			ntupleTree[ntuplaInteger]->SetBranchAddress("srcx", &srcx);
			ntupleTree[ntuplaInteger]->SetBranchAddress("srcy", &srcy);

			Int_t nEventsNtuppla=ntupleTree[ntuplaInteger]->GetEntries();

			cout << "Reading Ntuppla (" << xpoint << "," << xpoint <<"). Total number of events = " << nEventsNtuppla << " (before theta cut);" << endl;

			//Loop inside the events in the NTupla and do cut in theta. Also compute tau.
			Int_t counter=0;
			for (Int_t i = 0; i < (Int_t)nEventsNtuppla; i++){
				ntupleTree[ntuplaInteger]->GetEntry(i);
//				cout << "Event " << i << endl;

				//Time check
				if(currenttime < time) currenttime = time;
				else if(currenttime > time){
					cout << "WOBBLE EVENTS NOT ORDERED IN TIME" << endl;
					cout << "Im am in NTuple: " << xpoint << " at the event: " << i << endl;
					cout << "Last time was: " << currenttime << " and the time now is " << time << endl;
					return;
				}

				//Cuts in theta (thetaMin2Cut has to be < 0, except for Ring MC)
//				cout << "Energy of the event is: " << eest << endl;
//				cout << "Zenith of the event is: " << zenith << endl;
//				cout << "Azimuth of the event is: " << azimuth << endl;
//				cout << "Time of the event is: " << time << endl;
//				cout << "Dir x of the event is: " << dirx << endl;
//				cout << "Dir y of the event is: " << diry << endl;
//				cout << "Source x of the event is: " << srcx << endl;
//				cout << "Source y of the event is: " << srcy << endl;

				Float_t thetaMin2Cut = mHadAlphaCut->GetAlphaMinCut(eest,zenith);
				Float_t theta2Cut = mHadAlphaCut->GetAlphaCut(eest, zenith);

//				cout << "ThetaMinCut is " << thetaMin2Cut << " ThetaCut is " << theta2Cut << endl;
				if(thetaMin2Cut > 0){
					cout << "PROBLEM! THETHA MIN CUT IS BIGGER THAN 0" << endl;
					return;
				}
				if(theta2Cut < 0){
					cout << "PROBLEM! THETHA CUT IS NEGATIVE" << endl;
					continue;
				}

				//Cut in theta (For ON region)
				Double_t theta2 = mTheta2vsEst->GetTheta2(dirx, diry, srcx, srcy);
//				cout << "Theta for ON: " << theta2 << endl;
				if(theta2<0.) continue; //Wrong event, do not take it into account
				if (theta2 < theta2Cut){
//						cout << "SELECTED as ON" << endl;
//						cout << "Energy: " << eest << " Time: " << time << endl;
						EventEnergy_ON[xpoint].Set(numEvt_ON[xpoint] + 1);
						EventTime_ON[xpoint].Set(numEvt_ON[xpoint] + 1);
						EventEnergy_ON[xpoint].SetAt(eest, numEvt_ON[xpoint]);
						EventTime_ON[xpoint].SetAt(time, numEvt_ON[xpoint]);
//						cout << "Energy array has in position " << numEvt_ON[xpoint] << " the value " << EventEnergy_ON[xpoint].GetAt(numEvt_ON[xpoint]) << endl;
//						cout << "Time array has in position " << numEvt_ON[xpoint] << " the value " << EventTime_ON[xpoint].GetAt(numEvt_ON[xpoint]) << endl;
						numEvt_ON[xpoint]++;
						selected_events[xpoint]++;
					}

				//Cut in theta (For OFF regions)
				Float_t psi = TMath::ATan2(srcy, srcx);
				Float_t radius = sqrt(srcx*srcx+srcy*srcy);
				//Loop for the OFF regions
				for (Int_t ioff = 0; ioff < simultaneousBkgPositions; ioff++){
					Float_t psioff = psi + (ioff+1) * TMath::TwoPi()/((Float_t)simultaneousBkgPositions+1.);
					Float_t offsrcx = radius*TMath::Cos(psioff);
					Float_t offsrcy = radius*TMath::Sin(psioff);
					Float_t offtheta2 = mTheta2vsEst->GetTheta2(dirx, diry, offsrcx, offsrcy);
//					cout << "Theta for OFF: " << offtheta2 << endl;

					if (offtheta2 < theta2Cut){
//						cout << "SELECTED as OFF" << endl;
//						cout << "Energy: " << eest << " Time: " << time << endl;
						EventEnergy_OFF[xpoint].Set(numEvt_OFF[xpoint] + 1);
						EventTime_OFF[xpoint].Set(numEvt_OFF[xpoint] + 1);
						EventEnergy_OFF[xpoint].SetAt(eest, numEvt_OFF[xpoint]);
						EventTime_OFF[xpoint].SetAt(time, numEvt_OFF[xpoint]);
//						cout << "Energy array has in position " << numEvt_OFF[xpoint] << " the value " << EventEnergy_OFF[xpoint].GetAt(numEvt_OFF[xpoint]) << endl;
//						cout << "Time array has in position " << numEvt_OFF[xpoint] << " the value " << EventTime_OFF[xpoint].GetAt(numEvt_OFF[xpoint]) << endl;
						numEvt_OFF[xpoint]++;
						selected_events[xpoint]++;
					}
				}

			}//Events loop

			cout << "Selected events: " << selected_events[xpoint] << endl;
			cout << "As ON: " << numEvt_ON[xpoint] << endl;
			cout << "As OFF: " << numEvt_OFF[xpoint] << endl;

		}//xpoint loop
	} //Bkg1 condition

	else if (bckgMode!= 0 && bckgMode!= 1){
		cout << "Code not developed for this background Mode" << endl;
		return;
	}


	cout << "We display the selected events" << endl;

	Int_t numEvtON = 0, numEvtOFF = 0;
	for(Int_t i=0; i< numPointings; i++){
//		cout << "numEvt_ON[i]: " << numEvt_ON[i] << " numEvt_OFF[i]: " << numEvt_OFF[i] << endl;
		numEvtON+=numEvt_ON[i];
		numEvtOFF+=numEvt_OFF[i];
	}
	cout << "Number of ON events: " << numEvtON << endl;
	cout << "Number of OFF events: " << numEvtOFF << endl;
	Double_t numEvtOFF_tau = numEvtOFF*tau_Bkg1;
	cout << "Number of Normalized OFF events: " << numEvtOFF_tau << endl;



	// Display results
	TCanvas* canveventse = new TCanvas();
	gStyle->SetOptStat(1111);
	canveventse->Modified();
	canveventse->SetLogy();

	//Histograms for energy
	TH1D *h_on_energy  = new TH1D("h_on_energy","",100,TMath::Log10(1e1),TMath::Log10(1e5));
	TH1D *h_off_energy = new TH1D("h_off_energy","",100,TMath::Log10(1e1),TMath::Log10(1e5));


	//Fill the histograms with the ON and OFF events of all pointings.
	Int_t maxON = 0;
	Int_t maxOFF = 0;
	for(Int_t i=0; i< numPointings; i++){
		maxON = EventEnergy_ON[i].GetSize();
		maxOFF = EventEnergy_OFF[i].GetSize();

//		cout << "maxON: "  << maxON << endl;
//		cout << "maxOFF: "  << maxOFF << endl;

		for (int j=0; j<maxON; j++){
			h_on_energy->Fill(TMath::Log10(EventEnergy_ON[i].GetAt(j)));
		}
		for (int k=0; k<maxOFF; k++){
			h_off_energy->Fill(TMath::Log10(EventEnergy_OFF[i].GetAt(k)), tau_Bkg1);
		}
	}

	cout << "Entries in ON energy histo: " << h_on_energy->GetEntries() << endl;
	cout << "Entries in OFF energy histo: " << h_off_energy->GetEntries() << endl;

	h_off_energy->SetStats(0);
	h_off_energy->SetLineWidth(4);
	h_off_energy->SetLineColor(kBlack);
	h_off_energy->GetXaxis()->SetTitle("Log_{10}(E[GeV])");
	h_off_energy->GetYaxis()->SetTitle("number of events");
	h_off_energy->Draw("");
	h_on_energy->SetLineColor(kRed);
	h_on_energy->Draw("same");
	gPad->Modified();
	gPad->Update();


	TLegend *legendevents=new TLegend(0.6,0.8,0.8,0.9,"","NDC");
	legendevents->SetFillColor(0);
	legendevents->SetLineColor(0);
	legendevents->SetBorderSize(0);
	legendevents->SetTextFont(42);
	legendevents->SetTextSize(0.025);
	legendevents->AddEntry(h_on_energy,"On events","l");
	legendevents->AddEntry(h_off_energy,"Off events","l");
	legendevents->Draw("same");

//	gPad->SaveAs("Flute_energy.pdf");
//	gPad->SaveAs("Flute_energy_Kazuma.png");

	// Display results of time:
	TCanvas* canveventst = new TCanvas();
	gStyle->SetOptStat(1111);
	Double_t cutInE = 0; //(GeV)

	//Select limits on histogram
	Double_t min_ON=EventTime_ON[0].GetAt(0), max_ON=0;
	Double_t min_OFF=EventTime_OFF[0].GetAt(0), max_OFF=0;

	if(bckgMode == 0){ //Events are not saved in order in the wobble array.
		min_ON = mintime_bkg0;
		min_OFF = mintime_bkg0;
		max_ON = maxtime_bkg0;
		max_OFF = maxtime_bkg0;
	}

	if(bckgMode == 1){ //Events are saved in order in the wobble array
		min_ON=EventTime_ON[0].GetAt(0), max_ON=0;
		min_OFF=EventTime_OFF[0].GetAt(0), max_OFF=0;
		for(Int_t i=0; i< numPointings; i++){
			if(min_ON > EventTime_ON[i].GetAt(0)) min_ON = EventTime_ON[i].GetAt(0);
			if(min_OFF > EventTime_OFF[i].GetAt(0)) min_OFF = EventTime_OFF[i].GetAt(0);
			if(max_ON < EventTime_ON[i].GetAt(EventEnergy_ON[i].GetSize()-1)) max_ON = EventTime_ON[i].GetAt(EventEnergy_ON[i].GetSize()-1);
			if(max_OFF < EventTime_OFF[i].GetAt(EventEnergy_OFF[i].GetSize()-1)) max_OFF = EventTime_OFF[i].GetAt(EventEnergy_OFF[i].GetSize()-1);
		}
	}



	cout << "Maximum ON time event: " << max_ON << endl;
	cout << "Minimum ON time event: " << min_ON << endl;
	cout << "Maximum OFF time event: " << max_OFF << endl;
	cout << "Minimum OFF time event: " << min_OFF << endl;

	TH1D *h_on_time  = new TH1D("h_on_time","",80, min_ON, max_ON);
	TH1D *h_off_time  = new TH1D("h_off_time","",80, min_OFF, max_OFF);

	for(Int_t i=0; i< numPointings; i++){
		Int_t maxON = EventEnergy_ON[i].GetSize();
		Int_t maxOFF = EventEnergy_OFF[i].GetSize();
		for (int j=0; j<maxON; j++){
		  if(TMath::Abs(EventTime_ON[i].GetAt(j))>cutInE){
			  h_on_time->Fill(EventTime_ON[i].GetAt(j));
			}
		}
		for (int k=0; k<maxOFF; k++){
		  if(TMath::Abs(EventTime_OFF[i].GetAt(k))>cutInE){
			  h_off_time->Fill(EventTime_OFF[i].GetAt(k), tau_Bkg1);
			}
		}
	}

	cout << "Entries in ON time histo: " << h_on_time->GetEntries() << endl;
	cout << "Entries in OFF time histo: " << h_off_time->GetEntries() << endl;

	h_on_time->SetStats(0);
	h_on_time->SetLineWidth(1);
	h_on_time->SetLineColor(kBlack);
	h_on_time->GetXaxis()->SetTitle("Time(MJD)");
	h_on_time->GetYaxis()->SetTitle("number of events");
	h_on_time->Draw("");
	h_off_time->SetLineWidth(1);
	h_off_time->SetLineColor(14);
	h_off_time->Draw("same");
	gPad->Modified();
	gPad->Update();


	TLegend *legendeventstime=new TLegend(0.6,0.8,0.8,0.9,"","NDC");
	legendeventstime->SetFillColor(0);
	legendeventstime->SetLineColor(0);
	legendeventstime->SetBorderSize(0);
	legendeventstime->SetTextFont(42);
	legendeventstime->SetTextSize(0.025);
	legendeventstime->AddEntry(h_on_time,"On events","l");
	legendeventstime->AddEntry(h_off_time,"Off events","l");
	legendeventstime->Draw("same");

//	gPad->SaveAs("Flute_time.pdf");
//	gPad->SaveAs("Flute_time_Kazuma.png");

//	Create text files with the ON and OFF events
	FILE *fout_ON;
	fout_ON = fopen("GRB_Koji_ONevents.txt","wb");
	FILE *fout_OFF;
	fout_OFF = fopen("GRB_Koji_OFFevents.txt","wb");

	for(Int_t i=0; i< numPointings; i++){
		maxON = EventEnergy_ON[i].GetSize();
		maxOFF = EventEnergy_OFF[i].GetSize();

		for (int j=0; j<maxON; j++){
			fprintf(fout_ON,"%0.8lf    %0.8lf   \n", EventTime_ON[i].GetAt(j), EventEnergy_ON[i].GetAt(j));
		}
		for (int k=0; k<maxOFF; k++){
			fprintf(fout_OFF,"%0.8lf    %0.8lf   \n", EventTime_OFF[i].GetAt(k), EventEnergy_OFF[i].GetAt(k));
		}
	}

}




