/*
 * CombLik.C
 *
 *  Created on: May 12, 2016
 *      Author: lnogues
 *
 *      Code to combine the Likelihood analysis of two sources: Mrk 501 and PG1553.
 *      Everything will be done first in a separate way and then combined at the end.
 *
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
#include <TStyle.h>
using namespace std;


//Int_t simnumMrk;
//Int_t simnumPG;
Int_t simnum;
const Int_t totalsim = 1;//Only for the combination
Double_t mintau[totalsim] = 0;
Double_t downtau[totalsim] = 0; //Combined results
Double_t uptau[totalsim] = 0;
const Int_t Nevents = 3000;
Double_t EMrk[Nevents], tMrk[Nevents];
Double_t EPG[Nevents], tPG[Nevents];


TF1* fitMrk;
TF1* spectrumMrk;
TF1* fitPG;
TF1* spectrumPG;

TGraph* profileMrk = new TGraph();
TGraph* normprofileMrk = new TGraph();
TGraph* NprofileMrk = new TGraph();

TGraph* profilePG = new TGraph();
TGraph* normprofilePG = new TGraph();
TGraph* NprofilePG = new TGraph();

TGraph* Combprofile = new TGraph();
TGraph* normCombprofile = new TGraph();


Double_t tinitMrk = 0;
Double_t tinitPG = 0;
Double_t tfinMrk = 1610;
Double_t tfinPG = 8000;

Double_t EmaxMrk = 8.;
Double_t EmaxPG = 1.;
Double_t EminMrk = 0.3;
Double_t EminPG = 0.4;

Double_t doblegaus(Double_t* x, Double_t* par)
{

	Double_t bkg, bkg2, tmax1, sigma,tmax2, sigma2;
	bkg = 18.32;
	tmax1 = 2540.;
	sigma = 798;
	bkg2 = 15.62;
	tmax2 = 5898.;
	sigma2 = 1070.;

	TF1 *gaus1 = new TF1("gaus1", gaus, tinitPG, tfinPG);
	gaus1->SetParameters(bkg, tmax1, sigma);
	TF1 *gaus2 = new TF1("gaus2", gaus, tinitPG, tfinPG);
	gaus2->SetParameters(bkg2, tmax2, sigma2);

	Double_t returnValue = gaus1->Eval(x[0])+gaus2->Eval(x[0]);
	delete gaus1;
	delete gaus2;
	return returnValue;
}

Double_t gaus(Double_t* x, Double_t* par) //with background
{
	Double_t funci = par[0]*exp(-(pow(x[0]-par[1],2)/(2*pow(par[2],2))));
	return funci;
}




void inputs()
{
	//Fixed templates for both sources


	//Mrk501 templates
	fitMrk = new TF1("fitMrk","gaus", tinitMrk, tfinMrk);
	fitMrk->SetParameters(1, 805., 219.);
	spectrumMrk = new TF1("spectrumMrk", "TMath::Power(x,-[0])", EminMrk, EmaxMrk);
	spectrumMrk->SetParameter(0, 2.2);



	//PG1553 templates
	fitPG = new TF1("fitPG",doblegaus, tinitPG, tfinPG, 0);
	spectrumPG = new TF1("spectrumPG", "TMath::Power(x,-[0])", EminPG, EmaxPG);
	spectrumPG->SetParameter(0, 4.8);


}

void MrkLik(char *Mrkfname)
{
	//Events reading

		const Int_t binsMrk = 200;
		Double_t LMrk[binsMrk];
		Double_t L2Mrk[binsMrk];
		Double_t tauMrk[binsMrk];
		Double_t N1Mrk[binsMrk], N2Mrk[binsMrk];
		Int_t NLCMrk, NSPMrk = 0;
		Double_t CHIMrk = 0;

		ifstream in;
		in.open(Mrkfname);
		Double_t v1,v2 = 0;
		int NpointsMrk = 0;
		Int_t lineMrk = 0;
		if(in.fail()) // checks to see if file opened
		    {
		      cout << "Reading file failed" << endl;
		    }
		else{
		while(1)
	    	{
			if(lineMrk == 0) //Events in LC
			{
				in >> v1;
				NLCMrk = v1;
				lineMrk++;
			}
			if(lineMrk == 1) //Events in Fit
			{
				in >> v1;
				NSPMrk = v1;
				lineMrk++;
			}
			if(lineMrk == 2) //Chi value in PL fit
			{
				in >> v1;
				CHIMrk = v1;
				lineMrk++;
			}
			if(lineMrk>2) //Events
			{
				in >> v1 >> v2;
				if (!in.good()) break;
				tMrk[NpointsMrk] = v1;
				EMrk[NpointsMrk] = v2;
				NpointsMrk++;
				lineMrk++;
			}
	    	}
		}

		//cout << Npoints << endl;
		//cout << NLC << " " << NSP << " " << CHI << endl;


		//-----------------------------------------Fit Ranges------------------------------------

		Double_t sumMrk = 0;
		for(Int_t a = 0; a < NpointsMrk; a++)
		{
			sumMrk+= EMrk[a];
		}

		Double_t meanEMrk = sumMrk/NpointsMrk;

		//cout << "Mean E: " << meanE <<endl;

		Double_t initauMrk = -200;

		//-----------------------------------------LOOP ON TAU------------------------------------

		for(Int_t j=0; j<=binsMrk; j++)
		{
			//if(j==0){
			tauMrk[j] = initauMrk+j*2;
			LMrk[j] = 0;

			N1Mrk[j] = 1;
			N2Mrk[j] = fitMrk->Integral(tinitMrk-tauMrk[j]*meanEMrk, tfinMrk-tauMrk[j]*meanEMrk);

			//cout << N1[j] << " " << N2[j] << " "  << endl;


			//-----------------------------------------LOOP ON EVENTS------------------------------------

			for(Int_t i=0; i<NpointsMrk; i++)
			{
				if (EMrk[i]>=EminMrk && EMrk[i]<=EmaxMrk && tMrk[i]>=tinitMrk && tMrk[i]<=tfinMrk)
				{
					Double_t PDFMrk[Nevents];

					PDFMrk[i] = (spectrumMrk->Eval(EMrk[i]))*(fitMrk->Eval(tMrk[i]-tauMrk[j]*EMrk[i]))/(NLCMrk*NSPMrk*N2Mrk[j]);
					/*cout << "energy "<< E[i] << endl;
					cout << "spectrum "<< spectrum->Eval(E[i]) << endl;
					cout << "collar "<< collar->Eval(E[i]) << endl;
					cout << "LC "<< fit->Eval(t[i]-tau[j]*E[i]) << endl;
					cout << "NNN "<< NLC*NSP*N1[j] << endl;
					cout << "PDF "<< PDF[i] << endl;*/

					if (PDFMrk[i] < 1.) //We add the weights here
					{
						Double_t lnPDFMrk[Nevents];
						lnPDFMrk[i]= TMath::Log(PDFMrk[i]);
						//cout << "lnPDF " << lnPDF[i] << endl;

						if( lnPDFMrk[i] < 0.)
						{
							LMrk[j] += -2*lnPDFMrk[i];
							//cout << "L " << L[j] << endl;
						}
						else
						{
							cout << "lnPDF is positive" << endl;
						}
					}
					else
					{
						cout << "Over 1 PDF " << i << endl;
					}

				}
				else
				{
					continue;
				}
			}
			//cout << "N " << N2[j] << endl;
			//cout << "L " << L[j] << endl;
			//cout << "tau " << tau[j] << endl;
			profileMrk->SetPoint(j, tauMrk[j], LMrk[j]);
			NprofileMrk->SetPoint(j, tauMrk[j], N2Mrk[j]);

			//} //test tau
		}

		//-------------------------------- Profile, tau min and upper limits -----------------------

		/*TCanvas* rofilecanvasMrk = new TCanvas();
		profileMrk->GetXaxis()->SetTitle("tau");
		profileMrk->GetYaxis()->SetTitle("-2lnL");
		profileMrk->SetMarkerColor(2);
		profileMrk->Draw("AP*");*/


		Double_t minimumMrk = LMrk[0];
		Double_t tauminMrk = tauMrk[0];
		Int_t minbin = 0;
		for(Int_t j=1; j<=binsMrk; j++)
		{
			if(LMrk[j] < minimumMrk)
				{
				minimumMrk = LMrk[j];
				tauminMrk = tauMrk[j];
				minbin = 0;
				}
		}

		cout << "Mrk01: Minimum: " << minimumMrk << " for tau " << tauminMrk << endl;
		//cout << "Bin of the min: " << minbin << endl;
		mintau[simnum]= tauminMrk;

		for(Int_t j=0; j<=binsMrk; j++)
		{
			L2Mrk[j] = LMrk[j]-minimumMrk;
			normprofileMrk->SetPoint(j, tauMrk[j], L2Mrk[j]);
		}


		/*TCanvas* normprofilecanvasMrk = new TCanvas();
		normprofileMrk->GetXaxis()->SetTitle("tau");
		normprofileMrk->GetYaxis()->SetTitle("-2#Delta lnL");
		normprofileMrk->SetMarkerColor(2);
		normprofileMrk->Draw("AP*");*/

		/*TCanvas* profilesMrk = new TCanvas();
		profilesMrk->SetTitle("Mrk501");
		profilesMrk->Divide(2,1);
		profilesMrk->cd(1);
		profileMrk->SetTitle("Mrk501");
		profileMrk->GetXaxis()->SetTitle("tau");
		profileMrk->GetYaxis()->SetTitle("-2lnL");
		profileMrk->SetMarkerColor(2);
		profileMrk->Draw("AP*");
		profilesMrk->cd(2);
		normprofileMrk->GetXaxis()->SetTitle("tau");
		normprofileMrk->GetYaxis()->SetTitle("-2#Delta lnL");
		normprofileMrk->SetMarkerColor(2);
		normprofileMrk->Draw("AP*");*/

		/*TCanvas* Ncanvas = new TCanvas();
		Ncanvas->Divide(2,1);
		Ncanvas->cd(1);
		normprofile->GetXaxis()->SetTitle("tau");
		normprofile->GetYaxis()->SetTitle("-2#Delta lnL");
		normprofile->Draw("AP*");
		Ncanvas->cd(2);
		Nprofile->GetXaxis()->SetTitle("tau");
		Nprofile->GetYaxis()->SetTitle("Normalization");
		Nprofile->Draw("AP*");*/

		Double_t ULMrk, LLMrk;
		ULMrk=LLMrk=0;
		Double_t eMrk = 0.6;
		Double_t value1Mrk, value2Mrk;

		for(Int_t j=minbin; j<=binsMrk; j++)
		{
			//cout << "LL j: " << j << " tau[j]: " << tau[j] << endl;
			value1Mrk =L2Mrk[j];
			//cout <<"value for Mrk501: " << value1Mrk << " for tau: " << tauMrk[j] << endl;
			if(value1Mrk >= 2.71-eMrk && value1Mrk <= 2.71+eMrk)
			{
				LLMrk=tauMrk[j];
				break;
			}
		}

		for(Int_t j=binsMrk; j>=minbin; j--)
		{
			//cout << "UL j: " << j << " tau[j]: " << tau[j] << endl;
			value2Mrk = L2Mrk[j];
			if(value2Mrk >= 2.71-eMrk && value2Mrk <= 2.71+eMrk)
			{
				ULMrk=tauMrk[j];
				break;
			}
		}

		cout << "Mrk50: LL: " << LLMrk << " UL: " << ULMrk << endl;
		//cout << "LLeval: " << normprofile->Eval(LL) << " ULeval: " << normprofile->Eval(UL) << endl;

		downtau[simnum] = LLMrk;
		uptau[simnum] = ULMrk;


}

void PGLik(char *PGfname)
{
	//Events reading

		const Int_t binsPG = 300;
		Double_t LPG[binsPG];
		Double_t L2PG[binsPG];
		Double_t tauPG[binsPG];
		Double_t N1PG[binsPG], N2PG[binsPG];
		Int_t NLCPG, NSPPG = 0;
		Double_t CHIPG = 0;

		ifstream in;
		in.open(PGfname);
		Double_t v1,v2 = 0;
		int NpointsPG = 0;
		Int_t linePG = 0;
		if(in.fail()) // checks to see if file opened
		    {
		      cout << "Reading file failed" << endl;
		    }
		else{
		while(1)
	    	{
			if(linePG == 0) //Events in LC
			{
				in >> v1;
				NLCPG = v1;
				linePG++;
			}
			if(linePG == 1) //Events in Fit
			{
				in >> v1;
				NSPPG = v1;
				linePG++;
			}
			if(linePG == 2) //Chi value in PL fit
			{
				in >> v1;
				CHIPG = v1;
				linePG++;
			}
			if(linePG>2) //Events
			{
				in >> v1 >> v2;
				if (!in.good()) break;
				tPG[NpointsPG] = v1;
				EPG[NpointsPG] = v2;
				NpointsPG++;
				linePG++;
			}
	    	}
		}

		//cout << Npoints << endl;
		//cout << NLC << " " << NSP << " " << CHI << endl;


		//-----------------------------------------Fit Ranges------------------------------------

		Double_t sumPG = 0;
		for(Int_t a = 0; a < NpointsPG; a++)
		{
			sumPG+= EPG[a];
		}

		Double_t meanEPG = sumPG/NpointsPG;

		//cout << "Mean E: " << meanE <<endl;

		Double_t initauPG = -1500;

		//-----------------------------------------LOOP ON TAU------------------------------------

		for(Int_t j=0; j<=binsPG; j++)
		{
			//if(j==0){
			tauPG[j] = initauPG+j*10;
			LPG[j] = 0;

			N1PG[j] = 1;
			N2PG[j] = fitPG->Integral(tinitPG-tauPG[j]*meanEPG, tfinPG-tauPG[j]*meanEPG);

			//cout << N1[j] << " " << N2[j] << " "  << endl;


			//-----------------------------------------LOOP ON EVENTS------------------------------------

			for(Int_t i=0; i<NpointsPG; i++)
			{
				if (EPG[i]>=EminPG && EPG[i]<=EmaxPG && tPG[i]>=tinitPG && tPG[i]<=tfinPG)
				{
					Double_t PDFPG[Nevents];

					PDFPG[i] = (spectrumPG->Eval(EPG[i]))*(fitPG->Eval(tPG[i]-tauPG[j]*EPG[i]))/(NLCPG*NSPPG*N2PG[j]);
					/*cout << "energy "<< E[i] << endl;
					cout << "spectrum "<< spectrum->Eval(E[i]) << endl;
					cout << "collar "<< collar->Eval(E[i]) << endl;
					cout << "LC "<< fit->Eval(t[i]-tau[j]*E[i]) << endl;
					cout << "NNN "<< NLC*NSP*N1[j] << endl;
					cout << "PDF "<< PDF[i] << endl;*/

					if (PDFPG[i] < 1.) //We add the weights here
					{
						Double_t lnPDFPG[Nevents];
						lnPDFPG[i]= TMath::Log(PDFPG[i]);
						//cout << "lnPDF " << lnPDF[i] << endl;

						if( lnPDFPG[i] < 0.)
						{
							LPG[j] += -2*lnPDFPG[i];
							//cout << "L " << L[j] << endl;
						}
						else
						{
							cout << "lnPDF is positive" << endl;
						}
					}
					else
					{
						cout << "Over 1 PDF " << i << endl;
					}

				}
				else
				{
					continue;
				}
			}
			//cout << "N " << N2[j] << endl;
			//cout << "L " << L[j] << endl;
			//cout << "tau " << tau[j] << endl;
			profilePG->SetPoint(j, tauPG[j], LPG[j]);
			NprofilePG->SetPoint(j, tauPG[j], N2PG[j]);

			//} //test tau
		}

		//-------------------------------- Profile, tau min and upper limits -----------------------

		/*TCanvas* profilecanvasPG = new TCanvas();
		profilePG->GetXaxis()->SetTitle("tau");
		profilePG->GetYaxis()->SetTitle("-2lnL");
		profilePG->SetMarkerColor(4);
		profilePG->Draw("AP*");*/

		Double_t minimumPG = LPG[0];
		Double_t tauminPG = tauPG[0];
		for(Int_t j=1; j<=binsPG; j++)
		{
			if(LPG[j] < minimumPG)
				{
				minimumPG = LPG[j];
				tauminPG = tauPG[j];
				}
		}

		cout << "PG: Minimum: " << minimumPG << " for tau " << tauminPG << endl;
		mintau[simnum]= tauminPG;

		for(Int_t j=0; j<=binsPG; j++)
		{
			L2PG[j] = LPG[j]-minimumPG;
			normprofilePG->SetPoint(j, tauPG[j], L2PG[j]);
		}


		/*TCanvas* normprofilecanvasPG = new TCanvas();
		normprofilePG->GetXaxis()->SetTitle("tau");
		normprofilePG->GetYaxis()->SetTitle("-2#Delta lnL");
		normprofilePG->SetMarkerColor(4);
		normprofilePG->Draw("AP*");*/

		/*TCanvas* profilesPG = new TCanvas();
		profilesPG->SetTitle("PG1553+113");
		profilesPG->Divide(2,1);
		profilesPG->cd(1);
		profilePG->SetTitle("PG1553+113");
		profilePG->GetXaxis()->SetTitle("tau");
		profilePG->GetYaxis()->SetTitle("-2lnL");
		profilePG->SetMarkerColor(4);
		profilePG->Draw("AP*");
		profilesPG->cd(2);
		normprofilePG->GetXaxis()->SetTitle("tau");
		normprofilePG->GetYaxis()->SetTitle("-2#Delta lnL");
		normprofilePG->SetMarkerColor(4);
		normprofilePG->Draw("AP*");*/

		/*TCanvas* Ncanvas = new TCanvas();
		Ncanvas->Divide(2,1);
		Ncanvas->cd(1);
		normprofile->GetXaxis()->SetTitle("tau");
		normprofile->GetYaxis()->SetTitle("-2#Delta lnL");
		normprofile->Draw("AP*");
		Ncanvas->cd(2);
		Nprofile->GetXaxis()->SetTitle("tau");
		Nprofile->GetYaxis()->SetTitle("Normalization");
		Nprofile->Draw("AP*");*/

		Double_t ULPG, LLPG;
		ULPG=LLPG=0;
		Double_t ePG = 0.8;
		Double_t value1PG, value2PG;

		for(Int_t j=0; j<=binsPG; j++)
		{
			//cout << "LL j: " << j << " tau[j]: " << tau[j] << endl;
			value1PG =L2PG[j];
			if(value1PG >= 2.71-ePG && value1PG <= 2.71+ePG)
			{
				LLPG=tauPG[j];
				break;
			}
		}

		for(Int_t j=binsPG; j>=0; j--)
		{
			//cout << "UL j: " << j << " tau[j]: " << tau[j] << endl;
			value2PG = L2PG[j];
			if(value2PG >= 2.71-ePG && value2PG <= 2.71+ePG)
			{
				ULPG=tauPG[j];
				break;
			}
		}

		cout << "PG1553: LL: " << LLPG << " UL: " << ULPG << endl;
		//cout << "LLeval: " << normprofile->Eval(LL) << " ULeval: " << normprofile->Eval(UL) << endl;

		downtau[simnum] = LLPG;
		uptau[simnum] = ULPG;


}

void CombCurve()
{
	//Add two graphs by evaluation them
	//First, check what source has a wider range
	Int_t min, max;
	min = max = 0;
	/*cout << "Mrk min: " << profileMrk->GetXaxis()->GetXmin() << endl;
	cout << "Mrk max: " << profileMrk->GetXaxis()->GetXmax() << endl;
	cout << "PG min: " << profilePG->GetXaxis()->GetXmin() << endl;
	cout << "PG max: " << profilePG->GetXaxis()->GetXmax() << endl;*/

	if(profileMrk->GetXaxis()->GetXmin() < profilePG->GetXaxis()->GetXmin()) min = profileMrk->GetXaxis()->GetXmin();
	else min = profilePG->GetXaxis()->GetXmin();
	if(profileMrk->GetXaxis()->GetXmax() > profilePG->GetXaxis()->GetXmax()) max = profileMrk->GetXaxis()->GetXmax();
	else max = profilePG->GetXaxis()->GetXmax();
	cout << "min: " << min << " max: " << max << endl;

	static const Int_t combpoint = 0;
	Double_t comminval= profileMrk->Eval(min)+profilePG->Eval(min);
	Int_t commintau= min;
	Int_t pointingraph = 0;
	Double_t L[10000]; //Too large
	Double_t tau[10000]; //Too large
	Double_t a, b;
	//Addition by evaluation and searching for the minimum
	for(Int_t i = min; i <= max; i=i+5){
		L[combpoint]= profileMrk->Eval(i)+profilePG->Eval(i);
		tau[combpoint] = i;
		Combprofile->SetPoint(combpoint, i, L[combpoint]);
		//cout << "Min value: " << comminval <<  " L value: "  << L[combpoint] << endl;
		if(L[combpoint] < comminval){
			comminval = L[combpoint];
			commintau = tau[combpoint];
			pointingraph = combpoint;
			Combprofile->GetPoint(combpoint, a, b);
			//cout << "Point x: "  << a << " y: " << b << endl;
		}
		combpoint++;
	}
	//cout << "Number of points in graph: " << combpoint-- << endl;

	//cout << "Point in graph of minimum: " << pointingraph << endl;

	cout << "Combined min tau: "  << commintau << endl;

	TCanvas* allprofiles = new TCanvas();
	allprofiles->SetTitle("PG1553+113");
	allprofiles->Divide(3,1);
	allprofiles->cd(1);
	profilePG->SetTitle("PG1553+113 Curve");
	profilePG->GetXaxis()->SetTitle("tau");
	profilePG->GetYaxis()->SetTitle("-2lnL");
	profilePG->SetMarkerColor(4);
	profilePG->Draw("AP*");
	allprofiles->cd(2);
	Combprofile->SetTitle("Comb Curve");
	Combprofile->GetXaxis()->SetTitle("tau");
	Combprofile->GetYaxis()->SetTitle("-2lnL");
	Combprofile->SetMarkerColor(1);
	Combprofile->Draw("AP*");
	allprofiles->cd(3);
	profileMrk->SetTitle("Mrk501 Curve");
	profileMrk->GetXaxis()->SetTitle("tau");
	profileMrk->GetYaxis()->SetTitle("-2lnL");
	profileMrk->SetMarkerColor(2);
	profileMrk->Draw("AP*");

	//Normalization and limits for the combination
	Double_t Lnorm[combpoint];
	for(Int_t i = 0; i< combpoint; i++){
		normCombprofile->SetPoint(i, tau[i], L[i]- comminval);
		Lnorm[i]= L[i]- comminval;
	}

	cout << "Number of points in second graph: " << normCombprofile->GetN() << endl;

	TCanvas* normprof = new TCanvas();
	normCombprofile->GetXaxis()->SetTitle("tau");
	normCombprofile->GetYaxis()->SetTitle("-2#Delta lnL");
	normCombprofile->SetMarkerColor(1);
	normCombprofile->Draw("AP*");

	Double_t UL, LL;
	UL=LL=0;
	Double_t e = 0.8;
	Double_t value1, value2;

	for(Int_t j=0; j<=pointingraph; j++)
	{
		value1 =Lnorm[j];
		if(value1 >= 2.71-e && value1 <= 2.71+e)
		{
			LL=tau[j];
			break;
		}
	}

	for(Int_t j=combpoint-1; j>=pointingraph; j--)
	{
		value2 = Lnorm[j];
		if(value2 >= 2.71-e && value2 <= 2.71+e)
		{
			UL=tau[j];
			break;
		}
	}

	cout << "Comb: LL: " << LL << " UL: " << UL << endl;

}

void WCombCurve()
{
	//Add two graphs by evaluation them
	//First, check what source has a wider range
	Double_t zMkr, zPG;
	zMkr = 0.034;
	zPG = 0.49;
	Double_t wPG = zPG/zMkr;
	cout << "Weight: " << wPG << endl;

	//We weight the curve of PG
	for (Int_t i = 0; i<profilePG->GetN(); i++)
	{
		profilePG->GetY()[i]*= wPG;
	}
	Int_t min, max;
	min = max = 0;
	/*cout << "Mrk min: " << profileMrk->GetXaxis()->GetXmin() << endl;
	cout << "Mrk max: " << profileMrk->GetXaxis()->GetXmax() << endl;
	cout << "PG min: " << profilePG->GetXaxis()->GetXmin() << endl;
	cout << "PG max: " << profilePG->GetXaxis()->GetXmax() << endl;*/

	if(profileMrk->GetXaxis()->GetXmin() < profilePG->GetXaxis()->GetXmin()) min = profileMrk->GetXaxis()->GetXmin();
	else min = profilePG->GetXaxis()->GetXmin();
	if(profileMrk->GetXaxis()->GetXmax() > profilePG->GetXaxis()->GetXmax()) max = profileMrk->GetXaxis()->GetXmax();
	else max = profilePG->GetXaxis()->GetXmax();
	cout << "min: " << min << " max: " << max << endl;

	static const Int_t combpoint = 0;
	Double_t comminval= profileMrk->Eval(min)+profilePG->Eval(min);
	Int_t commintau= min;
	Int_t pointingraph = 0;
	Double_t L[10000]; //Too large
	Double_t tau[10000]; //Too large
	Double_t a, b;
	//Addition by evaluation and searching for the minimum
	for(Int_t i = min; i <= max; i=i+2){
		L[combpoint]= profileMrk->Eval(i)+profilePG->Eval(i);
		tau[combpoint] = i;
		Combprofile->SetPoint(combpoint, i, L[combpoint]);
		//cout << "Min value: " << comminval <<  " L value: "  << L[combpoint] << endl;
		if(L[combpoint] < comminval){
			comminval = L[combpoint];
			commintau = tau[combpoint];
			pointingraph = combpoint;
			Combprofile->GetPoint(combpoint, a, b);
			//cout << "Point x: "  << a << " y: " << b << endl;
		}
		combpoint++;
	}
	//cout << "Number of points in graph: " << combpoint-- << endl;

	//cout << "Point in graph of minimum: " << pointingraph << endl;

	cout << "Combined min tau: "  << commintau << endl;

	TCanvas* allprofiles = new TCanvas();
	allprofiles->SetTitle("PG1553+113");
	allprofiles->Divide(3,1);
	allprofiles->cd(1);
	profilePG->SetTitle("PG1553+113 Curve");
	profilePG->GetXaxis()->SetTitle("tau");
	profilePG->GetYaxis()->SetTitle("-2lnL");
	profilePG->SetMarkerColor(4);
	profilePG->Draw("AP*");
	allprofiles->cd(2);
	Combprofile->SetTitle("Comb Curve");
	Combprofile->GetXaxis()->SetTitle("tau");
	Combprofile->GetYaxis()->SetTitle("-2lnL");
	Combprofile->SetMarkerColor(1);
	Combprofile->Draw("AP*");
	allprofiles->cd(3);
	profileMrk->SetTitle("Mrk501 Curve");
	profileMrk->GetXaxis()->SetTitle("tau");
	profileMrk->GetYaxis()->SetTitle("-2lnL");
	profileMrk->SetMarkerColor(2);
	profileMrk->Draw("AP*");

	//Normalization and limits for the combination
	Double_t Lnorm[combpoint];
	for(Int_t i = 0; i< combpoint; i++){
		normCombprofile->SetPoint(i, tau[i], L[i]- comminval);
		Lnorm[i]= L[i]- comminval;
	}

	cout << "Number of points in second graph: " << normCombprofile->GetN() << endl;

	TCanvas* normprof = new TCanvas();
	normCombprofile->GetXaxis()->SetTitle("tau");
	normCombprofile->GetYaxis()->SetTitle("-2#Delta lnL");
	normCombprofile->SetMarkerColor(1);
	normCombprofile->Draw("AP*");

	Double_t UL, LL;
	UL=LL=0;
	Double_t e = 0.8;
	Double_t value1, value2;

	for(Int_t j=0; j<=pointingraph; j++)
	{
		value1 =Lnorm[j];
		if(value1 >= 2.71-e && value1 <= 2.71+e)
		{
			LL=tau[j];
			break;
		}
	}

	for(Int_t j=combpoint-1; j>=pointingraph; j--)
	{
		value2 = Lnorm[j];
		if(value2 >= 2.71-e && value2 <= 2.71+e)
		{
			UL=tau[j];
			break;
		}
	}

	cout << "Comb: LL: " << LL << " UL: " << UL << endl;

}


void CombLik()
{
	inputs();

		//for(Int_t s=0; s<totalsim; s++)
		//{
	Int_t s =0;
	simnum = s;
	MrkLik(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/Mrk501/New/tau0/flareNBbase%i.txt", simnum));
	PGLik(Form("/home/lnogues/workspace_cpp/LIVanalysis/Simulation/PG1553/tau0/flareNB%i.txt", simnum));
	//CombCurve();
	WCombCurve();
	//cout << simnum << " tau: " << mintau[simnum] << " LL: " << downtau[simnum] << " UL: " << uptau[simnum] << endl;

		//}
}



