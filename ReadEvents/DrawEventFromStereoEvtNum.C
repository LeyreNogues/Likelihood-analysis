/*
 * DrawEventFromStereoEvtNum.C
 *
 *  Created on: Oct 9, 2017
 *      Author: lnogues
 */


#include<iostream>
#include<fstream>
#include<istream>
#include <TChain.h>
#include <MTime.h>
#include <TArrayF.h>
#include <MCerPhotEvt.h>
#include <MArrivalTime.h>
#include <MRawEvtHeader.h>
#include <TLatex.h>
#include <MGeomCamMagicTwo.h>
#include <MHCamera.h>


using namespace std;

////////////////////////////////////////////////////////////////
//
//         Draw MAGIC camera from calibrated files for a given
//                      stereo event number
//
////////////////////////////////////////////////////////////////


void draw_geometry(MGeomCam geom, TArrayD pixels) {

    MHCamera *cam = new MHCamera(geom);
    cam->SetBit(TH1::kNoStats|MHCamera::kCanDelete);
    cam->AddCamContent(pixels);
    cam->SetMinimum(0);
    cam->SetMaximum(50);
    cam->Draw("pixelindex");

}

void DrawEventFromStereoEvtNum(UInt_t stereoEvtNum){

    TChain rh("RunHeaders");
    rh.Add("./*.root");

    MGeomCamMagicTwo mGeomCam;

    rh.SetBranchStatus("*", 1);
    rh.SetBranchStatus("MGeomCam.*",             1);
    rh.SetBranchAddress("MGeomCam.MGeomCam", &mGeomCam);

    rh.GetEvent(0);

    MTime* mTime = 0;
    MCerPhotEvt* mCerPhotEvt = 0;
    MArrivalTime* mArrivaltime = 0;
    MRawEvtHeader* mRawEvtHeader = 0;
//    MCerPhotPix *cpix = 0;

    TChain ch("Events");
    ch.Add("./*.root");
    cout << "Number of events: " << ch.Sizeof() << endl;
    ch.SetBranchStatus("*", 1);

    ch.SetBranchStatus("MTime.*",             1);
    ch.SetBranchStatus("MCerPhotEvt.*",    1);
    ch.SetBranchStatus("MArrivalTime.*",    1);
    ch.SetBranchStatus("MRawEvtHeader.*",    1);

    ch.SetBranchAddress("MTime.",         &mTime);
    ch.SetBranchAddress("MCerPhotEvt.", &mCerPhotEvt);
    ch.SetBranchAddress("MArrivalTime.", &mArrivaltime);
    ch.SetBranchAddress("MRawEvtHeader.", &mRawEvtHeader);

    TArrayD pixels, times;
    TArrayF times_raw;
    bool eventSelected=0;
    double prevMJD;

    Double_t cPixSignal=0;

    cout << "I will get inside the loop: Number of Entries " << ch.GetEntries() << endl;

    for (Int_t ievt = 0; ievt < ch.GetEntries(); ievt++) {
        pixels.Reset();
        times.Reset();
        times_raw.Reset();
        ch.GetEvent(ievt);
//        cout << "mRawEvtHeader->GetStereoEvtNumber() = " << mRawEvtHeader->GetStereoEvtNumber() << endl;
        if (mRawEvtHeader->GetStereoEvtNumber() != stereoEvtNum)    continue;

        pixels.Set(mGeomCam.GetNumPixels());
        times.Set(mArrivaltime->GetSize());
        times_raw.Set(mArrivaltime->GetSize());

        for (Int_t ipix = 0; ipix < mCerPhotEvt->GetNumPixels(); ipix++) {
            MCerPhotPix &cpix = (*mCerPhotEvt)[ipix];
//            cpix = (MCerPhotPix*) mCerPhotEvt[ipix];
            if (cpix.GetPixId() == 0) cPixSignal = cpix.GetNumPhotons();
            times_raw = mArrivaltime->GetData();
            Double_t signal = cpix.GetNumPhotons();
            if (signal < 0.)
                continue;
            pixels[cpix.GetPixId()] = signal;
//            pixels[cpix->GetPixId()] = log10(signal);
            times[cpix.GetPixId()] = times_raw[cpix.GetPixId()];
        }

        cout.precision(17);
        draw_geometry(mGeomCam, pixels);

        cout << "Plotting image for M1 at MTime = " << mTime->GetMjd() << endl;
        TLatex *texf = new TLatex(0.04,0.79,"M1");
        texf->SetTextColor(1);
        texf->SetTextAlign(22);
        texf->SetTextSize(0.07);
        texf->Draw();
//        TLatex *texfCpixText = new TLatex(0.04,0.73,Form("cpix",cPixSignal));
//        texfCpixText->SetTextColor(1);
//        texfCpixText->SetTextAlign(22);
//        texfCpixText->SetTextSize(0.03);
//        texfCpixText->Draw();
        TLatex *texfCpix = new TLatex(0.04,0.69,Form("%2.1f phe",cPixSignal));
        texfCpix->SetTextColor(1);
        texfCpix->SetTextAlign(22);
        texfCpix->SetTextSize(0.03);
        texfCpix->Draw();


//        delete mTime;
//        delete mCerPhotEvt;
//        delete mArrivaltime;
//        delete mRawEvtHeader;
        return;
    }

//    delete mGeomCam;
    delete mTime;
    delete mRawEvtHeader;
    delete mCerPhotEvt;
    delete mArrivaltime;
    return;
}

