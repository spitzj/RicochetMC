#include "OscillationFitter.h"
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <stdlib.h>
#include <TMath.h>
#include <iostream>
#include <TCanvas.h>

using namespace std;

OscillationFitter::OscillationFitter()
{ }

OscillationFitter::~OscillationFitter()
{
    delete theExperiment;
    delete GlobalScanContour;
}

void OscillationFitter::SetFitter(ReactorExperiment* theExpt, TH1F* eSpectrum)
{
    theExperiment = theExpt;
    energySpectrum = eSpectrum;
}

void OscillationFitter::GlobalScanFit(const Int_t nBinsDM2, Double_t minDM2, Double_t maxDM2, const Int_t nBinsSin22ThetaS, Double_t minSin22ThetaS, Double_t maxSin22ThetaS, Int_t seed)
{
    // Do high statistics MC estimate of the unoscillated signal
    // Note: this loses the original parameters of the experiment!!!
    // need to fix this bug in ReactorExperiment.cxx!!!
    // need to fix hardcoding!!!
    // Change this to be experiment passed to function--more streamlined
    Double_t highStatTime = 10000.0 * 365.0 * 24.0 * 3600.0;
    
    // Points to test in parameter space
    Double_t thetaPoints[nBinsSin22ThetaS+1];
    for (Int_t thisThetaBin = 0; thisThetaBin <= nBinsSin22ThetaS; thisThetaBin++)
        thetaPoints[thisThetaBin] = TMath::Power(10, thisThetaBin * (TMath::Log10(maxSin22ThetaS) - TMath::Log10(minSin22ThetaS)) / nBinsSin22ThetaS + TMath::Log10(minSin22ThetaS));
    
    Double_t massPoints[nBinsDM2+1];
    for (Int_t thisMassBin = 0; thisMassBin <= nBinsDM2; thisMassBin++)
        massPoints[thisMassBin] = TMath::Power(10, thisMassBin * (TMath::Log10(maxDM2) - TMath::Log10(minDM2)) / nBinsDM2 + TMath::Log10(minDM2));
    
    // Basic definitions
    GlobalScanContour = new TH2F("GlobalScanContour", "90% CL upper limit ('global scan')", nBinsSin22ThetaS, thetaPoints, nBinsDM2, massPoints);
    const Double_t chiSq90CL = 4.71;
    const Double_t pvalCL = 0.9;
    
    
    // Compute chi^2 at each point in parameter space
    // NOTE: Should implement some kind of more effcient
    // algorithm for scan to avoid computing in regions
    // of low significance.
    Double_t progress = 10;
    for (Int_t thisThetaBin = 1; thisThetaBin <= nBinsSin22ThetaS; thisThetaBin++)
    {
        for (Int_t thisMassBin = 1; thisMassBin <= nBinsDM2; thisMassBin++)
        {
            cout << "Scanning delta m^2 bin #" << thisMassBin << " and theta bin #" << thisThetaBin << endl;
            theExperiment->SetRun(5000.0, 5.5, 42, 32, highStatTime, thetaPoints[thisThetaBin], massPoints[thisMassBin], 400.0, "U235");
            theExperiment->Run(seed + thisThetaBin*nBinsDM2 + thisMassBin);
            TH1F* highStatESpectrum = theExperiment->GetEnergySpectrum();
            highStatESpectrum->Scale(0.0001);
            highStatESpectrum->Sumw2();   // set the errors
            
            Double_t pvalue = energySpectrum->Chi2Test(highStatESpectrum, "UW");
            GlobalScanContour->SetBinContent(thisThetaBin, thisMassBin, pvalue);
            
            /*TCanvas can1("can1");
            highStatESpectrum->SetLineColor(2);
            highStatESpectrum->Draw();
            energySpectrum->Draw("same");
            char filename[20];
            Int_t n = sprintf(filename, "testplot_%d_%d.pdf", thisThetaBin, thisMassBin);
            can1.SaveAs(filename);*/
        }
    }
}


TH2F* OscillationFitter::GetGlobalScanContour()
{
    return GlobalScanContour;
}

/*void OscillationFitter::RasterScanFit()
{
//stuff
}*/

/*void OscillationFitter::FeldmanCousinsFit()
{
//stuff
}*/