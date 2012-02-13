#ifndef ReactorExperiment_h
#define ReactorExperiment_h 1
// ReactorExperiment.h
//
// A basic class for generating Monte Carlo data for a reactor
// coherent neutrino experiment.  This can be used as input to
// the classes for sensitivity analysis.
//
// Last modified: Dec. 14, 2011
// Created by: Adam Anderson
// adama@mit.edu

#include <TH1F.h>
#include <TRandom3.h>

class ReactorExperiment
{
public:
    ReactorExperiment(const char* thename);
    ~ReactorExperiment();
    
    void SetRun(Double_t detMass, Double_t megaWatts, Int_t neutrons, Int_t protons, Double_t expTime, Double_t sin2thetaS2, Double_t deltaMassSquare, Double_t dist, const char* reactorIsotope);
    //void SetTime(Double_t expTime);
    //void SetSinThetaS(Double_t angle);
    //void SetDeltaMassSquare(Double_t deltaMassSquare);
    void RunTruth(Int_t seed);
    void Run(Int_t seed);
    TH1F* GetEnergySpectrum();
    TH1F* GetDetectedNuEnergySpectrum();
    TF1* GetTrueNuEnergySpectrum();
    TF1* GetRecoilFunc();
    Int_t GetNEvt();
    Double_t GetTime();
    
private:
    Double_t Gfermi;
    Double_t Sin2ThetaW;
    Double_t InverseMeVtoCm;
    Double_t elementaryCharge;
    Double_t mass;
    Double_t time;
    Double_t sin2theta2;
    Double_t deltaMsq;
    Double_t distance;
    Double_t nucleonMass;
    Double_t Qweak;
    Double_t NuSpectrumMinE;
    Double_t NuSpectrumMaxE;
    Double_t activity;
    Int_t nNeutrons;
    Int_t nProtons;
    Int_t nEvt;
    
    const char* isotope;
    const char* name;
    
    TH1F* energySpectrum;
    TH1F* NeutrinoEvtHistogram;
    TF1* NeutrinoEnergySpectrum;
    TF1* DiffRecoilSpectrumAtConstE;
    TF1* DiffRecoilSpectrumAtConstT;
    TF1* RecoilSpectrumFunc;
    TF1* DetectedNuSpectrum;
    TF1* ScatteringCS;
    
    TRandom3* randGen;
    
    Int_t GenerateNumOfNuRecoils(Double_t activity, TF1* RecoilSpectrumArg, Double_t SpectrumMin, Double_t SpectrumMax, UInt_t seed);
    void FillRecoilAndNuSpectra(TH1F* recoilHisto, Int_t numOfEvents, TF1* RecoilSpectrumArg, TF1* DiffRecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, Int_t seed);
    void FillRecoilSpectrum(TH1F* recoilHisto, Int_t numOfEvents, TF1* RecoilSpectrumArg, TF1* DiffRecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, Int_t seed);
    Double_t RecoilSpectrum(Double_t* recoilEnergy, Double_t*);
    Double_t MaxRecoilEnergy(Double_t NeutrinoEnergy, Double_t mTarget);
};
#endif
