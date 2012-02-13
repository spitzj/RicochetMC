// ReactorExperiment.cxx
//
// A basic class for generating Monte Carlo data for a reactor
// coherent neutrino experiment.  This can be used as input to
// the classes for sensitivity analysis.
//
// Last modified: Dec. 30, 2011
// Created by: Adam Anderson
// adama@mit.edu

#include "ReactorExperiment.h"
#include <iostream>
#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TCanvas.h>

using namespace std;

ReactorExperiment::ReactorExperiment(const char* thename)
{
    name = thename;
    Gfermi = 1.16637e-11;				// [MeV^-2]
	Sin2ThetaW = 0.2387;
	InverseMeVtoCm = 1.97e-11;
    elementaryCharge = 1.602176565e-19;
    NuSpectrumMinE = 1.0;  // [MeV]
	NuSpectrumMaxE = 10.0; // [MeV]
    
    DiffRecoilSpectrumAtConstE = new TF1("DiffRecoilSpectrumAtConstE", "(x<[4])*TMath::Power([0]*[1],2)/(4*TMath::Pi()) * [2] * (1 - ([2] * x)/(2 * TMath::Power([3],2))) + (x>[4])*0", 0.0, 0.01);
    DiffRecoilSpectrumAtConstT = new TF1("DiffRecoilSpectrumAtConstT", "([3] < x / (1 + [2]/(2*x)))*TMath::Power([0]*[1],2)/(4*TMath::Pi()) * [2] * (1 - ([2] * [3])/(2 * TMath::Power(x,2))) + ([3] > x / (1 + [2]/(2*x)))*0", NuSpectrumMinE, NuSpectrumMaxE);
    RecoilSpectrumFunc = new TF1("RecoilSpectrum", this, &ReactorExperiment::RecoilSpectrum, 0.0, 0.01, 0, "ReactorExperiment", "RecoilSpectrum");
    ScatteringCS = new TF1("ScatteringCS", "TMath::Power([0]*[1]*[3],2)/(4*TMath::Pi()) * [2] * (x/(1 + [2]/(2*x))) * (1 - ([2]/(4*x*x)) * (x/(1 + [2]/(2*x))))", NuSpectrumMinE, NuSpectrumMaxE);
    randGen = new TRandom3();
    NeutrinoEnergySpectrum = NULL;
    DetectedNuSpectrum = NULL;
    energySpectrum = NULL;
    NeutrinoEvtHistogram = NULL;
}


ReactorExperiment::~ReactorExperiment()
{
    delete energySpectrum;
    delete NeutrinoEvtHistogram;
    delete name;
    delete isotope;
    delete DiffRecoilSpectrumAtConstE;
    delete DiffRecoilSpectrumAtConstT;
    delete RecoilSpectrumFunc;
    delete ScatteringCS;
    delete randGen;
}


void ReactorExperiment::SetRun(Double_t detMass, Double_t megaWatts, Int_t neutrons, Int_t protons, Double_t expTime, Double_t sin2thetaS2, Double_t deltaMassSquare, Double_t dist, const char* reactorIsotope)
{
    mass = detMass;
    activity = megaWatts/200.0/elementaryCharge;			  // [sec^-1]
    time = expTime;
    sin2theta2 = sin2thetaS2;
    nNeutrons = neutrons;
    nProtons = protons;
    deltaMsq = deltaMassSquare;
    distance = dist;
    isotope = reactorIsotope;
    nucleonMass = (nNeutrons * 939.565) + (nProtons * 938.272);     // [MeV]
    Qweak = nNeutrons - (1 - 4*Sin2ThetaW) * nProtons;
    
    if (NeutrinoEnergySpectrum != NULL)
        delete NeutrinoEnergySpectrum;
    
    if (strcmp(isotope, "U235") == 0)
    {
        // high-energy component taken from table VI of 1101.2663
        // low-energy component is a rough estimate from PRD 39, 3378
        NeutrinoEnergySpectrum = new TF1("NeutrinoEnergySpectrum", "(x > 1.0)*TMath::Exp([0] + [1]*x + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4) + [5]*TMath::Power(x,5)) * (1 - [6] * pow(TMath::Sin(0.0127 * [7] * [8] / x),2)) + (x < 1.0) * 2.5 * (1 - [6] * pow(TMath::Sin(0.0127 * [7] * [8] / x),2))", NuSpectrumMinE, NuSpectrumMaxE);
        NeutrinoEnergySpectrum->SetParameters(3.217, -3.111, 1.395, -0.369, 0.04445, -0.002053, sin2thetaS2, dist, deltaMassSquare);
    }
    else if (strcmp(isotope, "Pu239") == 0)
    {
        NeutrinoEnergySpectrum = new TF1("NeutrinoEnergySpectrum", "(x > 2.0)*TMath::Exp([0] + [1]*x + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4) + [5]*TMath::Power(x,5)) + (x < 2.0 && x > 1.0)*(2.3 - (x - 1)*(2.3 - TMath::Exp([0] + [1]*2.0 + [2]*4.0 + [3]*8.0 + [4]*16.0 + [5]*32.0))) * (1 - [6] * pow(TMath::Sin(0.0127 * [7] * [8] / x),2)) + (x < 1.0) * 2.3 * (1 - [6] * pow(TMath::Sin(0.0127 * [7] * [8] / x),2))", NuSpectrumMinE, NuSpectrumMaxE);
        NeutrinoEnergySpectrum->SetParameters(6.413, -7.432, 3.535, -0.882, 0.1025, -0.00455, sin2thetaS2, dist, deltaMassSquare);
    }
    else
    {
        cout << "Error: You did not enter a valid isotope name" << endl;
        return;
    }
    
	DiffRecoilSpectrumAtConstE->SetParameters(Gfermi, Qweak, nucleonMass);
    DiffRecoilSpectrumAtConstT->SetParameters(Gfermi, Qweak, nucleonMass);
    ScatteringCS->SetParameters(Gfermi, Qweak, nucleonMass, InverseMeVtoCm);
    
    if (energySpectrum != NULL)
        delete energySpectrum;
    if (NeutrinoEvtHistogram != NULL)
        delete NeutrinoEvtHistogram;
    if (DetectedNuSpectrum != NULL)
        delete DetectedNuSpectrum;
}


/*void ReactorExperiment::SetTime(Double_t expTime)
{
    time = expTime;
}


void ReactorExperiment::SetSinThetaS(Double_t angle)
{
    sin2theta2 = angle;
    NeutrinoEnergySpectrum->SetParameter(6,angle);
}

void ReactorExperiment::SetDeltaMassSquare(Double_t deltaMassSquare)
{
    deltaMsq = deltaMassSquare;
    NeutrinoEnergySpectrum->SetParameter(8,deltaMassSquare);
}*/


void ReactorExperiment::RunTruth(Int_t seed)
{    
    DetectedNuSpectrum = new TF1("DetectedNuSpectrum","ScatteringCS * NeutrinoEnergySpectrum", 0.0, 10.0);
    energySpectrum = new TH1F(name,"recoil energy spectrum;Energy [eVnr];Events / 10 eV",50,1000000*MaxRecoilEnergy(1.0,nucleonMass),500.0);
    NeutrinoEvtHistogram = new TH1F("NeutrinoEvtHistogram","detected #bar{#nu}_{e} energy spectrum;Energy [MeV];Events / 0.5 MeV",50,0.0,10);
    randGen->SetSeed(seed);
    
    nEvt = GenerateNumOfNuRecoils(activity, DetectedNuSpectrum, NuSpectrumMinE, NuSpectrumMaxE, seed);
    FillRecoilAndNuSpectra(energySpectrum, nEvt, DetectedNuSpectrum, DiffRecoilSpectrumAtConstE, NuSpectrumMinE, NuSpectrumMaxE, seed);
}


void ReactorExperiment::Run(Int_t seed)
{    
    DetectedNuSpectrum = new TF1("DetectedNuSpectrum","ScatteringCS * NeutrinoEnergySpectrum", 0.0, 10.0);
    energySpectrum = new TH1F(name,"recoil energy spectrum;Energy [eVnr];Events / 10 eV",50,1000000*MaxRecoilEnergy(1.0,nucleonMass),500.0);
    NeutrinoEvtHistogram = new TH1F("NeutrinoEvtHistogram","detected #bar{#nu}_{e} energy spectrum;Energy [MeV];Events / 0.5 MeV",50,0.0,10);
    randGen->SetSeed(seed);
    
    nEvt = GenerateNumOfNuRecoils(activity, DetectedNuSpectrum, NuSpectrumMinE, NuSpectrumMaxE, seed);
    FillRecoilSpectrum(energySpectrum, nEvt, DetectedNuSpectrum, DiffRecoilSpectrumAtConstE, NuSpectrumMinE, NuSpectrumMaxE, seed);
}


TH1F* ReactorExperiment::GetEnergySpectrum()
{
    return energySpectrum;
}


TH1F* ReactorExperiment::GetDetectedNuEnergySpectrum()
{
    return NeutrinoEvtHistogram;
}


TF1* ReactorExperiment::GetTrueNuEnergySpectrum()
{
    return NeutrinoEnergySpectrum;
}


Int_t ReactorExperiment::GetNEvt()
{
    return nEvt;
}


Double_t ReactorExperiment::GetTime()
{
    return time;
}


TF1* ReactorExperiment::GetRecoilFunc()
{
    return RecoilSpectrumFunc;
}


// Generate the number of events for this run.
// This is just Poisson statistics and the cross section.
Int_t ReactorExperiment::GenerateNumOfNuRecoils(Double_t activity, TF1* RecoilSpectrumArg, Double_t SpectrumMin, Double_t SpectrumMax, UInt_t seed)
{
	// Set up the constants
	Double_t AvogadroConst = 6.022e23;

	Double_t molarMass = nNeutrons + nProtons;
	
	// Calculate mean number of events
	Double_t SpectrumWeightedCrossSection = RecoilSpectrumArg->Integral(SpectrumMin, SpectrumMax);
	Double_t MeanEvents = time * AvogadroConst * activity * mass * SpectrumWeightedCrossSection / (4 * TMath::Pi() * molarMass * pow(distance,2));
    
	// Generate a random number of events from Poisson distribution
	Int_t numOfEvents = randGen->Poisson(MeanEvents);
	return numOfEvents;
}


// Function to generate the spectrum from coherent
// neutrino scattering.
void ReactorExperiment::FillRecoilAndNuSpectra(TH1F* recoilHisto, Int_t numOfEvents, TF1* RecoilSpectrumArg, TF1* DiffRecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, Int_t seed)
{	
    Double_t SurvivalProbability;
    Int_t eventCount = 0;
    
	// Now let's generate actually generate the events and
	// fill a histogram
	for (Int_t evtNum = 0; evtNum < numOfEvents; evtNum++)
    {
        // Get the neutrino energy
        Double_t thisNuEnergy = RecoilSpectrumArg->GetRandom(SpectrumMin, SpectrumMax);
           
        NeutrinoEvtHistogram->Fill(thisNuEnergy);
        
        // Get the recoil energy
        DiffRecoilSpectrum->SetParameter(3, thisNuEnergy);
        DiffRecoilSpectrum->SetParameter(4, MaxRecoilEnergy(thisNuEnergy, nucleonMass));
        DiffRecoilSpectrum->SetRange(0.0,MaxRecoilEnergy(thisNuEnergy, nucleonMass));
        Double_t thisRecoilEnergy = DiffRecoilSpectrum->GetRandom(0.0,MaxRecoilEnergy(thisNuEnergy, nucleonMass));
        recoilHisto->Fill(1000000*thisRecoilEnergy);
        eventCount++;
	}
}


// Function to generate the spectrum from coherent
// neutrino scattering, just generating a histogram
// and not generating individual events
void ReactorExperiment::FillRecoilSpectrum(TH1F* recoilHisto, Int_t numOfEvents, TF1* RecoilSpectrumArg, TF1* DiffRecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, Int_t seed)
{
    Double_t normalization = RecoilSpectrumFunc->Integral(RecoilSpectrumFunc->GetXmin(), RecoilSpectrumFunc->GetXmax());
    Int_t nBins = recoilHisto->GetNbinsX();
    Int_t eventCount = 0;
    Int_t nEvtInBin;
    Double_t binLowEdge, binWidth, avgEventsInBin;
    
    for (Int_t thisBin = 1; thisBin <= nBins; thisBin++) {
        binLowEdge = recoilHisto->GetBinLowEdge(thisBin) / 1000000.0;
        binWidth = recoilHisto->GetBinWidth(thisBin) / 1000000.0;
        
        avgEventsInBin = RecoilSpectrumFunc->Integral(binLowEdge, binLowEdge+binWidth) * numOfEvents / normalization;
        
        nEvtInBin = randGen->Poisson(avgEventsInBin);
        recoilHisto->SetBinContent(thisBin, nEvtInBin);
        eventCount += nEvtInBin;
    }
}


Double_t ReactorExperiment::RecoilSpectrum(Double_t* recoilEnergy, Double_t*)
{
    TF1 theRecoilSpectrum("theRecoilSpectrum", "DiffRecoilSpectrumAtConstT*NeutrinoEnergySpectrum", NuSpectrumMinE, NuSpectrumMaxE);
    theRecoilSpectrum.SetParameter(3,*recoilEnergy);
    Double_t value = theRecoilSpectrum.Integral(1.0,10.0);

    return value;
}


// Just computes the maximum recoil due to a neutrino
// of the given energy.  Note that all energies are in MeV.
Double_t ReactorExperiment::MaxRecoilEnergy(Double_t NeutrinoEnergy, Double_t mTarget)
{
	return NeutrinoEnergy / (1.0 + mTarget / (2.0*NeutrinoEnergy));
}



