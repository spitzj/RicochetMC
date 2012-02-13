// ReactorNuAnalysis.C
//
// This is an analysis loop which takes the results of the 
// Geant4 background Monte Carlo and the Mathematica
// calculation of the neutrino signal and generates a
// random experiment to analyze.  At the moment, the
// analysis is very rudimentary and just fits a histogram
// to extract a signal.  In this future, a Bayesian
// analysis could be added.
// 
// Note: This is the ROOT/C++ version of this script.
// Another PyROOT/Python version will eventually be
// available and will supercede this version.  I maintain
// this version only because PyROOT can be tricky to set
// up, and users might find it useful to have a version
// that is guaranteed to be compatible with the
// precompiled binaries of ROOT.
//
// NOTE: Livetime effects still need to be added;
// exposure weighting of histograms should be more
// transparent when doing background subtraction
//
// Adam Anderson
// adama@mit.edu

#include <stdio.h>
#include <TH1.h>

void ReactorNuAnalysis()
{
	// -------- VARIABLE DEFINITIONS and SETUP --------
	// Style settings
	gStyle->SetOptStat("");
	
	// Constants
	Double_t NuSpectrumMinE = 0.0;  // [MeV]
	Double_t NuSpectrumMaxE = 10.0; // [MeV]
	Int_t seed = 43534;
	const Double_t Gfermi = 1.16637e-11;																			 // [MeV^-2]
	const Double_t Sin2ThetaW = 0.2387;
	const Double_t InverseMeVtoCm = 1.97e-11;
	const Double_t elementaryCharge = 1.602176565e-19;
	
	// Spectrum files (these are just samples for now)
	const char ReactorNeutronBackgroundFile[] = "SampleData.txt";
	
	// Define constants relevant for this run
	Double_t OnOffTimeRatio = 1.0/1.0;
	Double_t time = 100 * 24.0*3600.0;																				 // [sec]
	Double_t detMass = 5000.0;																								 // [g]
	Double_t distance = 200.0;																								 // [cm]
	Double_t activity = 6.0/200.0/elementaryCharge;														 // [sec^-1]
	const Int_t nNeutrons = 14;
	const Int_t nProtons = 14;
	const Double_t Qweak = nNeutrons - (1 - 4*Sin2ThetaW) * nProtons;
	const Double_t NucleonMass = (nNeutrons * 939.565) + (nProtons * 938.272);			 // [MeV]
	
	// Define the histograms
	TH1F* RecoilEvtHistogram = new TH1F("RecoilEvtHistogram","^{28}Si recoil energy spectrum;Energy [eV];Events / 10 eV",50,0.0,500.0);
	TH1F* NeutrinoEvtHistogram = new TH1F("NeutrinoEvtHistogram","#bar{#nu}_{e} energy spectrum;Energy [MeV];Events / 0.5 MeV",50,0.0,20);
	TH1F* TheoryRecoilEvtHistogram = new TH1F("TheoryRecoilEvtHistogram","High-statistics (\"theoretical\") Coherent Recoil Spectrum;Energy [eV];Events / 10 eV",50,0.0,500.0);
	TH1F* TheoryNeutrinoEvtHistogram = new TH1F("TheoryNeutrinoEvtHistogram","High-statistics (\"theoretical\") Neutrino Energy Spectrum;Energy [MeV] / 0.5 MeV;Events",50,0.0,20);
	TH1F* TheoryRecoilEvtHistogramFit = new TH1F("TheoryRecoilEvtHistogram","MC Fit;Events",50,0.0,0.002);
	TH1F* EMNoLukeBackground = new TH1F("EMNoLukeBackground","Background from EM recoils (before Luke Effect amplification;Energy [MeV];Events)",50,0.0,0.002);
	TH1F* ReactorNeutronBackground = new TH1F("ReactorNeutronBackground","Background from reactor neutrons;Energy [MeV];Events",50,0.0,0.002);
	TH1F* ReactorOnCosmoNeutronBackground = new TH1F("ReactorOnCosmoNeutronBackground","Reactor-on background from muon-induced neutrons;Energy [MeV];Events)",50,0.0,0.002);
	TH1F* ReactorOffCosmoNeutronBackground = new TH1F("ReactorOffCosmoNeutronBackground","Reactor-off background from muon-induced neutrons;Energy [MeV];Events)",50,0.0,0.002);	
	TH1F* ReactorOnHisto = new TH1F("ReactorOnHisto","Recoil Spectrum for Reactor-On Data;Energy [MeV];Events",50,0.0,0.002);
	TH1F* ReactorOffHisto = new TH1F("ReactorOffHisto","Recoil Spectrum for Reactor-Off Data;Energy [MeV];Events",50,0.0,0.002);
	TH1F* BackgroundSubtractedSignal = new TH1F("BackgroundSubtractedSignal","Recoil Spectrum for Reactor-Off Data;Energy [MeV];Events",50,0.0,0.002);
	
	// Energy spectra
	// (Spectral parameterizations from arXiv:1101.2663v3)
	TF1* NeutrinoEnergySpectrum = new TF1("NeutrinoSpectrum", "TMath::Exp([0] + [1]*x + [2]*TMath::Power(x,2) + [3]*TMath::Power(x,3) + [4]*TMath::Power(x,4) + [5]*TMath::Power(x,5))", NuSpectrumMinE, NuSpectrumMaxE);
	NeutrinoEnergySpectrum->SetParameters(3.217, -3.111, 1.395, -0.369, 0.04445, -0.002053);
	TF1* IntegratedRecoilSpectrum = new TF1("IntegratedRecoilSpectrum", "TMath::Power([0]*[1]*[3],2)/(4*TMath::Pi()) * [2] * (x/(1 + [2]/(2*x))) * (1 - ([2]/(4*x*x)) * (x/(1 + [2]/(2*x))))", NuSpectrumMinE, NuSpectrumMaxE);
	IntegratedRecoilSpectrum->SetParameters(Gfermi, Qweak, NucleonMass, InverseMeVtoCm);
	TF1* RecoilSpectrum = new TF1("RecoilSpectrum","IntegratedRecoilSpectrum * NeutrinoSpectrum", 0.0, 10.0);
	TF1* DiffRecoilSpectrumAtConstE = new TF1("DiffRecoilSpectrumAtConstE", "(x<[4])*TMath::Power([0]*[1],2)/(4*TMath::Pi()) * [2] * (1 - ([2] * x)/(2 * TMath::Power([3],2))) + (x>[4])*0", NuSpectrumMinE, 0.01);
	DiffRecoilSpectrumAtConstE->SetParameters(Gfermi, Qweak, NucleonMass);
	
	
	
	
	// -------- HISTOGRAM FILLING --------
	// Fill "experimental" histograms
	Int_t nEvt = GenerateNumOfNuRecoils(time, detMass, distance, activity, nNeutrons, nProtons, RecoilSpectrum, NuSpectrumMinE, NuSpectrumMaxE, seed);
	FillNuRecoilSpectrum(RecoilEvtHistogram, NeutrinoEvtHistogram, nEvt, nNeutrons, nProtons, RecoilSpectrum, DiffRecoilSpectrumAtConstE, NuSpectrumMinE, NuSpectrumMaxE, seed+4);
	FillRecoilSpectrumFromFile(ReactorNeutronBackground, 100.0, 10.0, ReactorNeutronBackgroundFile, seed+1);      // Testing purposes only.  CHANGE ME!!
	FillRecoilSpectrumFromFile(ReactorOnCosmoNeutronBackground, 50.0, 10.0, ReactorNeutronBackgroundFile, seed+2);		// Testing purposes only.  CHANGE ME!!
	FillRecoilSpectrumFromFile(ReactorOffCosmoNeutronBackground, 50.0, 10.0, ReactorNeutronBackgroundFile, seed+3);		// Testing purposes only.  CHANGE ME!!
	
	cout << "nEvt: " << nEvt << endl;
	
	// Fill high-statistics "theoretical" histograms
	FillNuRecoilSpectrum(TheoryRecoilEvtHistogram, TheoryNeutrinoEvtHistogram, 10000, nNeutrons, nProtons, RecoilSpectrum, DiffRecoilSpectrumAtConstE, NuSpectrumMinE, NuSpectrumMaxE, seed+5);
	TheoryNeutrinoEvtHistogram->Scale(nEvt/TheoryNeutrinoEvtHistogram->GetEntries());
	TheoryRecoilEvtHistogram->Scale(nEvt/TheoryNeutrinoEvtHistogram->GetEntries());
	
	// Combine the histograms into total
	ReactorOnHisto->Add(RecoilEvtHistogram);
	ReactorOnHisto->Add(ReactorNeutronBackground);
	ReactorOnHisto->Add(ReactorOnCosmoNeutronBackground);
	
	ReactorOffHisto->Add(ReactorOffCosmoNeutronBackground);
	
	
	
	
	// -------- HYPOTHESIS TESTING --------
	cout << "p-value between Reactor-On and Reactor-Off Data: " << ReactorOnHisto->Chi2Test(ReactorOffHisto) << endl;
	cout << "p-value between the simulated data and the Monte Carlo histogram: " << RecoilEvtHistogram->Chi2Test(TheoryRecoilEvtHistogram) << endl;
	
	
	
	
	// -------- BACKGROUND SUBTRACTION and FITTING --------
	// Normalize reactor-off data to reactor-on data by exposure
	BackgroundSubtractedSignal->Add(ReactorOnHisto, ReactorOffHisto, 1.0, -1.0*OnOffTimeRatio);
	
	// Use TFractionFitter to do fitting
	TObjArray *FractionFitData = new TObjArray(2);
	FractionFitData->Add(TheoryRecoilEvtHistogram);
	FractionFitData->Add(ReactorOffHisto);
	TFractionFitter* ffit = new TFractionFitter(ReactorOnHisto, FractionFitData);
	ffit->Constrain(0,0.1,10.0);
	ffit->Constrain(1,1.0,1.0);
	Int_t status = ffit->Fit();
	TH1F* result = (TH1F*) ffit->GetPlot();
	
	// Build a stacked histogram for plotting
	Double_t param, error;
	THStack *ReactorOnStackedFit = new THStack("ReactorOnStackedFit","Signal fits for Reactor-On data");
	ffit->GetResult(0,param,error);
	TheoryRecoilEvtHistogramFit->Add(TheoryRecoilEvtHistogram,param);
	ReactorOnStackedFit->Add(TheoryRecoilEvtHistogramFit);
	ReactorOnStackedFit->Add(ReactorOffHisto);
	
	
	
	// -------- MAKE PLOTS --------
	// Set drawing settings
	RecoilEvtHistogram->SetLineColor(1);
	NeutrinoEvtHistogram->SetLineColor(1);
	TheoryRecoilEvtHistogram->SetLineColor(2);
	TheoryNeutrinoEvtHistogram->SetLineColor(2);
	ReactorOnHisto->SetLineColor(1);
	result->SetLineColor(2);
	
	// Draw everything
	TCanvas* c1 = new TCanvas("c1");
	RecoilEvtHistogram->Draw("E1");
	TheoryRecoilEvtHistogram->Draw("same");
	legend1 = new TLegend(0.6,0.7,0.89,0.89);
  legend1->AddEntry(RecoilEvtHistogram,"Data","lep");
	legend1->AddEntry(TheoryRecoilEvtHistogram,"Monte Carlo","l");
	legend1->SetFillColor(0);
  legend1->Draw();
	
	TCanvas* c2 = new TCanvas("c2");
	NeutrinoEvtHistogram->Draw("E1");
	TheoryNeutrinoEvtHistogram->Draw("same");
	legend2 = new TLegend(0.6,0.7,0.89,0.89);
  legend2->AddEntry(RecoilEvtHistogram,"Data","lep");
	legend2->AddEntry(TheoryRecoilEvtHistogram,"Monte Carlo","l");
	legend2->SetFillColor(0);
  legend2->Draw();
	
	TCanvas* c3 = new TCanvas("c3");
	ReactorOnHisto->Draw("E1");
	ReactorOffHisto->Draw("E1,same");
	legend3 = new TLegend(0.6,0.7,0.89,0.89);
	legend3->AddEntry(ReactorOnHisto,"Reactor On","lep");
	legend3->AddEntry(ReactorOffHisto,"Reactor Off","lep");
	legend3->SetFillColor(0);
  legend3->Draw();
	
	TCanvas* c4 = new TCanvas("c4");
	ReactorNeutronBackground->Draw("E1");
	
	//This plot is broken: THStack is not being used correctly
	/*TCanvas* c5 = new TCanvas("c5");
	ReactorOnHisto->Draw("E1");
	TheoryRecoilEvtHistogramFit->SetFillColor(kRed);
	TheoryRecoilEvtHistogramFit->SetMarkerStyle(1);
	TheoryRecoilEvtHistogramFit->SetMarkerColor(kRed);
	ReactorOffHisto->SetFillColor(kBlue);
	ReactorOffHisto->SetMarkerStyle(1);
	ReactorOffHisto->SetMarkerColor(kBlue);
	ReactorOnStackedFit->Draw("same");
	legend5 = new TLegend(0.6,0.7,0.89,0.89);
	legend5->AddEntry(ReactorOnHisto,"Reactor On","lep");
	legend5->AddEntry(ReactorOffHisto,"Reactor Off");
	legend5->AddEntry(TheoryRecoilEvtHistogramFit,"Coherent Scattering");
	legend5->SetFillColor(0);
  legend5->Draw();*/
	
	WriteHistogramToFile(RecoilEvtHistogram, "histogramOutput.txt");
	WriteNuRecoilEvents("MonoenergeticEvents.txt", 0.811, nEvt, nNeutrons, nProtons, RecoilSpectrum, DiffRecoilSpectrumAtConstE, NuSpectrumMinE, NuSpectrumMaxE, seed+5);
	
}


// Function to generate the spectrum from coherent
// neutrino scattering.
void FillNuRecoilSpectrum(TH1F* recoilHisto, TH1F* nuEnergyHisto, Int_t nEvt, Int_t nNeutrons, Int_t nProtons, TF1* RecoilSpectrum, TF1* DiffRecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, Int_t seed)
{	
	Double_t NucleonMass = (nNeutrons * 939.565) + (nProtons * 938.272); // [MeV]
	
	// Now let's generate actually generate the events and
	// fill a histogram
	for (Int_t evtNum = 0; evtNum < nEvt; evtNum++) {
		// Get the neutrino energy
		Double_t thisNuEnergy = RecoilSpectrum->GetRandom(SpectrumMin, SpectrumMax);
		nuEnergyHisto->Fill(thisNuEnergy);
		
		// Get the recoil energy
		DiffRecoilSpectrum->SetParameter(3, thisNuEnergy);
		DiffRecoilSpectrum->SetParameter(4, MaxRecoilEnergy(thisNuEnergy, NucleonMass));
		Double_t thisRecoilEnergy = DiffRecoilSpectrum->GetRandom(0.0,MaxRecoilEnergy(thisNuEnergy, NucleonMass));
		recoilHisto->Fill(1000000*thisRecoilEnergy);
	}
}


// Function to generate the spectrum from coherent
// neutrino scattering for monoenergetic neutrinos
// and write events to file
void WriteNuRecoilEvents(const char* filename, Double_t nuEnergy, Int_t nEvt, Int_t nNeutrons, Int_t nProtons, TF1* RecoilSpectrum, TF1* DiffRecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, Int_t seed)
{
	FILE* outfile;
	outfile = fopen(filename, "w");
	Double_t NucleonMass = (nNeutrons * 939.565) + (nProtons * 938.272); // [MeV]
	
	// Now let's generate actually generate the events and
	// fill a histogram
	for (Int_t evtNum = 0; evtNum < nEvt; evtNum++)
	{
		// Get the recoil energy
		DiffRecoilSpectrum->SetParameter(3, nuEnergy);
		DiffRecoilSpectrum->SetParameter(4, MaxRecoilEnergy(nuEnergy, NucleonMass));
		Double_t thisRecoilEnergy = DiffRecoilSpectrum->GetRandom(0.0,MaxRecoilEnergy(nuEnergy, NucleonMass));
		fprintf(outfile, "%e\n", thisRecoilEnergy);
	}
}


// Generate the number of events for this run.
// This is just Poisson statistics and the cross section.
Int_t GenerateNumOfNuRecoils(Double_t time, Double_t detMass, Double_t distance, Double_t activity, Int_t nNeutrons, Int_t nProtons, TF1* RecoilSpectrum, Double_t SpectrumMin, Double_t SpectrumMax, UInt_t seed)
{
	// Set up the constants
	Double_t AvogadroConst = 6.022e23;
	Double_t NucleonMass = (nNeutrons * 939.565) + (nProtons * 938.272); // [MeV]
	Double_t molarMass = nNeutrons + nProtons;
	
	// Calculate mean number of events
	Double_t SpectrumWeightedCrossSection = RecoilSpectrum->Integral(SpectrumMin, SpectrumMax);
	Double_t MeanEvents = time * AvogadroConst * activity * detMass * SpectrumWeightedCrossSection / (4 * TMath::Pi() * molarMass * pow(distance,2));
	
	// Generate a random number of events from Poisson distribution
	TRandom3* randGen = new TRandom3(seed);
	Int_t nEvt = randGen->Poisson(MeanEvents);
	return nEvt;
}


// Read in an event spectrum from a text file containing
// spectral data, and then generate some random events.
// This assumes that the input data file contains two columns,
// the first with (evenly spaced!) energy bins in units of eV, and the second
// with a number of events per kg per day per keV
void FillRecoilSpectrumFromFile(TH1F* recoilHisto, Double_t time, Double_t mass, const char* filename, Int_t seed)
{
	FILE* dataFile;
	dataFile = fopen(filename,"r");
	if (dataFile == NULL) {
		cout << "The file " << filename << " is not found!!! Quitting now." << endl;
		exit(1);
	}
	
	// Count the number of lines to set the array sizes
	fstream lineCountStream;
	lineCountStream.open(filename, fstream::in);
	Int_t lineCount = 0;
	while (lineCountStream.peek() != EOF) {
		lineCountStream.ignore(128, '\n');
		lineCount++;
	}
	lineCountStream.close();
	cout << "The file has " << lineCount << " lines." << endl;
	
	// Should add a section here to double check the formatting of the file
	
	// Warning message about units
	cout << "Reading data and filling histogram from " << filename << "..." << endl;
	cout << "I am assuming that the first column in the file is in units of MeV "
			 << "and that the second column is in units of 1/(kg day keV)!" << endl << endl;
	
	// Read the data
	const Int_t nLines = lineCount;
	Float_t energies[nLines];
	Float_t rates[nLines];
	for (Int_t thisLine = 0; thisLine < nLines; thisLine++)
		fscanf(dataFile, "%E %E", &energies[thisLine], &rates[thisLine]);
	fclose(dataFile);
	
	// Set up the source histogram
	TH1F SpectralDataHistogram("SpectralDataHistogram", "Spectrum of Events from File", nLines-1, energies);
	for (Int_t thisBin = 0; thisBin < (nLines - 1); thisBin++) {
		SpectralDataHistogram.SetBinContent(thisBin,rates[thisBin]);
	}
	
	// Fill the recoil histogram with the correct number
	// of events, randomly chosen from the source histogram.
	TRandom3* randGen = new TRandom3(seed);
	Float_t meanNEvents = time * mass * SpectralDataHistogram.ComputeIntegral();
	Int_t nEvents = randGen->Poisson(meanNEvents);
	cout << "nEvents = " << nEvents << endl;
	for (Int_t thisEvent; thisEvent<nEvents; thisEvent++)
		recoilHisto->Fill(SpectralDataHistogram.GetRandom());
}


// Just computes the maximum recoil due to a neutrino
// of the given energy.  Note that all energies are in MeV.
Double_t MaxRecoilEnergy(Double_t NeutrinoEnergy, Double_t mTarget)
{
	return NeutrinoEnergy / (1 + mTarget / (2*NeutrinoEnergy));
}


// This writes the contents of a histogram to a text file,
// in case someone wants to analyze the data in another
// program/format.
void WriteHistogramToFile(TH1F* theHisto, const char* filename)
{
	FILE* outfile;
	outfile = fopen(filename, "w");
	
	Int_t nbins = theHisto->GetSize();
	
	Float_t binLowEdge, binContent;
	for(Int_t thisBin = 0; thisBin < nbins; thisBin++)
	{
		binLowEdge = theHisto->GetBinLowEdge(thisBin);
		binContent = theHisto->GetBinContent(thisBin);
		fprintf(outfile, "%d %d\n", binLowEdge, binContent);
		cout << binLowEdge << " " << binContent << endl;
	}
	fclose(outfile);
}