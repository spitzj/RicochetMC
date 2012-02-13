#ifndef OscillationFitter_h
#define OscillationFitter_h 1
#include "ReactorExperiment.h"
#include <TH1F.h>
#include <TH2F.h>

class OscillationFitter
{
public:
    OscillationFitter();
    ~OscillationFitter();
    
    void SetFitter(ReactorExperiment* theExpt, TH1F* eSpectrum);
    void GlobalScanFit(const Int_t nBinsDM2, Double_t minDM2, Double_t maxDM2, const Int_t nBinsSin22ThetaS, Double_t minSin22ThetaS, Double_t maxSin22ThetaS, Int_t seed);
    TH2F* GetGlobalScanContour();
    //void RasterScanFit();
    //void FeldmanCousinsFit();
    
private:
    ReactorExperiment* theExperiment;
    TH2F* GlobalScanContour;
    TH1F* energySpectrum;
};
#endif 