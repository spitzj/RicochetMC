////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVetoSD.hh                                        //
//  Description: sensitive detector class for veto scintillator       //
//               shield, modified from original RMCmain version      //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        25 January 2011                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef RMCVetoSD_h
#define RMCVetoSD_h 1

#include "G4VSensitiveDetector.hh"
#include "RMCgeometry/detectors/RMCVetoHit.hh"

class G4Step;
class G4HCofThisEvent;


class RMCVetoSD : public G4VSensitiveDetector
{
  public:
    RMCVetoSD(const G4String&);
    ~RMCVetoSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

  private:

    RMCVetoHitsCollection* vetoHitsCollection;
};

#endif
