////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCZipSD.hh                                         //
//  Description: sensitive detector class for zip detector, modified  //
//               from original RMCmini version                       //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        7 July 2010                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef RMCZipSD_h
#define RMCZipSD_h 1

#include "G4VSensitiveDetector.hh"
#include "RMCgeometry/detectors/RMCZipHit.hh"

class G4Step;
class G4HCofThisEvent;


class RMCZipSD : public G4VSensitiveDetector
{
  public:
    RMCZipSD(const G4String&);
    ~RMCZipSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

  private:

    RMCZipHitsCollection* zipHitsCollection;
};

#endif
