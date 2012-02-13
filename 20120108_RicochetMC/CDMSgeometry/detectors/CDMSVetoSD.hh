// $Id: CDMSVetoSD.hh,v 1.2 2011/05/04 20:17:22 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVetoSD.hh                                        //
//  Description: sensitive detector class for veto scintillator       //
//               shield, modified from original CDMSmain version      //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        25 January 2011                                      //
//                                                                    //
//  20110504  M. Kelsey -- Fix constness of ctor string.              //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSVetoSD_h
#define CDMSVetoSD_h 1

#include "G4VSensitiveDetector.hh"
#include "CDMSgeometry/detectors/CDMSVetoHit.hh"

class G4Step;
class G4HCofThisEvent;


class CDMSVetoSD : public G4VSensitiveDetector
{
  public:
    CDMSVetoSD(const G4String&);
    ~CDMSVetoSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

  private:

    CDMSVetoHitsCollection* vetoHitsCollection;
};

#endif
