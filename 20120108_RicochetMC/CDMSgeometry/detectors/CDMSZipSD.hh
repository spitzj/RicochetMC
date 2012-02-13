// $Id: CDMSZipSD.hh,v 1.4 2011/05/04 20:17:22 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSZipSD.hh                                         //
//  Description: sensitive detector class for zip detector, modified  //
//               from original CDMSmini version                       //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        7 July 2010                                          //
//                                                                    //
//  20110504  M. Kelsey -- Fix constness of ctor string.              //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSZipSD_h
#define CDMSZipSD_h 1

#include "G4VSensitiveDetector.hh"
#include "CDMSgeometry/detectors/CDMSZipHit.hh"

class G4Step;
class G4HCofThisEvent;


class CDMSZipSD : public G4VSensitiveDetector
{
  public:
    CDMSZipSD(const G4String&);
    ~CDMSZipSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

  private:

    CDMSZipHitsCollection* zipHitsCollection;
};

#endif
