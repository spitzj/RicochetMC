// $Id: CKGZipSD.hh,v 1.2 2011/06/29 23:12:37 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CKGZipSD.hh                                          //
//  Description: sensitive detector class for zip detector, for use   //
//               CKG configurable detector                            //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        4 November 2010                                      //
//                                                                    //
//  20110629  M. Kelsey -- Past const-ref string.                     //
////////////////////////////////////////////////////////////////////////

#ifndef CKGZipSD_h
#define CKGZipSD_h 1

#include "G4VSensitiveDetector.hh"
#include "CDMSgeometry/detectors/CDMSZipHit.hh"

class G4Step;
class G4HCofThisEvent;


class CKGZipSD : public G4VSensitiveDetector
{
  public:
    CKGZipSD(const G4String&, G4int, G4int);
    ~CKGZipSD();

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    void EndOfEvent(G4HCofThisEvent*);

  private:

    CDMSZipHitsCollection* zipHitsCollection;
    G4int NTowers;
    G4int NZipsPerTower;
};

#endif
