////////////////////////////////////////////////////////////////////////
// $Id: CDMS_SelectMISS.cc,v 1.2 2011/07/21 20:22:17 kevmc Exp $
//                                                                    //
//  File:        CDMS_SelectMISS.cc                                   //
//  Description: Find events containing multiple scatters per crystal //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        8 July 2011                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSactions/CDMS_SelectMISS.hh"
#include "CDMSgeometry/detectors/CDMSZipHit.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4VHitsCollection.hh"


// Select events with multiple hits from same track

G4bool CDMS_SelectMISS::accept(const G4Event* evt) const {
  if (verboseLevel > 1) G4cout << GetName() << "::accept" << G4endl;

  CDMSZipHitsCollection* ZHC = getHits(evt);
  if (!ZHC || ZHC->entries() == 0) return false;	// No hits, no work

  // Scan through hits to find multiples from same track
	G4int nZipHits = ZHC->entries();
	G4double bottomBoundary = 0.001;
	G4double topBoundary = 0.009;
	
	G4bool hasBottomHit = false;
	G4bool hasTopHit = false;
	
	for (G4int i = 0; i < nZipHits; i++) {
	  G4double z3 = (*ZHC)[i]->GetPostStepPosition().z()/meter;
		hasBottomHit = hasBottomHit || (z3 < bottomBoundary);
		hasTopHit = hasTopHit || (z3>topBoundary);
	}
	

  return (hasTopHit && hasBottomHit);		  // No duplicates found
}


// Extract hits from event for analysis

CDMSZipHitsCollection* CDMS_SelectMISS::getHits(const G4Event* evt) const {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4int zipCollID = SDman ? SDman->GetCollectionID("zipHitsCollection") : -1;

  if (zipCollID < 0) return 0;			// Avoid unnecessary work

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  return (HCE ? dynamic_cast<CDMSZipHitsCollection*>(HCE->GetHC(zipCollID)) : 0);
}
