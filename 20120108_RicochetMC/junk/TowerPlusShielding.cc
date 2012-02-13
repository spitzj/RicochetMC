// $Id: CDMSTowerConstruction.cc,v 1.15 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        TowerPlusShielding.hh                             //
//                                                                    //
//  Description: Construction of ZIP Tower with readout strings       //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        3 November 2010                                      //
//                                                                    //
//  20101130  M. Kelsey -- Override verbosity to pass through value.  //
//		Remove UseZipConstruction(), create in ctor.          //
//  20101203  M. Kelsey -- Make number of ZIPs per tower free param., //
//		Build housing as single structure here		      //
//  20101208  M. Kelsey -- Fix rotation of side strips along tower.   //
//  		Get housing parameters directly from ZIP.             //
//  20101209  M. Kelsey -- Fix bug in setting TowerLid radius arrays  //
//  20101210  M. Kelsey -- Fix calculation of towerClearanceR.        //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20110215  M. Kelsey -- Report mass at end of construction.        //
//  20110421  M. Kelsey -- Use Manager to access ZIP.                 //
//  20110629  M. Kelsey -- Initialize and use name arg in base class  //
//  20110722  M. Kelsey -- Discard matTable, using singleton instead. //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/detectors/TowerPlusShielding.hh"
#include "CDMSgeometry/detectors/CDMSShieldConstruction.hh"
#include "CDMSgeometry/detectors/CDMSZipConstruction.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSgeometry/interface/CDMSTowerMessenger.hh"
#include "CDMSg4base/CDMSMaterialTable.hh"
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyhedra.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include <cmath>

TowerPlusShielding::TowerPlusShielding()
  : theShield(new CDMSShieldConstruction()), theTower(new CDMSTowerConstruction()) {
}

TowerPlusShielding::~TowerPlusShielding() {
  delete theShield;
}


void TowerPlusShielding::SetVerboseLevel(G4int verbose) {
  CDMSVDetectorGeometry::SetVerboseLevel(verbose);
  if (theShield) theShield->SetVerboseLevel(verbose);
  if (theTower) theTower->SetVerboseLevel(verbose);
}


void TowerPlusShielding::FillExtraParameters() {
  if (!theShield || !theTower) return;		// If no ZIPs, no tower can be built

  theShield->FillExtraParameters();
  theTower->FillExtraParameters();
}


G4LogicalVolume* TowerPlusShielding::BuildGeometry() {
  if (verboseLevel)
    G4cout << "TowerPlusShielding::BuildGeometry()" << G4endl;

  FillExtraParameters();

  theShield->BuildGeometry();
  theTower->BuildGeometry();
}


void TowerPlusShielding::PrintParameters(std::ostream& os) const {
  theShield->PrintParameters(os);
  theTower->PrintParameters(os);
}
