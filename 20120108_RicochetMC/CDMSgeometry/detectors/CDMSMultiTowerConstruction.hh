#ifndef CDMSMultiTowerConstruction_hh
#define CDMSMultiTowerConstruction_hh 1
// $Id: CDMSMultiTowerConstruction.hh,v 1.9 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSMultiTowerConstruction.hh                        //     
//  Description: Full multiple tower detector with icebox             //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        22 November 2010                                     //
//                                                                    //
//  20101130  M. Kelsey -- Override verbosity to pass through value   //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20110105  M. Kelsey -- Add "shield" here as detector component,   //
//		add optional world-volume parameter to BuildVessel(), //
//		add Messenger in order to set shielding flag.         //
//  20110121  M. Kelsey -- Redo building so that towers are always    //
//		placed within cryostat volume; latter is optionally   //
//		placed within shielding.                              //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"

class G4LogicalVolume;
class CDMSGeometryManager;
class CDMSTowerPattern;
class CDMSTowerConstruction;
class CDMSVesselConstruction;
class CDMSShieldConstruction;
class CDMSMultiTowerMessenger;


class CDMSMultiTowerConstruction: public CDMSVDetectorGeometry {
public:
  CDMSMultiTowerConstruction();
  virtual ~CDMSMultiTowerConstruction();
  
  virtual G4double GetRadius() const;
  virtual G4double GetLength() const;
  virtual G4LogicalVolume* BuildGeometry();
  virtual void FillExtraParameters();

  virtual void SetVerboseLevel(G4int verbose=0);

  virtual void UseShield(G4bool shield=true);	// Build veto shield if true

private:
  G4LogicalVolume* BuildShield();
  G4LogicalVolume* BuildVessel();

  void PlaceVesselInShield(G4LogicalVolume* vessel, G4LogicalVolume* shield);
  void PlaceTowers(G4LogicalVolume* vessel, G4LogicalVolume* shield);

  CDMSGeometryManager* theManager;
  CDMSTowerPattern* towerPattern;
  CDMSTowerConstruction* towerBuilder;
  CDMSVesselConstruction* iceboxBuilder;
  CDMSShieldConstruction* shieldBuilder;

  CDMSMultiTowerMessenger* messenger;
};

#endif	/* CDMSMultiTowerConstruction_hh */
