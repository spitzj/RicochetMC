#ifndef CDMSGeometryManager_hh
#define CDMSGeometryManager_hh 1
// $Id: CDMSGeometryManager.hh,v 1.4 2011/05/05 00:30:27 dwright Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGeometryManager.hh                               //
//  Description: Factory to create "singleton" instances of geometry  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        14 February 2011                                     //
//                                                                    //
//  20110421  M. Kelsey -- Add Am-241 source holder                   //
//////////////////////////////////////////////////////////////////////// 

// NOTE:  Subclasses are defined here so clients have access to inheritance
#include "CDMSgeometry/detectors/CDMSMiniDetConstruction.hh"
#include "CDMSgeometry/detectors/CDMSMultiTowerConstruction.hh"
#include "CDMSgeometry/detectors/CDMSShieldConstruction.hh"
#include "CDMSgeometry/detectors/CDMSTowerConstruction.hh"
#include "CDMSgeometry/detectors/CDMSTowerPattern.hh"
#include "CDMSgeometry/detectors/CDMSTowerSupportConstruction.hh"
#include "CDMSgeometry/detectors/CDMSVesselConstruction.hh"
#include "CDMSgeometry/detectors/CDMSZipConstruction.hh"
#include "CDMSgeometry/detectors/CKGDetConstruction.hh"
#include "CDMSgeometry/detectors/Am241SourceHolder.hh"
#include "CDMSgeometry/labs/NoLab.hh"
#include "CDMSgeometry/labs/SnoLab.hh"
#include "CDMSgeometry/labs/MITReactor.hh"


class CDMSGeometryManager {
public:
  // Singleton constructor (naming parallels GEANT4 conventions)
  static CDMSGeometryManager* Instance();
  static CDMSGeometryManager* GetGeometryManager() { return Instance(); }

  virtual ~CDMSGeometryManager();

  // Detector components
  CDMSMultiTowerConstruction* GetMultiTower();
  CDMSShieldConstruction* GetShield();
  CDMSTowerConstruction* GetTower();
  CDMSTowerPattern* GetTowerPattern();
  CDMSTowerSupportConstruction* GetTowerSupport();
  CDMSVesselConstruction* GetVessel();
  CDMSZipConstruction* GetZip();
  Am241SourceHolder* GetAm241Holder();
  CDMSMiniDetConstruction* GetMiniDet();
  CKGDetConstruction* GetCKGDet();

  // Caverns and test facilities
  NoLab* GetNoLab();
  SnoLab* GetSnoLab();
  MITReactor* GetMITReactor();

private:
  CDMSGeometryManager();	// No one can create this except itself
  CDMSGeometryManager(const CDMSGeometryManager&);

  // Detector components
  CDMSMultiTowerConstruction* theMultiTower;
  CDMSShieldConstruction* theShield;
  CDMSTowerConstruction* theTower;
  CDMSTowerPattern* theTowerPattern;
  CDMSTowerSupportConstruction* theTowerSupport;
  CDMSVesselConstruction* theVessel;
  CDMSZipConstruction* theZip;
  CDMSMiniDetConstruction* theMiniDet;
  CKGDetConstruction* theCKGDet;
  Am241SourceHolder* theAm241;

  // Caverns and test facilities
  NoLab* theNoLab;
  SnoLab* theSnoLab;
  MITReactor* theMITReactor;
};

#endif /* CDMSGeometryManager_hh */
