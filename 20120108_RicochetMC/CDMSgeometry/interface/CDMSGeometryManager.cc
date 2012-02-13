// $Id: CDMSGeometryManager.cc,v 1.5 2011/06/17 20:03:55 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGeometryManager.cc                               //
//  Description: Factory to create "singleton" instances of geometry  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        14 February 2011                                     //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/interface/CDMSGeometryManager.hh"


// Constructors and destructor

CDMSGeometryManager::CDMSGeometryManager()
  : theMultiTower(0), theShield(0), theTower(0), theTowerPattern(0),
    theTowerSupport(0), theVessel(0), theZip(0), theMiniDet(0),
    theCKGDet(0), theAm241(0), theNoLab(0), theSnoLab(0),
    theMITReactor(0) {}

CDMSGeometryManager::~CDMSGeometryManager() {
  delete theMultiTower; theMultiTower=0;
  delete theShield; theShield=0;
  delete theTower; theTower=0;
  delete theTowerPattern; theTowerPattern=0;
  delete theTowerSupport; theTowerSupport=0;
  delete theVessel; theVessel=0;
  delete theZip; theZip=0;
  delete theMiniDet; theMiniDet=0;
  delete theCKGDet; theCKGDet=0;
  delete theAm241; theAm241=0;
  delete theNoLab; theNoLab=0;
  delete theSnoLab; theSnoLab=0;
  delete theMITReactor; theMITReactor=0;
}

CDMSGeometryManager* CDMSGeometryManager::Instance() {
  static CDMSGeometryManager* theManager = 0;
  if (!theManager) theManager = new CDMSGeometryManager;
  return theManager;
}


// Create geometry components on demand

CDMSMultiTowerConstruction* CDMSGeometryManager::GetMultiTower() {
  if (!theMultiTower) theMultiTower = new CDMSMultiTowerConstruction;
  return theMultiTower;
}

CDMSShieldConstruction* CDMSGeometryManager::GetShield() {
  if (!theShield) theShield = new CDMSShieldConstruction;
  return theShield;
}

CDMSTowerConstruction* CDMSGeometryManager::GetTower() {
  if (!theTower) theTower = new CDMSTowerConstruction;
  return theTower;
}

CDMSTowerPattern* CDMSGeometryManager::GetTowerPattern() {
  if (!theTowerPattern) theTowerPattern = new CDMSTowerPattern(theTower);
  return theTowerPattern;
}

CDMSTowerSupportConstruction* CDMSGeometryManager::GetTowerSupport() {
  if (!theTowerSupport) theTowerSupport = new CDMSTowerSupportConstruction;
  return theTowerSupport;
}

CDMSVesselConstruction* CDMSGeometryManager::GetVessel() {
  if (!theVessel) theVessel = new CDMSVesselConstruction;
  return theVessel;
}

CDMSZipConstruction* CDMSGeometryManager::GetZip() {
  if (!theZip) theZip = new CDMSZipConstruction;
  return theZip;
}

Am241SourceHolder* CDMSGeometryManager::GetAm241Holder() {
  if (!theAm241) theAm241 = new Am241SourceHolder;
  return theAm241;
}


// These detector structures are obsolete

CDMSMiniDetConstruction* CDMSGeometryManager::GetMiniDet() {
  if (!theMiniDet) theMiniDet = new CDMSMiniDetConstruction;
  return theMiniDet;
}

CKGDetConstruction* CDMSGeometryManager::GetCKGDet() {
  if (!theCKGDet) theCKGDet = new CKGDetConstruction;
  return theCKGDet;
}


// Create caverns and test facilities on demand

NoLab* CDMSGeometryManager::GetNoLab() {
  if (!theNoLab) theNoLab = new NoLab;
  return theNoLab;
}

SnoLab* CDMSGeometryManager::GetSnoLab() {
  if (!theSnoLab) theSnoLab = new SnoLab;
  return theSnoLab;
}

MITReactor* CDMSGeometryManager::GetMITReactor() {
  G4cout << "CDMSGeometryManager::GetMITReactor()" << G4endl;
  if (!theMITReactor) theMITReactor = new MITReactor;
  return theMITReactor;
}
