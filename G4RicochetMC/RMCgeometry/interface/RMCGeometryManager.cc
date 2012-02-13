////////////////////////////////////////////////////////////////////////
//  File:        RMCGeometryManager.cc                                //
//  Description: Factory to create "singleton" instances of geometry  //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCgeometry/interface/RMCGeometryManager.hh"


// Constructors and destructor

RMCGeometryManager::RMCGeometryManager()
: theShield(0), theTower(0), theTowerSupport(0),
 theVessel(0), theZip(0), theNoLab(0),
 theMITReactor(0) {}

RMCGeometryManager::~RMCGeometryManager() {
    delete theShield; theShield=0;
    delete theTower; theTower=0;
    delete theTowerSupport; theTowerSupport=0;
    delete theVessel; theVessel=0;
    delete theZip; theZip=0;
    delete theNoLab; theNoLab=0;
    delete theMITReactor; theMITReactor=0;
}

RMCGeometryManager* RMCGeometryManager::Instance() {
    static RMCGeometryManager* theManager = 0;
    if (!theManager) theManager = new RMCGeometryManager;
    return theManager;
}


// Create geometry components on demand
RMCShieldConstruction* RMCGeometryManager::GetShield() {
    if (!theShield) theShield = new RMCShieldConstruction;
    return theShield;
}

RMCTowerConstruction* RMCGeometryManager::GetTower() {
    if (!theTower) theTower = new RMCTowerConstruction;
    return theTower;
}

RMCTowerSupportConstruction* RMCGeometryManager::GetTowerSupport() {
    if (!theTowerSupport) theTowerSupport = new RMCTowerSupportConstruction;
    return theTowerSupport;
}

RMCVesselConstruction* RMCGeometryManager::GetVessel() {
    if (!theVessel) theVessel = new RMCVesselConstruction;
    return theVessel;
}

RMCZipConstruction* RMCGeometryManager::GetZip() {
    if (!theZip) theZip = new RMCZipConstruction;
    return theZip;
}

RMCExperimentConstruction* RMCGeometryManager::GetExperiment() {
    if (!theExperiment) theExperiment = new RMCExperimentConstruction;
    return theExperiment;
}


// Create labs on demand

NoLab* RMCGeometryManager::GetNoLab() {
    G4cout << "RMCGeometryManager::GetNoLab()" << G4endl;
    if (!theNoLab)
        theNoLab = new NoLab;
    
    return theNoLab;
}

RMCMITReactor* RMCGeometryManager::GetMITReactor() {
    G4cout << "RMCGeometryManager::GetMITReactor()" << G4endl;
    if (!theMITReactor)
        theMITReactor = new RMCMITReactor;
    
    return theMITReactor;
}
