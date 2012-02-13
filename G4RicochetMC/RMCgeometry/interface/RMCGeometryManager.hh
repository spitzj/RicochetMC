#ifndef RMCGeometryManager_hh
#define RMCGeometryManager_hh 1
////////////////////////////////////////////////////////////////////////
//  File:        RMCGeometryManager.hh                                //
//  Description: Factory to create "singleton" instances of geometry  //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

// NOTE:  Subclasses are defined here so clients have access to inheritance
#include "RMCgeometry/detectors/RMCShieldConstruction.hh"
#include "RMCgeometry/detectors/RMCTowerConstruction.hh"
#include "RMCgeometry/detectors/RMCTowerSupportConstruction.hh"
#include "RMCgeometry/detectors/RMCVesselConstruction.hh"
#include "RMCgeometry/detectors/RMCZipConstruction.hh"
#include "RMCgeometry/detectors/RMCExperimentConstruction.hh"
#include "RMCgeometry/labs/NoLab.hh"
#include "RMCgeometry/labs/RMCMITReactor.hh"


class RMCGeometryManager {
public:
    // Singleton constructor (naming parallels GEANT4 conventions)
    static RMCGeometryManager* Instance();
    static RMCGeometryManager* GetGeometryManager() { return Instance(); }
    
    virtual ~RMCGeometryManager();
    
    // Detector components
    RMCShieldConstruction* GetShield();
    RMCTowerConstruction* GetTower();
    RMCTowerSupportConstruction* GetTowerSupport();
    RMCVesselConstruction* GetVessel();
    RMCZipConstruction* GetZip();
    RMCExperimentConstruction* GetExperiment();
    
    // Lab facilities
    NoLab* GetNoLab();
    RMCMITReactor* GetMITReactor();
    
private:
    RMCGeometryManager();	// No one can create this except itself
    RMCGeometryManager(const RMCGeometryManager&);
    
    // Detector components
    RMCShieldConstruction* theShield;
    RMCTowerConstruction* theTower;
    RMCTowerSupportConstruction* theTowerSupport;
    RMCVesselConstruction* theVessel;
    RMCZipConstruction* theZip;
    RMCExperimentConstruction* theExperiment;
    
    // Caverns and test facilities
    NoLab* theNoLab;
    RMCMITReactor* theMITReactor;
};

#endif /* RMCGeometryManager_hh */
