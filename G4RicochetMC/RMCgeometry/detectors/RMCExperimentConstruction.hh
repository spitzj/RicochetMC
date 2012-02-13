#ifndef RMCExperimentConstruction_hh
#define RMCExperimentConstruction_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCExperimentConstruction.hh                         //     
//  Description: Full multiple tower detector with icebox             //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        12 February 2012                                     //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVDetectorGeometry.hh"

class G4LogicalVolume;
class RMCGeometryManager;
class RMCTowerConstruction;
class RMCVesselConstruction;
class RMCShieldConstruction;
class RMCExperimentMessenger;


class RMCExperimentConstruction: public RMCVDetectorGeometry {
public:
    RMCExperimentConstruction();
    virtual ~RMCExperimentConstruction();
    
    virtual G4double GetRadius() const;
    virtual G4double GetLength() const;
    virtual G4LogicalVolume* BuildGeometry();
    virtual void FillExtraParameters();
    
    virtual void SetVerboseLevel(G4int verbose=0);
    
    virtual void UseShield(G4bool shield=true);	// Build veto shield if true
    virtual void UseCryostat(G4bool shield=true);	// Build veto shield if true
    virtual void UseTower(G4bool shield=true);	// Build veto shield if true
    
private:
    G4LogicalVolume* BuildShield();
    G4LogicalVolume* BuildVessel();
    
    void PlaceVesselInShield(G4LogicalVolume* vessel, G4LogicalVolume* shield);
    void PlaceTowers(G4LogicalVolume* vessel, G4LogicalVolume* shield);
    
    RMCGeometryManager* theManager;
    RMCTowerConstruction* towerBuilder;
    RMCVesselConstruction* iceboxBuilder;
    RMCShieldConstruction* shieldBuilder;
    
    RMCExperimentMessenger* messenger;
};

#endif	/* RMCExperimentConstruction_hh */
