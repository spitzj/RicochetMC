////////////////////////////////////////////////////////////////////////
//  File:        RMCExperimentConstruction.cc                         //     
//  Description: Full multiple tower detector with icebox             //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        12 February 2012                                     //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCgeometry/detectors/RMCExperimentConstruction.hh"
#include "RMCgeometry/detectors/RMCTowerConstruction.hh"
#include "RMCgeometry/detectors/RMCVesselConstruction.hh"
#include "RMCgeometry/detectors/RMCShieldConstruction.hh"
#include "RMCgeometry/interface/RMCExperimentMessenger.hh"
#include "RMCgeometry/interface/RMCGeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


// Constructor and destructor

RMCExperimentConstruction::RMCExperimentConstruction()
: RMCVDetectorGeometry("Experiment"),
theManager(RMCGeometryManager::Instance()),
towerBuilder(0), iceboxBuilder(0), shieldBuilder(0),
messenger(new RMCExperimentMessenger(this))
{
    //iceboxBuilder = theManager->GetVessel();
    //towerBuilder = theManager->GetTower();
}

RMCExperimentConstruction::~RMCExperimentConstruction() {
    delete messenger; messenger=0;
}


void RMCExperimentConstruction::SetVerboseLevel(G4int verbose) {
    RMCVDetectorGeometry::SetVerboseLevel(verbose);
    
    if (towerBuilder) towerBuilder->SetVerboseLevel(verbose);
    if (iceboxBuilder) iceboxBuilder->SetVerboseLevel(verbose);
    if (shieldBuilder) shieldBuilder->SetVerboseLevel(verbose);
}


void RMCExperimentConstruction::UseShield(G4bool shield) {
    if (verboseLevel)
        G4cout << "RMCExperimentConstruction::UseShield " << shield << G4endl;
    
    shieldBuilder = shield ? theManager->GetShield() : 0;
}


void RMCExperimentConstruction::UseCryostat(G4bool cryostat) {
    if (verboseLevel)
        G4cout << "RMCExperimentConstruction::UseCryostat " << cryostat << G4endl;
    
    iceboxBuilder = cryostat ? theManager->GetVessel() : 0;
}


void RMCExperimentConstruction::UseTower(G4bool tower) {
    if (verboseLevel)
        G4cout << "RMCExperimentConstruction::UseTower " << tower << G4endl;
    
    towerBuilder = tower ? theManager->GetTower() : 0;
}


void RMCExperimentConstruction::FillExtraParameters() {
    if (verboseLevel>1)
        G4cout << "RMCExperimentConstruction::FillExtraParameters()" << G4endl;
    
    SetVerboseLevel(verboseLevel);	// Make sure children have verbosity
    
    if(towerBuilder) towerBuilder->FillExtraParameters();
    if(iceboxBuilder) iceboxBuilder->FillExtraParameters();
    
    if (shieldBuilder) {
        shieldBuilder->UseIcebox();
        shieldBuilder->FillExtraParameters();
    }
}


G4LogicalVolume* RMCExperimentConstruction::BuildGeometry() {
    if (verboseLevel)
        G4cout << "RMCExperimentConstruction::BuildGeometry()" << G4endl;
    
    FillExtraParameters();	// Ensure constituents are up to date
    
    // Basic detector is cryostat enclosing a set of towers
    G4LogicalVolume* icebox = BuildVessel();
    //PlaceTowers(icebox);
    
    // If shielding requested, basic detector is placed inside it
    G4LogicalVolume* logicalFullDetector = 0;
    
    if (shieldBuilder) {
        logicalFullDetector = BuildShield();
        PlaceVesselInShield(icebox, logicalFullDetector);
    } else {
        logicalFullDetector = icebox;
    }
    PlaceTowers(icebox, logicalFullDetector);
    
    if (verboseLevel)
        G4cout << "Complete detector mass " << logicalFullDetector->GetMass()/kg
        << " kg" << G4endl;
    
    
    // Return either shielding or basic detector as "world volume"
    return logicalFullDetector;
}


G4LogicalVolume* RMCExperimentConstruction::BuildShield() {
    if (!shieldBuilder) return 0;
    
    if (verboseLevel>1)
        G4cout << "RMCExperimentConstruction::BuildShield()" << G4endl;
    
    return shieldBuilder->BuildGeometry();	// The shield is the world
}


G4LogicalVolume* RMCExperimentConstruction::BuildVessel() {
    if (!iceboxBuilder) return 0;
    
    if (verboseLevel>1)
        G4cout << "RMCExperimentConstruction::BuildVessel()" << G4endl;
    
    return iceboxBuilder->BuildGeometry();	// The vesselis the void
}


void RMCExperimentConstruction::PlaceVesselInShield(G4LogicalVolume* vessel,
                                                     G4LogicalVolume* shield) {
    if (!iceboxBuilder) return;
    
    if (verboseLevel>1)
        G4cout << "RMCExperimentConstruction::PlaceVessel()" << G4endl;
    
    new G4PVPlacement(0, iceboxBuilder->GetPosition(), vessel,
                      "PhysicalVessel", shield, false, 0);
}


void RMCExperimentConstruction::PlaceTowers(G4LogicalVolume* vessel,
                                             G4LogicalVolume* shield) {
    if (!towerBuilder) return;
    
    if (verboseLevel>1)
        G4cout << "RMCExperimentConstruction::PlaceTowers()" << G4endl;
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0)),
                      towerBuilder->BuildGeometryCopy(1),
                      "PhysicalTower", shield, false, 1);
}


// Vacuum vessel or shield is considered boundary of detector

G4double RMCExperimentConstruction::GetRadius() const {
    return shieldBuilder ? shieldBuilder->GetRadius() : iceboxBuilder->GetRadius();
}

G4double RMCExperimentConstruction::GetLength() const {
    return shieldBuilder ? shieldBuilder->GetLength() : iceboxBuilder->GetLength();
}
