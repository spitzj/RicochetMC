// $Id: CDMSMultiTowerConstruction.cc,v 1.15 2011/06/29 23:12:37 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSMultiTowerConstruction.cc                        //     
//  Description: Full multiple tower detector with icebox             //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        22 November 2010                                     //
//                                                                    //
//  20101130  M. Kelsey -- Override verbosity to pass through value   //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20101211  M. Kelsey -- Follow TowerPattern change to set origin   //
//  20110105  M. Kelsey -- Add "shield" here as detector component,   //
//		to be used as world volume enclosing vacuum vessel,   //
//		add Messenger in order to set shielding flag.         //
//  20110121  M. Kelsey -- Redo building so that towers are always    //
//		placed within cryostat volume; latter is optionally   //
//		placed within shielding.                              //
//  20110204  M. Kelsey -- Register Tower with TowerPattern in ctor   //
//  20110215  M. Kelsey -- Report mass at end of construction.        //
//  20110328  M. Kelsey -- Fetch components from Manager.             //
//  20110629  M. Kelsey -- Initialize and use name arg in base class  //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/detectors/CDMSMultiTowerConstruction.hh"
#include "CDMSgeometry/detectors/CDMSTowerPattern.hh"
#include "CDMSgeometry/detectors/CDMSTowerConstruction.hh"
#include "CDMSgeometry/detectors/CDMSVesselConstruction.hh"
#include "CDMSgeometry/detectors/CDMSShieldConstruction.hh"
#include "CDMSgeometry/interface/CDMSMultiTowerMessenger.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"


// Constructor and destructor

CDMSMultiTowerConstruction::CDMSMultiTowerConstruction()
: CDMSVDetectorGeometry("MultiTower"),
theManager(CDMSGeometryManager::Instance()),
towerPattern(0), towerBuilder(0), iceboxBuilder(0), shieldBuilder(0),
messenger(new CDMSMultiTowerMessenger(this)) {
    towerPattern = theManager->GetTowerPattern();
    towerBuilder = theManager->GetTower();
    iceboxBuilder = theManager->GetVessel();
}

CDMSMultiTowerConstruction::~CDMSMultiTowerConstruction() {
    delete messenger; messenger=0;
}


void CDMSMultiTowerConstruction::SetVerboseLevel(G4int verbose) {
    CDMSVDetectorGeometry::SetVerboseLevel(verbose);
    towerBuilder->SetVerboseLevel(verbose);
    towerPattern->SetVerboseLevel(verbose);
    iceboxBuilder->SetVerboseLevel(verbose);
    
    if (shieldBuilder) shieldBuilder->SetVerboseLevel(verbose);
}


void CDMSMultiTowerConstruction::UseShield(G4bool shield) {
    if (verboseLevel)
        G4cout << "CDMSMultiTowerConstruction::UseShield " << shield << G4endl;
    
    shieldBuilder = shield ? theManager->GetShield() : 0;
}


void CDMSMultiTowerConstruction::FillExtraParameters() {
    if (verboseLevel>1)
        G4cout << "CDMSMultiTowerConstruction::FillExtraParameters()" << G4endl;
    
    SetVerboseLevel(verboseLevel);	// Make sure children have verbosity
    
    towerBuilder->FillExtraParameters();
    iceboxBuilder->FillExtraParameters();
    
    if (shieldBuilder) {
        shieldBuilder->UseIcebox();
        shieldBuilder->FillExtraParameters();
    }
}


G4LogicalVolume* CDMSMultiTowerConstruction::BuildGeometry() {
    if (verboseLevel)
        G4cout << "CDMSMultiTowerConstruction::BuildGeometry()" << G4endl;
    
    FillExtraParameters();	// Ensure constituents are up to date
    
    // Basic detector is cryostat enclosing a set of towers
    G4LogicalVolume* icebox = BuildVessel();
    //PlaceTowers(icebox);
    
    // If shielding requested, basic detector is placed inside it
    G4LogicalVolume* logicalFullDetector = 0;
    
    if (shieldBuilder) {
        logicalFullDetector = BuildShield();
        //PlaceVesselInShield(icebox, logicalFullDetector);  //TEMPORARY!!
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


G4LogicalVolume* CDMSMultiTowerConstruction::BuildShield() {
    if (!shieldBuilder) return 0;
    
    if (verboseLevel>1)
        G4cout << "CDMSMultiTowerConstruction::BuildShield()" << G4endl;
    
    return shieldBuilder->BuildGeometry();	// The shield is the world
}


G4LogicalVolume* CDMSMultiTowerConstruction::BuildVessel() {
    if (verboseLevel>1)
        G4cout << "CDMSMultiTowerConstruction::BuildVessel()" << G4endl;
    
    return iceboxBuilder->BuildGeometry();	// The vesselis the void
}


void CDMSMultiTowerConstruction::PlaceVesselInShield(G4LogicalVolume* vessel,
                                                     G4LogicalVolume* shield) {
    if (verboseLevel>1)
        G4cout << "CDMSMultiTowerConstruction::PlaceVessel()" << G4endl;
    
    new G4PVPlacement(0, iceboxBuilder->GetPosition(), vessel,
                      "PhysicalVessel", shield, false, 0);
}


void CDMSMultiTowerConstruction::PlaceTowers(G4LogicalVolume* vessel,
                                             G4LogicalVolume* shield) {
    if (verboseLevel>1)
        G4cout << "CDMSMultiTowerConstruction::PlaceTowers()" << G4endl;
    
    //towerPattern->SetCenter(GetPosition());	// Towers at center of mother
    towerPattern->SetCenter(G4ThreeVector(0.,0.,0.));
    
    G4int nTowers = towerPattern->GetNumberOfTowers();

    for (G4int i=0; i<nTowers; i++) {
        if (verboseLevel > 2)
            G4cout << " placing tower #" << i << " of " << nTowers
            << " at " << towerPattern->GetPosition(i) << G4endl;
        
        new G4PVPlacement(towerPattern->GetTransform3D(i),
                          towerBuilder->BuildGeometryCopy(i),
                          "PhysicalTower", shield, false, i);
    }
}


// Vacuum vessel or shield is considered boundary of detector

G4double CDMSMultiTowerConstruction::GetRadius() const {
    return shieldBuilder ? shieldBuilder->GetRadius() : iceboxBuilder->GetRadius();
}

G4double CDMSMultiTowerConstruction::GetLength() const {
    return shieldBuilder ? shieldBuilder->GetLength() : iceboxBuilder->GetLength();
}
