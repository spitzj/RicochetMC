#ifndef CosmogenicSource_h
#define CosmogenicSource_h 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSource.hh                                  //
//  Description: Primary track generator using CosmicMuon generator   //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Dennis Wright (SLAC)                   //
//  Date:        13 February 2012                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////
#include "RMCg4base/RMCVSourceConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4ParticleDefinition;
class G4ParticleGun;
class CosmogenicSourceMessenger;
class CosmicMuon;


class CosmogenicSource : public RMCVSourceConstruction
{
  public:
    CosmogenicSource(G4String name = "cosmoSource");
    virtual ~CosmogenicSource();

    void GeneratePrimaries(G4Event* anEvent);

    virtual G4LogicalVolume* BuildGeometry();

  private:
    void GenerateTime();
    void GeneratePosition();
    void GenerateBeam();
    void FillParticleGun();

    G4double R_active;        // Radius of active region of source

    static G4int eventNumber;

    // Type of generated muon (mu+ or mu-)
    G4ParticleDefinition* partDef;

    // Initial time and space coordinates of generated muon 
    G4double eventTime;
    G4ThreeVector beamPos;

    // Initial direction and energy of generated muon
    G4ThreeVector beamDir;
    G4double beamEnergy;

    G4bool paramChange;
    G4double smearingXwidth;
    G4double smearingYwidth;
    G4double Esum;
    G4int NEsum;

    CosmicMuon* muonSpectrum;
    CosmogenicSourceMessenger* theMessenger;
};

#endif
