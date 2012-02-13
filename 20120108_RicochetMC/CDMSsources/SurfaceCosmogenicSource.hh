#ifndef SurfaceCosmogenicSource_h
#define SurfaceCosmogenicSource_h 1

// $Id: CosmogenicSource.hh,v 1.4 2011/07/22 21:08:36 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSource.hh                                  //
//  Description: Primary track generator using the CosmicMuon spectra //
//               from Cassiday et al., (PRD 7, 2022 (1973) )          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        27 May 2011                                          //
//                                                                    //
//  20110722  M. Kelsey -- Minor rewrite of generating functions to   //
//		make better use of CDMSVParticleGenerator interface.  //
//  20110722  M. Kelsey -- Remove materials table; use singleton.     //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4ParticleDefinition;
class G4ParticleGun;
class CosmogenicSourceMessenger;
class SurfaceCosmicMuon;


class SurfaceCosmogenicSource : public CDMSVSourceConstruction
{
  public:
    SurfaceCosmogenicSource(G4String name = "cosmoSource");
    virtual ~SurfaceCosmogenicSource();

    void GeneratePrimaries(G4Event* anEvent);

    virtual G4double GetRadius() const {return R_active;}
    virtual G4double GetLength() const {return L_active;}
    void SetRadius(G4double rad) {R_active = rad;} 

    virtual G4LogicalVolume* BuildGeometry();

    void ResetGenerator();

    G4int GetEventNumber();
    G4double GetEventTime();

    void SetLabDepth(G4double); 
    void TestGenerator();

    void SetRockDensity(G4double dens) {rockDensity = dens;}
    void SetOverBurden(G4double ovb) {overBurden = ovb;}
    void SetCavernLength(G4double len) {hallX = len;}
    void SetCavernWidth(G4double wid) {hallY = wid;}
    void SetCavernHeight(G4double hei) {hallZ = hei;}

  private:
    void GenerateTime();
    void GeneratePosition();
    void GenerateBeam();
    void FillParticleGun();
    G4bool InvalidDirection();

    G4double R_active;        // Radius of active region of source
    G4double L_active;        // Length of active region of source

    static G4int eventNumber;

    // Type of generated muon (mu+ or mu-)
    G4ParticleDefinition* partDef;

    // Initial time and space coordinates of generated muon 
    G4double eventTime;
    G4ThreeVector beamPos;

    // Initial direction and energy of generated muon
    G4ThreeVector beamDir;
    G4double beamEnergy;

    // Full length, width and height of cavern
    G4double hallX;
    G4double hallY; 
    G4double hallZ;
 
    G4double rockDensity; 
    G4double overBurden;      // thickness of rock above cavern in simulation
                              // != actual lab depth
    G4double muEffectRad;
    G4bool paramChange;
    G4double GenerationArea;

    // Effective range for passing muons
    G4double EffectiveXRange;
    G4double EffectiveYRange;
    G4double EffectiveZTravel;

    G4double labDepth;

    SurfaceCosmicMuon* muonSpectrum;
    CosmogenicSourceMessenger* theMessenger;

    G4double Esum;
    G4int NEsum;
};

#endif




/*#ifndef SurfaceCosmogenicSource_h
#define SurfaceCosmogenicSource_h 1

// $Id: CosmogenicSource.hh,v 1.4 2011/07/22 21:08:36 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSource.hh                                  //
//  Description: Primary track generator using the CosmicMuon spectra //
//               from Cassiday et al., (PRD 7, 2022 (1973) )          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        27 May 2011                                          //
//                                                                    //
//  20110722  M. Kelsey -- Minor rewrite of generating functions to   //
//		make better use of CDMSVParticleGenerator interface.  //
//  20110722  M. Kelsey -- Remove materials table; use singleton.     //
////////////////////////////////////////////////////////////////////////

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;
class G4ParticleDefinition;
class G4ParticleGun;
class CosmogenicSourceMessenger;
class SurfaceCosmicMuon;


class SurfaceCosmogenicSource : public CDMSVSourceConstruction
{
  public:
    SurfaceCosmogenicSource(G4String name = "cosmoSource");
    virtual ~SurfaceCosmogenicSource();

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

    SurfaceCosmicMuon* muonSpectrum;
    CosmogenicSourceMessenger* theMessenger;
};

#endif*/

