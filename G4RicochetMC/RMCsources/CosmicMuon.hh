////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        CosmicMuon.hh                                                //
//  Description: generator of cosmic ray muons                                //
//                                                                            //
//  Author:      Adam Anderson (MIT)                                          //
//               Adapted from: D.H. Wright (SLAC)                             //
//  Date:        13 February 2012                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef SurfaceCosmicMuon_h
#define SurfaceCosmicMuon_h 1

#include "RMCsources/RMCVParticleGenerator.hh"
#include "CLHEP/Random/RandGeneral.h"
#include "globals.hh"


class CosmicMuon : public RMCVParticleGenerator {
public:
  CosmicMuon();
  virtual ~CosmicMuon() {;}
  
  
  // Generate kinematic distributions
  virtual G4double shootKineticEnergy() const;
  virtual G4ThreeVector shootDirection() const;
  virtual G4ParticleDefinition* shootParticle() const;

  G4double RandomTime(G4double area) const;
  
protected:
  G4double RandomPhi() const;
  G4double RandomTheta() const;
  G4double RandomEnergy() const;

  G4double MomentumDistribution(const G4double energy) const;
  G4double ThetaDistribution(const G4double theta) const;

private:
  CLHEP::RandGeneral* randGenEnergy;
  G4double binsize;
  G4double maxbin;
  G4double minbin;
  G4double* PDistFunc;
  G4double muonMass;
};

#endif

