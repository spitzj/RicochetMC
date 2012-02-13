// $Id: CosmicMuon.hh,v 1.4 2011/07/22 19:33:58 kelsey Exp $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        CosmicMuon.hh                                                //
//  Description: generator of cosmic ray muons according to the distributions //
//               in Cassiday (PRD 7, 2022 (1973) ) which uses the Groom       //
//               parameterization                                             //
//  Author:      unknown, adapted for use in CDMS framework by                //
//               D.H. Wright (SLAC)                                           //
//  Date:        24 May 2011                                                  //
//                                                                            //
//  20110716  M. Kelsey -- Drop unnecessary functions, improve const-ness     //
////////////////////////////////////////////////////////////////////////////////

#ifndef SurfaceCosmicMuon_h
#define SurfaceCosmicMuon_h 1

#include "CDMSsources/CDMSVParticleGenerator.hh"
#include "CLHEP/Random/RandGeneral.h"
#include "globals.hh"


class CosmicMuon : public CDMSVParticleGenerator {
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

