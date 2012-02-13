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
#include "globals.hh"


class SurfaceCosmicMuon : public CDMSVParticleGenerator {
public:
  SurfaceCosmicMuon();
  virtual ~SurfaceCosmicMuon() {;}
  
  void Test() const;
  void Initialize();
  
  // Generate kinematic distributions
  virtual G4double shootKineticEnergy() const;
  virtual G4ThreeVector shootDirection() const;
  virtual G4ParticleDefinition* shootParticle() const;

  G4double RandomTime(G4double area) const;

  virtual void SetMinEnergy(G4double muEMin) {muEnergyMin = muEMin;}
  virtual void SetMaxEnergy(G4double muEMax) {muEnergyMax = muEMax;}
  virtual void SetEnergyBins(G4int muEBins) {muEnergyBins = muEBins;}
  virtual void SetMWEDepth(G4double mwe) {MWEDepth = mwe;}
  
protected:
  G4double RandomPhi() const;
  G4double RandomTheta() const;
  G4double RandomEnergy() const;

  G4double ThetaEnergyDist(const G4double theta, const G4double energy) const;
  G4double Flux(const G4double min_energy, const G4double max_energy) const;
  G4double StoppingDistance(const G4double energy) const;
  G4double EnergyLossRate(const G4double energy) const;
  G4double ExponentialIntegral(const G4double x) const;
  
private:
  G4double FluxIndefiniteIntegral(const G4double energy, const G4double a,
				  const G4double c, const G4double d) const;
  G4double FunctionFMaximum(const G4double a) const;
  
  // To describe vertical flux of cosmic rays
  const G4double VerticalFlux_A;
  const G4double VerticalFlux_B;
  const G4double VerticalFlux_C;
  const G4double VerticalFlux_D;
  
  // To describe theta and energy distribution of cosmic rays
  const G4double ThetaEnergyDist_A;
  const G4double ThetaEnergyDist_B;
  G4double ThetaEnergyDist_Bh;
  const G4double ThetaEnergyDist_C;
  const G4double ThetaEnergyDist_D;
  G4double ThetaEnergyDist_Dh;
  
  // To describe theta distribution of cosmic rays
  const G4double ThetaDist_A;
  G4double ThetaDist_B;
  const G4double ThetaDist_C;
  G4double ThetaDist_D;
  
  // To describe energy distribution of cosmic rays
  G4double EnergyDist_A;
  const G4double EnergyDist_B;
  G4double EnergyDist_C;
  const G4double EnergyDist_D;
  
  // Loss rate is modeled piecewise linearly:  dE/dx = a(E + c)
  // For 1 < E < 10	(and E < 1)
  const G4double EnergyLossRate_1a;
  const G4double EnergyLossRate_1c;
  // For 10 < E < 100
  const G4double EnergyLossRate_10a;
  const G4double EnergyLossRate_10c;
  // For 100 < E < 1000
  const G4double EnergyLossRate_100a;
  const G4double EnergyLossRate_100c;
  // For 100 < E < 10000+
  const G4double EnergyLossRate_1000a;
  const G4double EnergyLossRate_1000c;
  
  // Additive constants for piecewise defined stopping distance
  const G4double StoppingDistance_0GeV;
  const G4double StoppingDistance_D1;
  const G4double StoppingDistance_10GeV;
  const G4double StoppingDistance_D10;
  const G4double StoppingDistance_100GeV;
  const G4double StoppingDistance_D100;
  const G4double StoppingDistance_1000GeV;
  const G4double StoppingDistance_D1000;
  
  const G4double muThetaMin;
  const G4double muThetaMax;
  const G4int muThetaBins;
  
  G4double muEnergyMin;
  G4double muEnergyMax;
  G4int muEnergyBins;
  
  // MWEDepth is absolute depth times ratio of rock density to water density
  G4double MWEDepth;
  
  // Total muon rate per unit horizontal area (flux)
  G4double TotalMuonFlux;
};

#endif	/* CosmicMuon_h */




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

/*#ifndef SurfaceCosmicMuon_h
#define SurfaceCosmicMuon_h 1

#include "CDMSsources/CDMSVParticleGenerator.hh"
#include "CLHEP/Random/RandGeneral.h"
#include "globals.hh"


class SurfaceCosmicMuon : public CDMSVParticleGenerator {
public:
  SurfaceCosmicMuon();
  virtual ~SurfaceCosmicMuon() {;}
  
  
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
  G4double* EDistFunc;
};

#endif
*/
