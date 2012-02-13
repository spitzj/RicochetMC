////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        CosmicMuon.cc                                                //
//  Description: generator of cosmic ray muons                                //
//                                                                            //
//  Author:      Adam Anderson (MIT)                                          //
//               Adapted from: D.H. Wright (SLAC)                             //
//  Date:        13 February 2012                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include "CosmicMuon.hh"
#include "G4ThreeVector.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "CLHEP/Random/RandGeneral.h"
#include "Randomize.hh"
#include <math.h>




CosmicMuon::CosmicMuon() : RMCVParticleGenerator("CosmicMuon")
{
  G4int nbins = 2000;
  binsize = 0.5;
  minbin = 0.5;
  maxbin = 1000.;
  muonMass = 105.6584*MeV;
  G4double thisbin = minbin;  
  G4int binCount = 0;
  PDistFunc = new G4double[2000];

  while(thisbin < maxbin)
  {
    // My random generator computes the momentum of the muon
    // but we need the kinetic energy; do the conversion here
    PDistFunc[binCount] = MomentumDistribution(thisbin);
    thisbin = thisbin + binsize;
    binCount++;
  }

  randGenEnergy = new CLHEP::RandGeneral(PDistFunc,nbins);
}


// Generate random energy, direction, and muon type
G4double CosmicMuon::shootKineticEnergy() const {
  return RandomEnergy();
}

G4ThreeVector CosmicMuon::shootDirection() const {
  static G4ThreeVector dir;	// Buffer to avoid some memory churn
  dir.setRThetaPhi(1., RandomTheta(), RandomPhi());
  return dir;
}

G4ParticleDefinition* CosmicMuon::shootParticle() const {
  static G4ParticleDefinition* muplus  = G4MuonPlus::Definition();
  static G4ParticleDefinition* muminus = G4MuonMinus::Definition();

  return (G4UniformRand() < 0.5) ? muplus : muminus;
}


// This might be useful some day, when the measured
// rate is implemented.  Until then, this is useless.
// (FIXME!!)
G4double CosmicMuon::RandomTime(G4double area) const
{
  return 1.;
}


// Assume uniformly distributed polar angle
G4double CosmicMuon::RandomPhi() const
{
  // Returns a random phi from -pi to pi.
  return (G4UniformRand()-0.5) * twopi;
}


// Assume a simple cos^2(theta) distribution;
// use rejection sampling
G4double CosmicMuon::RandomTheta() const
{
  G4double thetaMin = 0;
  G4double thetaMax = twopi/4.;
  G4double Fmax = 1;
  G4double Fmin = 0;

  G4double theta = G4UniformRand() * (thetaMax - thetaMin);
  G4double F = G4UniformRand() * (Fmax - Fmin);
  G4int count = 0;
  while(F > ThetaDistribution(theta) && count < 1000)
  {
    theta = G4UniformRand() * (thetaMax - thetaMin);
    count++;
  }
  if(count == 1000)
    G4cout << "CosmicMuon::RandomTheta: "
	   << "Ran rejection sampling for 1000 iterations "
	   << "but did not find suitable point!!" << G4endl;
  return theta;
}



// Returns a random energy from parameterization in hep-ph/0102042
// Note that this does not take into account the zenith-angle
// dependence of the muon energy spectrum.  This should be added
// at some point. (FIXME!)
G4double CosmicMuon::RandomEnergy() const
{
  // Draw random number from energy distribution.
  // Note that this gives a number on [0,1] that we
  // must convert to an energy in MeV.
  G4double randNumber = randGenEnergy->fire();
  G4double randMomentum = ((randNumber * (maxbin - minbin)) + minbin)*GeV;

  //G4cout << "p,E: " << randMomentum/MeV << ", " << muonMass << ", " << sqrt(pow(randMomentum,2) + pow(muonMass,2)) - muonMass << G4endl;

  return sqrt(pow(randMomentum,2) + pow(muonMass,2)) - muonMass;
}


// This is equation (10) of hep-ph/0102042
G4double CosmicMuon::MomentumDistribution(const G4double momentum) const
{
  // Define parameters
  G4double H1 = 0.135;
  G4double H2 = -2.529;
  G4double H3 = -5.76;
  G4double S2 = -2.10;

  G4double y = log10(momentum);

  G4double H = H1 * (pow(y,3.)/2. - 5.*pow(y,2.)/2. + 3.*y) +
               H2 * (-2.*pow(y,3.)/3. + 3.*pow(y,2.) - 10.*y/3. + 1.) +
               H3 * (pow(y,3.)/6. - pow(y,2.)/2 + y/3.) +
               S2 * (pow(y,3.)/3. - 2.*pow(y,2.) + 11.*y/3. - 2.);

  return pow(10.,H);
}

G4double CosmicMuon::ThetaDistribution(const G4double theta) const
{
  return pow(cos(theta),2);
}

