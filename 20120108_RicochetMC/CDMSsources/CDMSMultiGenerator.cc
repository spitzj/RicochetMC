////////////////////////////////////////////////////////////////////////
// $Id: CDMSMultiGenerator.cc,v 1.6 2011/07/07 05:04:47 kelsey Exp $
//  File:        CDMSMultiGenerator.hh                                //
//  Description: Base class to generate events from mixtures of       //
//		 sources, including both line spectra and decays.     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 July 2011                                          //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/CDMSMultiGenerator.hh"
#include "CDMSsources/CDMSMultiGeneratorMessenger.hh"
#include "CDMSsources/CDMSGammaLines.hh"
#include "G4Gamma.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>


// Constructor and destructor

CDMSMultiGenerator::CDMSMultiGenerator(const G4String& name)
  : CDMSVParticleGenerator(name), messenger(0) {
  messenger = new CDMSMultiGeneratorMessenger(this);
}

CDMSMultiGenerator::~CDMSMultiGenerator() {
  for (size_t i=0; i<GetNSpectra(); i++) delete spectra[i];
  spectra.clear();
  isotopes.clear();
  weights.clear();
  delete messenger;
}


// Build list of source spectra to be used (takes ownership)

void CDMSMultiGenerator::AddSpectrum(CDMSGammaLines* lines, G4double weight) {
  if (lines) {
    lines->SetVerboseLevel(verboseLevel);
    spectra.push_back(lines);
    weights.push_back(weight);
    normalized = false;
  }
}

// Build list of nuclei to be use with RadioactiveDecay production

void CDMSMultiGenerator::AddNucleus(G4int Z, G4int A, G4double weight) {
  static G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  AddNucleus(table->GetIon(Z, A, 0.*keV), weight);
}

void CDMSMultiGenerator::AddNucleus(G4ParticleDefinition* isotope,
				    G4double weight) {
  if (isotope) {
    isotopes.push_back(isotope);
    weights.push_back(weight);
    normalized = false;
  }
}


// Configure collimation with optional cone angle; origin means isotropic

void CDMSMultiGenerator::SetDirection(const G4ThreeVector& dir, G4double hAng) {
  if (dir == origin) {			// Isotropic generation
    direction.set(0.,0.,1.);
    halfAngle = 180.*deg;
  } else {
    direction = dir;
    halfAngle = hAng;
  }
}


// Generate kinematics and particle type by selection from sources

const G4LorentzVector& CDMSMultiGenerator::shoot() const {
  thisSource = chooseSource();  // Select which source to use for current event

  if (verboseLevel>1) 
    G4cout << GetName() << "::shoot() chooses source " << thisSource
	   << G4endl;

  return CDMSVParticleGenerator::shoot();	// Call through to base
}

G4double CDMSMultiGenerator::shootKineticEnergy() const {
  if (verboseLevel>1) 
    G4cout << GetName() << "::shootKineticEnergy" << G4endl;

  if (!useSpectrum(thisSource)) return 0.;		// RDM is at rest

  return spectra[thisSource]->shootKineticEnergy();
}

G4ThreeVector CDMSMultiGenerator::shootDirection() const {
  if (verboseLevel>1) G4cout << GetName() << "::shootDirection" << G4endl;

  // Generate direction vector within user-specified cone
  G4ThreeVector dir = direction;
  if (halfAngle > 0.) {
    // Generate uniform direction around central axis
    G4double phi = 2.*pi*G4UniformRand();
    G4double cosMin = std::cos(halfAngle);
    G4double cosTheta = (1.-cosMin)*G4UniformRand() + cosMin;	// [cosMin,1.)
    
    dir.setPhi(dir.phi()+phi);
    dir.setTheta(dir.theta()+std::acos(cosTheta));
  }

  return dir;
}

G4ParticleDefinition* CDMSMultiGenerator::shootParticle() const {
  if (verboseLevel>1) G4cout << GetName() << "::shootParticle" << G4endl;

  if (useSpectrum(thisSource)) return G4Gamma::Definition();

  return isotopes[thisSource-GetNSpectra()];
}


// Select source from lists for generating next event

G4int CDMSMultiGenerator::chooseSource() const {
  G4int n = GetNSources();
  if (0 == n) return -1;		// No decays?  No work to do

  if (!normalized) const_cast<CDMSMultiGenerator*>(this)->Normalize();

  if (verboseLevel>1) G4cout << GetName() << "::chooseSource" << G4endl;

  G4double rand = G4UniformRand();	// Choose source from cumulative dist

  if (verboseLevel>2) G4cout << " generated rand = " << rand;

  std::vector<G4double>::const_iterator src
    = std::find_if(wtSum.begin(), wtSum.end(),
		   std::bind2nd(std::greater_equal<G4double>(), rand));

  G4int index = src - wtSum.begin();	// Convert iterator to index offset
  if (verboseLevel>2) G4cout << " selected source " << index << G4endl;

  return index;
}


// Rescale source weights to sum to unity

void CDMSMultiGenerator::Normalize() {
  if (normalized) return;		// Avoid duplicate work
  else normalized = true;

  size_t nWeights = weights.size();
  if (nWeights == 0) {
    G4cerr << GetName() << ": ERROR: No sources provided!" << G4endl;
    return;
  }

  if (verboseLevel>1) 
    G4cout << GetName() << "::Normalize() with " << nWeights << G4endl;

  if (nWeights != GetNSources()) {
    G4cerr << GetName() << ": ERROR: List of weights different from"
	   << " list of sources." << G4endl;
    if (nWeights < GetNSources()) {
      G4cerr << " Truncating list of intensities." << G4endl;
    }

    weights.resize(nWeights, 0.);	// Add zeroes or truncate excess
  }

  wtSum.resize(nWeights, 0.);		// Cumulative intensity distribution
  std::partial_sum(weights.begin(), weights.end(), wtSum.begin());

  if (verboseLevel>1) {
    G4cout << " cumulative distribution has " << wtSum.size() << " entries: ";
    std::copy(wtSum.begin(), wtSum.end(),
	      std::ostream_iterator<G4double>(G4cout, " "));
    G4cout << G4endl;
  }

  G4double sum = wtSum.back();		// Normalize everything to unity
  if (verboseLevel>1) G4cout << " got cumulative intensity " << sum << G4endl;

  std::transform(wtSum.begin(), wtSum.end(), wtSum.begin(),
		 std::bind2nd(std::divides<G4double>(), sum));

  std::transform(weights.begin(), weights.end(), weights.begin(),
		 std::bind2nd(std::divides<G4double>(), sum));

  if (verboseLevel>1)			// Sanity check
    G4cout << " after normalizing, sum is " << wtSum.back() << G4endl;
}


// Report source configuration

void CDMSMultiGenerator::PrintSources(std::ostream& os) const {
  if (GetNSpectra() > 0) {
    os << GetName() << " " << GetNSpectra() << " discrete spectra:" << G4endl;
    for (size_t i=0; i<GetNSpectra(); i++) {
      os << " " << spectra[i]->GetName() << " (" << weights[i] << ")\n"
	 << spectra[i] << G4endl;
    }
  }

  if (GetNNuclei() > 0) {
    os << GetName() << " " << GetNNuclei() << " radionuclides:" << G4endl;
    for (size_t i=0; i<GetNNuclei(); i++) {
      os << " " << isotopes[i]->GetParticleName() << " ("
	 << weights[GetNSpectra()+i] << ")" << G4endl;
    }
  }
}
