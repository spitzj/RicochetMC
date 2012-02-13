////////////////////////////////////////////////////////////////////////
// $Id: CDMSGammaLines.cc,v 1.3 2011/06/17 05:07:12 kelsey Exp $
//  File:        CDMSGammaLines.cc                                    //
//  Description: Generate specificially listed gammas for CDMS        //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        20 May 2011                                          //
//                                                                    //
//  20110520  M. Kelsey -- Base class inspired by Am241Lines          //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/CDMSGammaLines.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>


// Constructor and destructor

CDMSGammaLines::CDMSGammaLines(const G4String& name, const G4ThreeVector& dir,
			       G4int verbose)
  : CDMSVParticleGenerator(name, G4Gamma::Definition()), direction(dir) {
  SetVerboseLevel(verbose);
}

CDMSGammaLines::CDMSGammaLines(const G4ThreeVector& dir, G4int verbose)
  : CDMSVParticleGenerator("CDMSGammaLines", G4Gamma::Definition()),
    direction(dir) {
  SetVerboseLevel(verbose);
}

CDMSGammaLines::~CDMSGammaLines() {}


// Add line to spectrum (not necessarily normalized)

void CDMSGammaLines::AddLine(G4double energy, G4double strength) {
  if (verboseLevel>2) {
    G4cout << GetName() << "::AddLine(" << energy/keV << " keV, "
	   << strength << ")" << G4endl;
  }

  lineEnergy.push_back(energy);
  lineStrength.push_back(strength);
}


// Generate cumulative intensity distribution, normalized to unity

void CDMSGammaLines::NormalizeLines() {
  if (normalized()) return;			// Normalization already done

  size_t NGLines = lineEnergy.size();
  if (NGLines == 0) {
    G4cerr << GetName() << ": ERROR: No gamma lines provided!" << G4endl;
    return;
  }

  if (verboseLevel>1) 
    G4cout << GetName() << "::NormalizeLines() with " << NGLines << " gammas"
	   << G4endl;

  if (NGLines != lineStrength.size()) {
    G4cerr << GetName() << ": ERROR: List of intensities different from"
	   << " list of lines." << G4endl;
    if (NGLines < lineStrength.size()) {
      G4cerr << " Truncating list of intensities." << G4endl;
    }

    lineStrength.resize(NGLines, 0.);	// Add zeroes or truncate excess
  }

  lineSum.resize(NGLines, 0.);		// Cumulative intensity distribution
  std::partial_sum(lineStrength.begin(), lineStrength.end(), lineSum.begin());

  if (verboseLevel>1) {
    G4cout << " cumulative distribution has " << lineSum.size() << " entries:"
	   << G4endl;
    std::copy(lineSum.begin(), lineSum.end(),
	      std::ostream_iterator<G4double>(G4cout, "\n"));
  }

  G4double sum = lineSum.back();	// Normalize everything to unity
  if (verboseLevel>1) G4cout << " got cumulative intensity " << sum << G4endl;

  std::transform(lineSum.begin(), lineSum.end(), lineSum.begin(),
		 std::bind2nd(std::divides<G4double>(), sum));

  std::transform(lineStrength.begin(), lineStrength.end(), lineStrength.begin(),
		 std::bind2nd(std::divides<G4double>(), sum));

  if (verboseLevel>1)			// Sanity check
    G4cout << " after normalizing, sum is " << lineSum.back() << G4endl;

  if (verboseLevel) PrintLines(G4cout);
}


// Generate random gamma energy by selecting a line

G4double CDMSGammaLines::shootKineticEnergy() const {
  if (verboseLevel>1) G4cout << GetName() << "::shootKineticEnergy()" << G4endl;

  if (!normalized()) {			// Ensure that strengths sum to unity
    const_cast<CDMSGammaLines*>(this)->NormalizeLines();
  }

  G4double eRand = G4UniformRand();	// Choose gamma from cumulative dist

  if (verboseLevel>2) G4cout << " generated eRand = " << eRand;

  std::vector<G4double>::const_iterator theLine
    = std::find_if(lineSum.begin(), lineSum.end(),
		   std::bind2nd(std::greater_equal<G4double>(), eRand));

  if (verboseLevel>2)
    G4cout << " selected entry " << (theLine-lineSum.begin()) << G4endl;

  // Get gamma energy by applying iterator offset to list of energies
  return *(lineEnergy.begin() + (theLine-lineSum.begin()));
}


// Generate a fixed or random direction for the gamma

G4ThreeVector CDMSGammaLines::shootDirection() const {
  if (verboseLevel > 1) G4cout << GetName() << "::shootDirection()" << G4endl;

  if (direction != origin) return direction;	// Avoid unnecessary work

  G4ThreeVector dir;
  G4double phi = 2.*pi*G4UniformRand();
  G4double theta = std::acos(2.*G4UniformRand() - 1.);	// -90,+90 deg

  if (verboseLevel>2)
    G4cout << " direction theta " << theta << " phi " << phi << G4endl;

  dir.setRThetaPhi(1., theta, phi);

  return dir;
}


// Report spectrum

void CDMSGammaLines::PrintLines(std::ostream& os) const {
  if (!normalized()) {			// Ensure that strengths sum to unity
    const_cast<CDMSGammaLines*>(this)->NormalizeLines();
  }

  os << GetName() << " with " << lineEnergy.size() << " gammas: " << std::endl;
  for (size_t i=0; i<lineEnergy.size(); i++) {
    os << " Line " << i << " : " << lineEnergy[i]/keV << " keV, "
       << lineStrength[i]/perCent << "% (" << lineSum[i] << ")" << std::endl;
  }
}
