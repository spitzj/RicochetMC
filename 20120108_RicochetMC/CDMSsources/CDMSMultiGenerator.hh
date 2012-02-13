////////////////////////////////////////////////////////////////////////
// $Id: CDMSMultiGenerator.hh,v 1.5 2011/07/07 05:04:47 kelsey Exp $
//  File:        CDMSMultiGenerator.hh                                //
//  Description: Base class to generate events from mixtures of       //
//		 sources, including both line spectra and decays.     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 July 2011                                          //
//                                                                    //
//  20110706  M. Kelsey -- Add weighting, "gammaLines" -> "spectra"   //
//  20110707  M. Kelsey -- Add separate SetHalfAngle() function, and  //
//		new Messenger.                                        //
//////////////////////////////////////////////////////////////////////// 

#ifndef CDMSMultiGenerator_hh
#define CDMSMultiGenerator_hh 1

#include "CDMSsources/CDMSVParticleGenerator.hh"
#include "globals.hh"
#include <vector>
#include <iosfwd>

class CDMSGammaLines;
class G4ParticleDefinition;
class CDMSMultiGeneratorMessenger;


class CDMSMultiGenerator : public CDMSVParticleGenerator {
public:
  CDMSMultiGenerator(const G4String& name="CDMSMultiGenerator");

  virtual ~CDMSMultiGenerator();

  // Build list of source spectra to be used (takes ownership)
  virtual void AddSpectrum(CDMSGammaLines* lines, G4double weight=1.);
  virtual void AddNucleus(G4int Z, G4int A, G4double weight=1.);
  virtual void AddNucleus(G4ParticleDefinition* isotope, G4double weight=1.);

  virtual size_t GetNSpectra() const { return spectra.size(); }
  virtual size_t GetNNuclei() const { return isotopes.size(); }
  virtual size_t GetNSources() const { return GetNSpectra()+GetNNuclei(); }

  // Configure collimation
  virtual void SetDirection(const G4ThreeVector& dir, G4double hAng=0.*deg);
  virtual void SetHalfAngle(G4double hAng) { halfAngle = hAng; }
  virtual const G4ThreeVector& GetDirection() const { return direction; }
  virtual G4double GetHalfAngle() const { return halfAngle; }

  // Functions to generate kinematic distributions
  virtual G4double shootKineticEnergy() const;
  virtual G4ThreeVector shootDirection() const;
  virtual G4ParticleDefinition* shootParticle() const;

  // Overload base-class version to do multi-source selection
  virtual const G4LorentzVector& shoot() const;

  // Redefine other accessor functions as call-throughs (to avoid "hiding")
  virtual void shoot(G4LorentzVector& mom, G4ParticleDefinition*& type) const {
    return CDMSVParticleGenerator::shoot(mom, type);
  }

  virtual void shoot(G4double& ekin, G4ThreeVector& dir,
		     G4ParticleDefinition*& type) const {
    return CDMSVParticleGenerator::shoot(ekin, dir, type);
  }

  // Report configuration
  virtual void PrintSources(std::ostream& os) const;

private:
  G4int chooseSource() const;	// Select which source to use for generation
  mutable G4int thisSource;

  G4ThreeVector direction;		// (0,0,0) for isotropic generation
  G4double halfAngle;			// ... for cone around direction

  G4bool useSpectrum(G4int i) const {
    return (i>=0 && (size_t)i<spectra.size());
  }

  void Normalize();			// Rescale weights to sum to unity

  std::vector<CDMSGammaLines*> spectra;		// List of spectra
  std::vector<G4ParticleDefinition*> isotopes;	// List of radionuclides
  std::vector<G4double> weights;		// Selection weights (all)
  std::vector<G4double> wtSum;			// Cumulative weights
  G4bool normalized;				// Flag if weights sum to 1.

  CDMSMultiGeneratorMessenger* messenger;
};

#endif	/* CDMSMultiGenerator_hh */
