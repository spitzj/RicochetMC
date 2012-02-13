#ifndef RMCVParticleGenerator_hh
#define RMCVParticleGenerator_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVParticleGenerator.hh                            //
//  Description: base class for generating background particles with  //
//               specific characteristics (radiogenic or cosmogenic)  //
//                                                                    //
//  Client code: Configure particle type in ctor or with SetParticle. //
//               On each event, get back 4-vector and particle type.  //
//                                                                    //
//  Subclasses:  Implement functions to throw energy, direction.      //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        16 June 2011                                         //
//                                                                    //
//  20110705  M. Kelsey -- Make shootXXX() functions public.          //
//////////////////////////////////////////////////////////////////////// 

#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleDefinition;


class RMCVParticleGenerator {
public:
  RMCVParticleGenerator(const G4String& name="RMCVParticleGenerator",
			 G4ParticleDefinition* type=0);
  virtual ~RMCVParticleGenerator();

  // Subclasses MUST implement functions to generate kinematic distributions
  virtual G4double shootKineticEnergy() const = 0;
  virtual G4ThreeVector shootDirection() const = 0;

  // Subclasses MAY implement function to call SetParticle at generation
  virtual G4ParticleDefinition* shootParticle() const { return particle; }

  virtual void SetVerboseLevel(G4int verbose) { verboseLevel = verbose; }

  // Configure particle used for next generated event
  virtual void SetParticle(G4ParticleDefinition* type) { particle = type; }
  virtual G4ParticleDefinition* GetParticle() const { return particle; }

  // Main function to generate a random particle; subclass may overload
  // provided they call through to base
  virtual const G4LorentzVector& shoot() const;

  // Accessor functions to get values for configuring G4ParticleGun
  virtual void shoot(G4LorentzVector& mom, G4ParticleDefinition*& type) const;

  virtual void shoot(G4double& ekin, G4ThreeVector& dir,
		     G4ParticleDefinition*& type) const;

  virtual const G4LorentzVector& GetMomentum() const { return lastShoot; }

  virtual const G4String& GetName() const { return sourceName; }

private:
  G4String sourceName;

  mutable G4ParticleDefinition* particle;	// Can overwrite each event
  mutable G4LorentzVector lastShoot;

protected:
  G4int verboseLevel;			// For subclass convenience

  static const G4ThreeVector origin;	// For convenience reference (0,0,0)
};

#endif	/* RMCVParticleGenerator_hh */
