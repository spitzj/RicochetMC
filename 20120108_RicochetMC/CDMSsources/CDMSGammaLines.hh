#ifndef CDMSGammaLines_hh
#define CDMSGammaLines_hh 1
////////////////////////////////////////////////////////////////////////
// $Id: CDMSGammaLines.hh,v 1.4 2011/07/06 04:10:44 kelsey Exp $
//  File:        CDMSGammaLines.hh                                    //
//  Description: Generate specificially listed gammas for CDMS        //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        20 May 2011                                          //
//                                                                    //
//  20110520  M. Kelsey -- Base class inspired by Am241Lines          //
//  20110524  M. Kelsey -- Add operator<<() for printing              //
//  20110616  M. Kelsey -- Inherit from new CDMSVParticleGenerator    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/CDMSVParticleGenerator.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <vector>
#include <iosfwd>


class CDMSGammaLines : public CDMSVParticleGenerator {
public:
  CDMSGammaLines(const G4String& name, const G4ThreeVector& dir=origin,
		 G4int verbose=0);

  CDMSGammaLines(const G4ThreeVector& dir=origin, G4int verbose=0);    

  virtual ~CDMSGammaLines();

  virtual G4double shootKineticEnergy() const;
  virtual G4ThreeVector shootDirection() const;

  // Unit vector for pencil-beam tests; set to (0,0,0) for isotropic outward
  void SetDirection(const G4ThreeVector& dir) { direction = dir.unit(); }
  const G4ThreeVector& GetDirection() const { return direction; }

  void PrintLines(std::ostream& os) const;	// Report spectrum data

protected:
  // Subclass must use this from ctor to fill spectrum
  void AddLine(G4double energy, G4double strength);

  void NormalizeLines();	// Set line strengths to sum to unity
  G4bool normalized() const { return lineSum.size() == lineEnergy.size(); }

private:
  G4ThreeVector direction;	// For generating pencil-beam sources

  std::vector<G4double> lineEnergy;	// Must be filled by subclass ctor
  std::vector<G4double> lineStrength;	// Must be filled by subclass ctor
  std::vector<G4double> lineSum;
};

// Global reporting (inlined) as call-throughs

inline std::ostream& 
operator<<(std::ostream& os, const CDMSGammaLines& theDet) {
  theDet.PrintLines(os);
  return os;
}

inline std::ostream& 
operator<<(std::ostream& os, const CDMSGammaLines* theDet) {
  if (theDet) theDet->PrintLines(os);
  return os;
}

#endif	/*  CDMSGammaLines_hh */
