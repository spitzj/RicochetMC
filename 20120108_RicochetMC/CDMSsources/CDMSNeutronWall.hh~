// $Id: CDMSGammaSphere.hh,v 1.6 2011/07/07 05:04:47 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGammaSphere.hh                                   //
//  Description: Inward spherical source for detector response        //
//                                                                    //
//  Author:      Michael Kelsey                                       //
//  Date:        23 May 2011                                          //
//                                                                    //
//  20110524  M. Kelsey -- Add configuration reporting.               //
//  20110705  M. Kelsey -- Use CDMSMultiGenerator for production, add //
//		opening angle around radial direction.                //
//  20110706  M. Kelsey -- Follow "GammaLines"->"Spectrum" change.    //
//  20110707  M. Kelsey -- Drop AddXXX() functions; no longer needed  //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSGammaSphere_hh
#define CDMSGammaSphere_hh 1

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "CDMSsources/CDMSMultiGenerator.hh"
#include "globals.hh"
#include <iosfwd>

class CDMSGammaSphereMessenger;
class G4Event;
class G4RadioactiveDecay;


class CDMSGammaSphere : public CDMSVSourceConstruction {
public:
  CDMSGammaSphere();
  virtual ~CDMSGammaSphere();

  // Pass verbosity down to generator
  virtual void SetVerboseLevel(G4int verbose=0) {
    CDMSVSourceConstruction::SetVerboseLevel(verbose);
    generator.SetVerboseLevel(verbose);
  }

  // Number of (independent) gammas to generate per event
  virtual void SetParticlesPerEvent(G4int n=1) { particlesPerEvent = n; }
  virtual G4int GetParticlesPerEvent() const   { return particlesPerEvent; }

  // Specify whether gammas are emitted inward (default) or outward
  virtual void SetRadialDirection(G4bool in=true)  { inwardGammas = in; }
  virtual void SetRadialDirection(G4double in=-1.) { inwardGammas = (in<0.); }
  virtual G4bool GetRadialDirection() const        { return inwardGammas; }

  virtual void SetHalfAngle(G4double val=0.*deg) { halfAngle = val; }
  virtual G4double GetHalfAngle() const { return halfAngle; }

  // Gammas are emitted uniformly from spherical surface
  virtual void SetRadius(G4double r=1.*meter) { radius = r; }
  virtual G4double GetRadius() const { return radius; }
  virtual G4double GetLength() const { return 2.*GetRadius(); }

  // Event generation, called automatically by RunManager
  virtual void GeneratePrimaries(G4Event* event);

  // Use CDMSVDetectorGeometry for diagnostic output
  virtual void PrintParameters(std::ostream& os) const;

protected:
  void GeneratePrimary(G4Event* event);		// Generate one gamma/decay
  void setGunPosition();			// Position on source sphere
  G4ThreeVector getGunDirection() const;	// Return radial unit vector

  // Acquire pointer to decay physics process in order to set collimation
  G4RadioactiveDecay* findRadioactiveDecayProcess() const;

private:
  CDMSMultiGenerator generator;

  G4bool inwardGammas;		// TRUE (default) if gammas point inward
  G4double radius;		// Radius of emitting sphere
  G4double halfAngle;		// Cone angle around radial direction

  G4int particlesPerEvent;

  CDMSGammaSphereMessenger* messenger;
};

#endif	/* CDMSGammaSphere_hh */
