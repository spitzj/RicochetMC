#ifndef CDMSPhysicsList_hh
#define CDMSPhysicsList_hh 1
// $Id: CDMSPhysicsList.hh,v 1.4 2011/05/25 22:15:46 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSPhysicsList.hh                                   //
//  Description: Physics processes relevant to CDMS backgrounds       //
//                                                                    //
//  Author:      Michael Kelsey                                       //
//  Date:        14 February 2011                                     //
//                                                                    //
//  20110322  M. Kelsey -- Add optical-photon physics.                //
//  20110427  M. Kelsey -- Add verbose as constructor argument.       //
//  20110525  M. Kelsey -- Add lists' instantiation function.         //
////////////////////////////////////////////////////////////////////////

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class CDMSPhysicsMessenger;

class CDMSPhysicsList: public G4VModularPhysicsList {
public:
  CDMSPhysicsList(G4int verbose=0);
  virtual ~CDMSPhysicsList();

  virtual void ConstructParticle();
  virtual void ConstructProcess();

  void SetVerboseLevel(G4int vb=0) { verboseLevel = vb; }

  // These functions and implementations copied from Hadr01/PhysicsList
  void SetCuts();
  void SetCutForGamma(G4double cut)    { SetCutValue(cutForGamma=cut, "gamma"); }
  void SetCutForElectron(G4double cut) { SetCutValue(cutForElectron=cut, "e-"); }
  void SetCutForPositron(G4double cut) { SetCutValue(cutForPositron=cut, "e+"); }
  void SetCutForProton(G4double cut)   { SetCutValue(cutForProton=cut, "proton"); }

private:
  G4int verboseLevel;

  // Copied from Hadr01/PhysicsList, but don't know if necessary
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForProton;

  // Need to carry around these pointers so they can be deleted at end
  void ConstructPhysicsLists();

  G4VPhysicsConstructor*  particleList;
  G4VPhysicsConstructor*  emPhysicsList;
  G4VPhysicsConstructor*  hadPhysicsList;
  G4VPhysicsConstructor*  decayPhysicsList;
  G4VPhysicsConstructor*  opticalPhysicsList;

  CDMSPhysicsMessenger* messenger;
};

#endif	/* CDMSPhysicsList_hh */
