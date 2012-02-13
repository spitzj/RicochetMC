#ifndef RMCPhysicsList_hh
#define RMCPhysicsList_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCPhysicsList.hh                                    //
//  Description: Physics processes relevant to backgrounds            //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               (adapted from Mike Kelsey (SLAC))                    //
//  Date:        8 January 2012                                       //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class RMCPhysicsMessenger;

class RMCPhysicsList: public G4VModularPhysicsList {
public:
  RMCPhysicsList(G4int verbose=0);
  virtual ~RMCPhysicsList();

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

  RMCPhysicsMessenger* messenger;
};

#endif	/* RMCPhysicsList_hh */
