// $Id: RMCPhysicsList.cc,v 1.5 2011/05/25 22:15:46 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCPhysicsList.cc                                   //
//  Description: Physics processes relevant to RMC backgrounds       //
//                                                                    //
//  Author:      Michael Kelsey                                       //
//  Date:        14 February 2011                                     //
//                                                                    //
//  20110322  M. Kelsey -- Add optical-photon physics.                //
//  20110428  M. Kelsey -- Per T. Koi, move builders to ConstructXXX  //
//  20110525  M. Kelsey -- Create lists separately, call              //
//		ConstructParticle for all of them.                    //
////////////////////////////////////////////////////////////////////////

#include "RMCsources/RMCPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "HadronPhysicsShielding.hh"
/*** No messenger yet; everything is hardwired
#include "RMCsources/RMCPhysicsMessenger.hh"
***/


// Constructor and destructor

RMCPhysicsList::RMCPhysicsList(G4int verbose)
  : G4VModularPhysicsList(), verboseLevel(verbose), 
    cutForGamma(defaultCutValue), cutForElectron(defaultCutValue),
    cutForPositron(defaultCutValue), cutForProton(defaultCutValue),
    particleList(0), emPhysicsList(0), hadPhysicsList(0),
    decayPhysicsList(0), opticalPhysicsList(0),
    messenger(0) {}

RMCPhysicsList::~RMCPhysicsList() {
  delete opticalPhysicsList;
  delete decayPhysicsList;
  delete hadPhysicsList;
  delete emPhysicsList;
  delete particleList;
  /*** No messenger yet; everything is hardwired
  delete messenger;
  ***/
}

// Must be created during running (NeutronHP), not in ctor above

void RMCPhysicsList::ConstructPhysicsLists() {
  if (!particleList)
    particleList = new G4DecayPhysics(verboseLevel);

  if (!emPhysicsList)
    emPhysicsList = new G4EmStandardPhysics(verboseLevel);

  if (!hadPhysicsList)
    hadPhysicsList = new HadronPhysicsShielding("Shielding",true);

  if (!decayPhysicsList)
    decayPhysicsList = new G4RadioactiveDecayPhysics(verboseLevel);

  if (!opticalPhysicsList) 
    opticalPhysicsList = new G4OpticalPhysics(verboseLevel);
}


// Build physics lists from components

void RMCPhysicsList::ConstructParticle() {
  ConstructPhysicsLists();

  particleList->ConstructParticle();
  emPhysicsList->ConstructParticle();
  hadPhysicsList->ConstructParticle();
  decayPhysicsList->ConstructParticle();
  opticalPhysicsList->ConstructParticle();
}

void RMCPhysicsList::ConstructProcess() {
  ConstructPhysicsLists();

  AddTransportation();				// Required!
  particleList->ConstructProcess();
  emPhysicsList->ConstructProcess();
  hadPhysicsList->ConstructProcess();
  decayPhysicsList->ConstructProcess();
  opticalPhysicsList->ConstructProcess();
}


// These functions and implementations copied from Hadr01/PhysicsList

void RMCPhysicsList::SetCuts() {
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  SetCutValue(cutForProton, "proton");

  if (verboseLevel) DumpCutValuesTable();
}
