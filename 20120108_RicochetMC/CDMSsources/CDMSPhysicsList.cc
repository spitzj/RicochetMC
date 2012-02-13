// $Id: CDMSPhysicsList.cc,v 1.5 2011/05/25 22:15:46 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSPhysicsList.cc                                   //
//  Description: Physics processes relevant to CDMS backgrounds       //
//                                                                    //
//  Author:      Michael Kelsey                                       //
//  Date:        14 February 2011                                     //
//                                                                    //
//  20110322  M. Kelsey -- Add optical-photon physics.                //
//  20110428  M. Kelsey -- Per T. Koi, move builders to ConstructXXX  //
//  20110525  M. Kelsey -- Create lists separately, call              //
//		ConstructParticle for all of them.                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/CDMSPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "HadronPhysicsShielding.hh"
/*** No messenger yet; everything is hardwired
#include "CDMSsources/CDMSPhysicsMessenger.hh"
***/


// Constructor and destructor

CDMSPhysicsList::CDMSPhysicsList(G4int verbose)
  : G4VModularPhysicsList(), verboseLevel(verbose), 
    cutForGamma(defaultCutValue), cutForElectron(defaultCutValue),
    cutForPositron(defaultCutValue), cutForProton(defaultCutValue),
    particleList(0), emPhysicsList(0), hadPhysicsList(0),
    decayPhysicsList(0), opticalPhysicsList(0),
    messenger(0) {}

CDMSPhysicsList::~CDMSPhysicsList() {
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

void CDMSPhysicsList::ConstructPhysicsLists() {
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

void CDMSPhysicsList::ConstructParticle() {
  ConstructPhysicsLists();

  particleList->ConstructParticle();
  emPhysicsList->ConstructParticle();
  hadPhysicsList->ConstructParticle();
  decayPhysicsList->ConstructParticle();
  opticalPhysicsList->ConstructParticle();
}

void CDMSPhysicsList::ConstructProcess() {
  ConstructPhysicsLists();

  AddTransportation();				// Required!
  particleList->ConstructProcess();
  emPhysicsList->ConstructProcess();
  hadPhysicsList->ConstructProcess();
  decayPhysicsList->ConstructProcess();
  opticalPhysicsList->ConstructProcess();
}


// These functions and implementations copied from Hadr01/PhysicsList

void CDMSPhysicsList::SetCuts() {
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  SetCutValue(cutForProton, "proton");

  if (verboseLevel) DumpCutValuesTable();
}
