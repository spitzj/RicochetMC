////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSourceMessenger.cc                         //
//  Description: User interface for cosmic muon generators            //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Dennis Wright (SLAC)                   //
//  Date:        13 February 2012                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "CosmogenicSourceMessenger.hh"
#include "CosmogenicSource.hh"

#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIdirectory.hh"

CosmogenicSourceMessenger::
CosmogenicSourceMessenger(CosmogenicSource* fPtclGun)
 : fParticleGun(fPtclGun)
{
  // create directory
  gunDirectory = new G4UIdirectory("/RMC/Cosmics/");
  gunDirectory->SetGuidance("Cosmogenic source control commands.");

  // depth of lab
  depthCmd = new G4UIcmdWithADoubleAndUnit("/RMC/Cosmics/depth",this);
  depthCmd->SetGuidance("Set depth of laboratory");
  depthCmd->SetParameterName("Z",true,true);
  depthCmd->SetDefaultUnit("m");

  // verbosity
  verbosityCmd = new G4UIcmdWithAnInteger("/RMC/Cosmics/verbose",this);
  verbosityCmd->SetGuidance("Set Verbose level for gun");
  verbosityCmd->SetGuidance(" 0 : Silent");
  verbosityCmd->SetGuidance(" 1 : Limited information");
  verbosityCmd->SetGuidance(" 2 : Detailed information");
  verbosityCmd->SetParameterName("level",false);
  verbosityCmd->SetRange("level>=0 && level <=2");

  // generator test
  genTestCmd = new G4UIcmdWithoutParameter("/RMC/Cosmics/gentest",this);
  genTestCmd->SetGuidance("Test cosmic muon generator");
}


CosmogenicSourceMessenger::~CosmogenicSourceMessenger() {
  delete depthCmd;
  delete verbosityCmd;
  delete genTestCmd;
  delete gunDirectory;
}


/*void CosmogenicSourceMessenger::SetNewValue(G4UIcommand* command,
                                            G4String newValues)
{
  if (command == depthCmd) {
    fParticleGun->SetLabDepth(depthCmd->GetNewDoubleValue(newValues));

  } else if(command == verbosityCmd) {
    fParticleGun->SetVerboseLevel(verbosityCmd->GetNewIntValue(newValues));

  } else if (command == genTestCmd) {
    fParticleGun->TestGenerator();

  } else {
    G4cout << "Error entering command" << G4endl;
  }
  }*/




