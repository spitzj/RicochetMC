// $Id: CosmogenicSourceMessenger.cc,v 1.2 2011/06/30 21:35:22 dwright Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSourceMessenger.cc                         //
//  Description: User interface for cosmic muon generators            //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        28 May 2011                                          //
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
  gunDirectory = new G4UIdirectory("/CDMS/Cosmics/");
  gunDirectory->SetGuidance("Cosmogenic source control commands.");

  // depth of lab
  depthCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Cosmics/depth",this);
  depthCmd->SetGuidance("Set depth of laboratory");
  depthCmd->SetParameterName("Z",true,true);
  depthCmd->SetDefaultUnit("m");

  // verbosity
  verbosityCmd = new G4UIcmdWithAnInteger("/CDMS/Cosmics/verbose",this);
  verbosityCmd->SetGuidance("Set Verbose level for gun");
  verbosityCmd->SetGuidance(" 0 : Silent");
  verbosityCmd->SetGuidance(" 1 : Limited information");
  verbosityCmd->SetGuidance(" 2 : Detailed information");
  verbosityCmd->SetParameterName("level",false);
  verbosityCmd->SetRange("level>=0 && level <=2");

  // generator test
  genTestCmd = new G4UIcmdWithoutParameter("/CDMS/Cosmics/gentest",this);
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




