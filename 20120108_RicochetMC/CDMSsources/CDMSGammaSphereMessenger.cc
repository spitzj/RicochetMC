////////////////////////////////////////////////////////////////////////
//  File:        CDMSGammaSphereMessenger.cc                          //
//  Description: messenger class for spherical gamma source           //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        24 May 2011                                          //
//                                                                    //
//  20110524  M. Kelsey -- Forgot to call base SetNewValue()          //
//  20110525  M. Kelsey -- Add some diagnostic messages.              //
//  20110707  M. Kelsey -- Move source addition to MultiGenMessenger  //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/CDMSGammaSphereMessenger.hh"
#include "CDMSsources/CDMSGammaSphere.hh"
#include "globals.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

// Constructor and destructor

CDMSGammaSphereMessenger::
CDMSGammaSphereMessenger(CDMSGammaSphere* theSource)
  : CDMSMessengerBase("/CDMS/GammaSphere/", "Uniform spherical gammas"),
    source(theSource), RadiusCmd(0), CenterCmd(0), InwardCmd(0),
    HalfAngleCmd(0), NParticlesCmd(0) {
  RadiusCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("Radius",
	"Radius of source sphere");
  RadiusCmd->SetParameterName("Size",false);
  RadiusCmd->SetRange("Size>=0.");
  RadiusCmd->SetUnitCategory("Length");  
  RadiusCmd->SetDefaultValue(source->GetRadius());
  RadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CenterCmd = CreateCommand<G4UIcmdWith3VectorAndUnit>("Position",
	"Center of source sphere in world coordinates");
  CenterCmd->SetParameterName("X", "Y", "Z", false);
  CenterCmd->SetUnitCategory("Length");
  CenterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  InwardCmd = CreateCommand<G4UIcmdWithABool>("Direction",
	"Emit gammas inward (true, 1) or outward (false, 0)");
  InwardCmd->SetParameterName("inward", true, false);
  InwardCmd->SetDefaultValue(source->GetRadialDirection());
  InwardCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HalfAngleCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("HalfAngle",
	"Half angle of cone around radial direction");
  HalfAngleCmd->SetParameterName("angle", true, false);
  HalfAngleCmd->SetUnitCategory("Angle");
  HalfAngleCmd->SetDefaultValue(source->GetHalfAngle());
  HalfAngleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  NParticlesCmd = CreateCommand<G4UIcmdWithAnInteger>("Gammas",
	"Number of gammas to be generated per event");
  NParticlesCmd->SetParameterName("N",true,true);
  NParticlesCmd->SetDefaultValue(source->GetParticlesPerEvent());
  NParticlesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

CDMSGammaSphereMessenger::~CDMSGammaSphereMessenger() {
  delete RadiusCmd;
  delete CenterCmd;
  delete InwardCmd;
  delete HalfAngleCmd;
  delete NParticlesCmd;
}


// Apply user command inputs

void CDMSGammaSphereMessenger::ActionAfterSetVerbose() {
  if (source) source->SetVerboseLevel(verboseLevel);
}

void 
CDMSGammaSphereMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (verboseLevel > 2) {
    G4cout << "CDMSGammaSphereMessenger::SetNewValue " << cmd->GetCommandPath()
	   << " " << value << G4endl;
  }

  CDMSMessengerBase::SetNewValue(cmd, value);	// To process global commands

  if (cmd == RadiusCmd)
    source->SetRadius(RadiusCmd->GetNewDoubleValue(value));

  if (cmd == CenterCmd)
    source->SetPosition(CenterCmd->GetNew3VectorValue(value));

  if (cmd == InwardCmd)
    source->SetRadialDirection(InwardCmd->GetNewBoolValue(value));

  if (cmd == HalfAngleCmd)
    source->SetHalfAngle(HalfAngleCmd->GetNewDoubleValue(value));

  if (cmd == NParticlesCmd)
    source->SetParticlesPerEvent(NParticlesCmd->GetNewIntValue(value));
}
