// $Id: CDMSTowerPatternMessenger.cc,v 1.3 2011/02/05 05:54:26 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSTowerPatternMessenger.cc                         //
//                                                                    //
//  Description: Messenger class to allow setting CDMS multi-tower    //
//               layout parameters                                    //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
//  20110121  M. Kelsey -- Use symbolic names for layout codes        //
//  20110204  M. Kelsey -- Add new command for Optimize(mass)         //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/interface/CDMSTowerPatternMessenger.hh"
#include "CDMSgeometry/detectors/CDMSTowerPattern.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

CDMSTowerPatternMessenger::CDMSTowerPatternMessenger(CDMSTowerPattern* pattern)
  : towerPattern(pattern) {
  cmdDir = new G4UIdirectory("/CDMS/Layout/");

  LayoutCodeCmd = new G4UIcmdWithAString("/CDMS/Layout/Code", this);
  LayoutCodeCmd->SetGuidance("Positons and rotation of towers in cryostat");
  LayoutCodeCmd->SetGuidance("ring    : circular pattern with central tower");
  LayoutCodeCmd->SetGuidance("        : tower sides face center");
  LayoutCodeCmd->SetGuidance("ringRot : same as 0, tower corners face center");
  LayoutCodeCmd->SetGuidance("grid    : diamond pattern of zips");
  LayoutCodeCmd->SetGuidance("hex     : close packed hexagonal array");
  LayoutCodeCmd->SetParameterName("LayoutCode",true,true);
  LayoutCodeCmd->SetCandidates("ring ringRot grid hex");
  LayoutCodeCmd->SetDefaultValue(towerPattern->GetLayoutName());

  NumberOfTowersCmd = new G4UIcmdWithAnInteger("/CDMS/Layout/Towers", this);
  NumberOfTowersCmd->SetGuidance("Number of towers in full detector");
  NumberOfTowersCmd->SetParameterName("NumberOfTowers",true,true);
  NumberOfTowersCmd->SetDefaultValue(towerPattern->GetNumberOfTowers());

  OptimizeCmd = new G4UIcmdWithADoubleAndUnit("/CDMS/Layout/Optimize", this);
  OptimizeCmd->SetGuidance("Generate tower pattern for detector");
  OptimizeCmd->SetGuidance("with specified active mass (default is 100 kg)");
  OptimizeCmd->SetParameterName("Mass",true,true);
  OptimizeCmd->SetUnitCategory("Mass");
  OptimizeCmd->SetDefaultValue(100.*kg);
}

CDMSTowerPatternMessenger::~CDMSTowerPatternMessenger() {
  delete LayoutCodeCmd;
  delete NumberOfTowersCmd;
  delete OptimizeCmd;
  delete cmdDir;
}


// Propagate user input to geometry model

void CDMSTowerPatternMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (cmd == LayoutCodeCmd) towerPattern->SetLayout(value);

  if (cmd == NumberOfTowersCmd)
    towerPattern->SetNumberOfTowers(NumberOfTowersCmd->GetNewIntValue(value));

  if (cmd == OptimizeCmd)
    towerPattern->Optimize(OptimizeCmd->GetNewDoubleValue(value));
}

