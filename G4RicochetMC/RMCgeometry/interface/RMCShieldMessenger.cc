////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCShieldMessenger.cc                               //
//                                                                    //
//  Description: Messenger class to allow setting of RMC shielding   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        22 December 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCShieldMessenger.hh"
#include "RMCgeometry/detectors/RMCShieldConstruction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

RMCShieldMessenger::RMCShieldMessenger(RMCShieldConstruction* builder)
  : RMCMessengerBase("/RMC/Shield/",
		      "Configuration of RMC shielding and veto structure"),
    theBuilder(builder) {

  ScintThickCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("ScintThick",
	"Thickness of all scintillator panels");
  ScintThickCmd->SetParameterName("ScintThick", true, true);
  ScintThickCmd->SetUnitCategory("Length");
  ScintThickCmd->SetDefaultValue(theBuilder->GetScintThick());

  TopPanelsXCmd = CreateCommand<G4UIcmdWithAnInteger>("TopPanelsX",
	"Number of end panels along X axis");
  TopPanelsXCmd->SetParameterName("N", false, true);
  TopPanelsXCmd->SetRange("N>=0");
  TopPanelsXCmd->SetDefaultValue(theBuilder->GetNTopPanelsX());

  TopPanelsYCmd = CreateCommand<G4UIcmdWithAnInteger>("TopPanelsY",
	"Number of end panels along Y axis");
  TopPanelsYCmd->SetParameterName("N", false, true);
  TopPanelsYCmd->SetRange("N>=0");
  TopPanelsYCmd->SetDefaultValue(theBuilder->GetNTopPanelsY());

  OverhangCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("Overhang",
	"Extension of end panels beyond outer shield radius");
  OverhangCmd->SetParameterName("Overhang", false, true);
  OverhangCmd->SetUnitCategory("Length");
  OverhangCmd->SetDefaultValue(theBuilder->GetOverhang());

  OverlapXCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("OverlapX",
	"Overlap of end panels along X axis");
  OverlapXCmd->SetParameterName("OverlapX", false, true);
  OverlapXCmd->SetUnitCategory("Length");
  OverlapXCmd->SetDefaultValue(theBuilder->GetOverlapX());

  OverlapYCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("OverlapY",
	"Overlap of end panels along Y axis");
  OverlapYCmd->SetParameterName("OverlapY", false, true);
  OverlapYCmd->SetUnitCategory("Length");
  OverlapYCmd->SetDefaultValue(theBuilder->GetOverlapY());

  TopVetoSHCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("TopVetoSupportHeight",
	"Gap between top of shielding and end panels");
  TopVetoSHCmd->SetParameterName("TopVetoSH", false, true);
  TopVetoSHCmd->SetUnitCategory("Length");
  TopVetoSHCmd->SetDefaultValue(theBuilder->GetTopVetoSupportHeight());

  BottomVetoSHCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("BottomVetoSupportHeight",
	"Gap between top of shielding and end panels");
  BottomVetoSHCmd->SetParameterName("BottomVetoSH", false, true);
  BottomVetoSHCmd->SetUnitCategory("Length");
  BottomVetoSHCmd->SetDefaultValue(theBuilder->GetBottomVetoSupportHeight());

  SidePanelsCmd = CreateCommand<G4UIcmdWithAnInteger>("SidePanels",
	"Number of azimuthal panel faces surrounding shielding");
  SidePanelsCmd->SetParameterName("N", false, true);
  SidePanelsCmd->SetRange("N==0||N>=3");		// Must be true polygon
  SidePanelsCmd->SetDefaultValue(theBuilder->GetNSidePanels());

  SideClearanceCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("SideClearance",
	"Gap between side of shielding and azimuthal panels");
  SideClearanceCmd->SetParameterName("SideClearance", false, true);
  SideClearanceCmd->SetUnitCategory("Length");
  SideClearanceCmd->SetDefaultValue(theBuilder->GetSidePanelClearance());

  SideOverlapCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("SideOverlap",
	"Vertical overlap of azimuthal panels");
  SideOverlapCmd->SetParameterName("SideOverlap", false, true);
  SideOverlapCmd->SetUnitCategory("Length");
  SideOverlapCmd->SetDefaultValue(theBuilder->GetSidePanelOverlap());

  SideCornerCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("SideCorner",
	"Coverage of corner extensions of azimuthal panels");
  SideCornerCmd->SetParameterName("SideCorner", false, true);
  SideCornerCmd->SetUnitCategory("Length");
  SideCornerCmd->SetDefaultValue(theBuilder->GetSideCornerOverlap());
}

RMCShieldMessenger::~RMCShieldMessenger() {
  delete ScintThickCmd;
  delete TopPanelsXCmd;
  delete TopPanelsYCmd;
  delete OverhangCmd;
  delete OverlapXCmd;
  delete OverlapYCmd;
  delete TopVetoSHCmd;
  delete BottomVetoSHCmd;
  delete SidePanelsCmd;
  delete SideClearanceCmd;
  delete SideOverlapCmd;
  delete SideCornerCmd;
}


// Propagate user input to geometry model

void RMCShieldMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  RMCMessengerBase::SetNewValue(cmd, value);	// To process global commands

  if (cmd == ScintThickCmd)
    theBuilder->SetScintThick(ScintThickCmd->GetNewDoubleValue(value));

  if (cmd == TopPanelsXCmd)
    theBuilder->SetNTopPanelsX(TopPanelsXCmd->GetNewIntValue(value));
  if (cmd == TopPanelsYCmd)
    theBuilder->SetNTopPanelsY(TopPanelsYCmd->GetNewIntValue(value));
  if (cmd == OverhangCmd)
    theBuilder->SetOverhang(OverhangCmd->GetNewDoubleValue(value));
  if (cmd == OverlapXCmd)
    theBuilder->SetOverlapX(OverlapXCmd->GetNewDoubleValue(value));
  if (cmd == OverlapYCmd)
    theBuilder->SetOverlapY(OverlapYCmd->GetNewDoubleValue(value));
  if (cmd == TopVetoSHCmd)
    theBuilder->SetTopVetoSupportHeight(TopVetoSHCmd->GetNewDoubleValue(value));
  if (cmd == BottomVetoSHCmd)
    theBuilder->SetBottomVetoSupportHeight(BottomVetoSHCmd->GetNewDoubleValue(value));

  if (cmd == SidePanelsCmd)
    theBuilder->SetNSidePanels(SidePanelsCmd->GetNewIntValue(value));
  if (cmd == SideClearanceCmd)
    theBuilder->SetSidePanelClearance(SideClearanceCmd->GetNewDoubleValue(value));
  if (cmd == SideOverlapCmd)
    theBuilder->SetSidePanelOverlap(SideOverlapCmd->GetNewDoubleValue(value));
  if (cmd == SideCornerCmd)
    theBuilder->SetSideCornerOverlap(SideCornerCmd->GetNewDoubleValue(value));
}

