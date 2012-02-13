// $Id: CDMSShieldMessenger.cc,v 1.8 2011/05/25 05:10:27 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSShieldMessenger.cc                               //
//                                                                    //
//  Description: Messenger class to allow setting of CDMS shielding   //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        22 December 2010                                     //
//                                                                    //
//  20110112  M. Kelsey -- Use new CDMS_UIcmdShieldLayer, base class, //
//		populate all necessary commands including guidance.   //
//  20110524  M. Kelsey -- Forgot to call base SetNewValue(), use new //
//		CreateCommand interface.                              //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/interface/CDMSShieldMessenger.hh"
#include "CDMSgeometry/detectors/CDMSShieldConstruction.hh"
#include "CDMSgeometry/interface/CDMS_UIcmdShieldLayer.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


// Constructor and destructor

CDMSShieldMessenger::CDMSShieldMessenger(CDMSShieldConstruction* builder)
  : CDMSMessengerBase("/CDMS/Shield/",
		      "Configuration of CDMS shielding and veto structure"),
    theBuilder(builder) {
  MuMetalCmd = CreateCommand<CDMS_UIcmdShieldLayer>("MuMetal", "Magnetic shield outside cryostat");
  MuMetalCmd->SetParameterName("MuMetal", false, true);
  MuMetalCmd->SetDefaultValue(theBuilder->GetMuMetalParams());

  InnerPolyCmd = CreateCommand<CDMS_UIcmdShieldLayer>("InnerPoly", "Polyethylene (neutron) shield");
  InnerPolyCmd->SetParameterName("InnerPoly", false, true);
  InnerPolyCmd->SetDefaultValue(theBuilder->GetInnerPolyParams());

  InnerLeadCmd = CreateCommand<CDMS_UIcmdShieldLayer>("InnerLead",
	"Lead (cosmic ray) shield");
  InnerLeadCmd->SetParameterName("InnerLead", false, true);
  InnerLeadCmd->SetDefaultValue(theBuilder->GetInnerLeadParams());

  OuterPolyCmd = CreateCommand<CDMS_UIcmdShieldLayer>("OuterPoly",
	"Polyethylene (neutron) shield");
  OuterPolyCmd->SetParameterName("OuterPoly", false, true);
  OuterPolyCmd->SetDefaultValue(theBuilder->GetOuterPolyParams());

  OuterLeadCmd = CreateCommand<CDMS_UIcmdShieldLayer>("OuterLead",
	"Lead (cosmic ray) shield");
  OuterLeadCmd->SetParameterName("OuterLead", false, true);
  OuterLeadCmd->SetDefaultValue(theBuilder->GetOuterLeadParams());

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

CDMSShieldMessenger::~CDMSShieldMessenger() {
  delete MuMetalCmd;
  delete InnerPolyCmd;
  delete InnerLeadCmd;
  delete OuterPolyCmd;
  delete OuterLeadCmd;
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

void CDMSShieldMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  CDMSMessengerBase::SetNewValue(cmd, value);	// To process global commands

  if (cmd == MuMetalCmd)
    theBuilder->SetMuMetalParams(MuMetalCmd->GetNewLayerData(value));
  if (cmd == InnerPolyCmd)
    theBuilder->SetInnerPolyParams(InnerPolyCmd->GetNewLayerData(value));
  if (cmd == InnerLeadCmd)
    theBuilder->SetInnerPolyParams(InnerLeadCmd->GetNewLayerData(value));
  if (cmd == OuterLeadCmd)
    theBuilder->SetInnerPolyParams(OuterLeadCmd->GetNewLayerData(value));
  if (cmd == OuterPolyCmd)
    theBuilder->SetInnerPolyParams(OuterPolyCmd->GetNewLayerData(value));

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

