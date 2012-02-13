////////////////////////////////////////////////////////////////////////
// $Id: Am241HolderMessenger.cc,v 1.5 2011/06/24 03:37:31 kelsey Exp $
//  File:        Am241HolderMessenger.hh                              //
//  Description: Configuration commands for Am241SourceHolder         //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        26 April 2011                                        //
//                                                                    //
//  20110427  M. Kelsey -- Add command to use hardwired gamma lines   //
//  20110524  M. Kelsey -- Forgot to call base SetNewValue(), use new //
//		CreateCommand interface.                              //
//  20110623  M. Kelsey -- Add configurable foil material             //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/interface/Am241HolderMessenger.hh"
#include "CDMSgeometry/detectors/Am241SourceHolder.hh"
#include "CDMSg4base/CDMSMessengerBase.hh"
#include "CDMSg4base/CDMS_UIcmdDoublesListAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"


// Constructor and destructor

Am241HolderMessenger::Am241HolderMessenger(Am241SourceHolder* holder)
  : CDMSMessengerBase("/CDMS/Am241/", "Configure Am-241 source holder"),
    theHolder(holder) {
  MakeLocatorCommands();
  MakePlateCommands();
  MakeCanisterCommands();
  MakeSourceCommands();
  MakeShieldCommands();
}


void Am241HolderMessenger::MakeLocatorCommands() {
  SrcPosCmd = CreateCommand<G4UIcmdWith3VectorAndUnit>("Position",
	"Set source position");
  SrcPosCmd->SetParameterName("X", "Y", "Z", false);
  SrcPosCmd->SetUnitCategory("Length");  
  SrcPosCmd->SetDefaultValue(theHolder->GetPosition());
  SrcPosCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SrcDirCmd = CreateCommand<G4UIcmdWith3Vector>("Direction",
	"Set source direction");
  SrcDirCmd->SetParameterName("VX", "VY", "VZ", false);
  SrcDirCmd->SetDefaultValue(theHolder->GetDirection());
  SrcDirCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

void Am241HolderMessenger::MakePlateCommands() {
  UseZipCmd = CreateCommand<G4UIcmdWithABool>("UseZip",
	"Get plate dimensions from ZIP geometry?");
  UseZipCmd->SetParameterName("UseZip",true,false);
  UseZipCmd->SetDefaultValue(true);
  UseZipCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  DrawSolidCmd = CreateCommand<G4UIcmdWithABool>("DrawSolid",
	"Draw source holder opaque (true) or wire-frame");
  DrawSolidCmd->SetParameterName("DrawSolid",true,false);
  DrawSolidCmd->SetDefaultValue(true);
  DrawSolidCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  MaterialCmd = CreateCommand<G4UIcmdWithAString>("Material",
	"Material used for holder plate and canister");
  MaterialCmd->SetParameterName("Mat", false);
  MaterialCmd->SetDefaultValue(theHolder->GetPlateMaterial());
  MaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SidesCmd = CreateCommand<G4UIcmdWithAnInteger>("Sides",
	"Number of sides on holder plate (same as ZIP housing)");
  SidesCmd->SetParameterName("NSides", false);
  SidesCmd->SetDefaultValue(theHolder->GetPlateSides());
  SidesCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PlateRCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("Radius",
	"Radius (center to side) of holder plate");
  PlateRCmd->SetParameterName("PlateR", false);
  PlateRCmd->SetUnitCategory("Length");
  PlateRCmd->SetDefaultValue(theHolder->GetPlateRadius());
  PlateRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PlateLCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("Thickness",
	"Thickness of holder plate");
  PlateLCmd->SetParameterName("PlateL", false);
  PlateLCmd->SetUnitCategory("Length");
  PlateLCmd->SetDefaultValue(theHolder->GetPlateThickness());
  PlateLCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PlateHoleRCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("PlateHoleRadius",
	"Radius of collimator holes in holder plate");
  PlateHoleRCmd->SetParameterName("PlateHoleR", false);
  PlateHoleRCmd->SetUnitCategory("Length");
  PlateHoleRCmd->SetDefaultValue(theHolder->GetPlateHoleRadius());
  PlateHoleRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

void Am241HolderMessenger::MakeCanisterCommands() {
  NumberOfCansCmd = CreateCommand<G4UIcmdWithAnInteger>("NumberOfCans",
	"Number of mounted source canisters");
  NumberOfCansCmd->SetParameterName("NCans", false);
  NumberOfCansCmd->SetDefaultValue(theHolder->GetNumberOfCans());
  NumberOfCansCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CanHeightCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("CanHeight",
	"Height of all source canisters");
  CanHeightCmd->SetParameterName("CanH", false);
  CanHeightCmd->SetUnitCategory("Length");
  CanHeightCmd->SetDefaultValue(theHolder->GetCanHeight());
  CanHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CanThickCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("CanThickness",
	"Thickness of source canister material");
  CanThickCmd->SetParameterName("CanThick", false);
  CanThickCmd->SetUnitCategory("Length");
  CanThickCmd->SetDefaultValue(theHolder->GetCanThickness());
  CanThickCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CanPosRCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("CanPositionR",
	"Radial position of canister on holder plate");
  CanPosRCmd->SetParameterName("CanPosR", false);
  CanPosRCmd->SetUnitCategory("Length");
  CanPosRCmd->SetDefaultValue(theHolder->GetCanPositionR());
  CanPosRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  CanPosPhiCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("CanPositionPhi",
	"Azimuthal position of canister on holder plate");
  CanPosPhiCmd->SetParameterName("CanPosPhi", false);
  CanPosPhiCmd->SetUnitCategory("Angle");
  CanPosPhiCmd->SetDefaultValue(theHolder->GetCanPositionPhi());
  CanPosPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

void Am241HolderMessenger::MakeSourceCommands() {
  UseGammasCmd = CreateCommand<G4UIcmdWithABool>("UseCDMSGammas",
	"Used hardwired gamma lines from CDMS?");
  UseGammasCmd->SetParameterName("UseGammas",true,false);
  UseGammasCmd->SetDefaultValue(true);
  UseGammasCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ActivityCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("Activity",
	"Activity of each source");
  ActivityCmd->SetParameterName("Strength", false);
  ActivityCmd->SetUnitCategory("Activity");
  ActivityCmd->SetDefaultValue(theHolder->GetSourceActivity());
  ActivityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PuckHeightCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("PuckHeight",
	"Height of all source disks");
  PuckHeightCmd->SetParameterName("PuckH", false);
  PuckHeightCmd->SetUnitCategory("Length");
  PuckHeightCmd->SetDefaultValue(theHolder->GetSourcePuckHeight());
  PuckHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PuckRCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("PuckRadius",
	"Radius of source disk (hockey puck)");
  PuckRCmd->SetParameterName("PuckR", false);
  PuckRCmd->SetUnitCategory("Length");
  PuckRCmd->SetDefaultValue(theHolder->GetSourcePuckRadius());
  PuckRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ActiveRCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("ActiveRadius",
	"Radius of active region within source disk");
  ActiveRCmd->SetParameterName("ActiveR", false);
  ActiveRCmd->SetUnitCategory("Length");
  ActiveRCmd->SetDefaultValue(theHolder->GetActiveRadius());
  ActiveRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  ActiveLCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("ActiveHeight",
	"Height of active region within source disk");
  ActiveLCmd->SetParameterName("ActiveL", false);
  ActiveLCmd->SetUnitCategory("Length");
  ActiveLCmd->SetDefaultValue(theHolder->GetActiveHeight());
  ActiveLCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

void Am241HolderMessenger::MakeShieldCommands() {
  LeadThicknessCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("LeadThickness",
	"Thickness of lead shields");
  LeadThicknessCmd->SetParameterName("LeadThick", false);
  LeadThicknessCmd->SetUnitCategory("Length");
  LeadThicknessCmd->SetDefaultValue(theHolder->GetLeadDiskThickness());
  LeadThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FoilThicknessCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("FoilThickness",
	"Thickness of alpha-filtering foils");
  FoilThicknessCmd->SetParameterName("FoilThick", false);
  FoilThicknessCmd->SetUnitCategory("Length");
  FoilThicknessCmd->SetDefaultValue(theHolder->GetAlphaFoilThickness());
  FoilThicknessCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  FoilMaterialCmd = CreateCommand<G4UIcmdWithAString>("FoilMaterial", "Material used for alpha foils");
  FoilMaterialCmd->SetParameterName("Mat", false);
  FoilMaterialCmd->SetDefaultValue(theHolder->GetAlphaFoilMaterial());
  FoilMaterialCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  FoilMaterialCmd->SetCandidates("G4_Al G4_Cu G4_KAPTON");

  LeadHoleRCmd = CreateCommand<G4UIcmdWithADoubleAndUnit>("LeadHoleRadius",
	"Radius of collimator holes in lead");
  LeadHoleRCmd->SetParameterName("LeadHoleR", false);
  LeadHoleRCmd->SetUnitCategory("Length");
  LeadHoleRCmd->SetDefaultValue(theHolder->GetLeadHoleRadius());
  LeadHoleRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HolePosRCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("LeadHoleR",
	"Radial position of collimator hole in lead shield");
  HolePosRCmd->SetParameterName("HolePosR", false);
  HolePosRCmd->SetUnitCategory("Length");
  HolePosRCmd->SetDefaultValue(theHolder->GetHolePositionR());
  HolePosRCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  HolePosPhiCmd = CreateCommand<CDMS_UIcmdDoublesListAndUnit>("LeadHolePhi",
	"Azimuthal position of collimator hole in lead shield");
  HolePosPhiCmd->SetParameterName("HolePosPhi", false);
  HolePosPhiCmd->SetUnitCategory("Angle");
  HolePosPhiCmd->SetDefaultValue(theHolder->GetHolePositionPhi());
  HolePosPhiCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

Am241HolderMessenger::~Am241HolderMessenger() {
  delete SrcPosCmd;
  delete SrcDirCmd;
  delete UseZipCmd;
  delete UseGammasCmd;
  delete DrawSolidCmd;
  delete MaterialCmd;
  delete SidesCmd;
  delete PlateRCmd;
  delete PlateLCmd;
  delete CanHeightCmd;
  delete CanThickCmd;
  delete PuckHeightCmd;
  delete LeadThicknessCmd;
  delete FoilThicknessCmd;
  delete FoilMaterialCmd;
  delete LeadHoleRCmd;
  delete PlateHoleRCmd;
  delete NumberOfCansCmd;
  delete ActivityCmd;
  delete PuckRCmd;
  delete ActiveRCmd;
  delete ActiveLCmd;
  delete CanPosRCmd;
  delete CanPosPhiCmd;
  delete HolePosRCmd;
  delete HolePosPhiCmd;
}


// Process commands

void Am241HolderMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  CDMSMessengerBase::SetNewValue(cmd, value);	// To process global commands

  if (cmd == SrcPosCmd) theHolder->SetPosition(SrcPosCmd->GetNew3VectorValue(value) );
  if (cmd == SrcDirCmd) theHolder->SetDirection(SrcDirCmd->GetNew3VectorValue(value) );

  if (cmd == UseZipCmd) theHolder->UseZipBuilder(UseZipCmd->GetNewBoolValue(value) );
  if (cmd == UseGammasCmd) theHolder->UseCDMSGammas(UseZipCmd->GetNewBoolValue(value) );
  if (cmd == DrawSolidCmd) theHolder->DrawSolid(DrawSolidCmd->GetNewBoolValue(value) );
  if (cmd == MaterialCmd) theHolder->SetPlateMaterial(value);
  if (cmd == SidesCmd) theHolder->SetPlateSides(SidesCmd->GetNewIntValue(value) );
  if (cmd == PlateRCmd) theHolder->SetPlateRadius(PlateRCmd->GetNewDoubleValue(value) );
  if (cmd == PlateLCmd) theHolder->SetPlateThickness(PlateLCmd->GetNewDoubleValue(value) );
  if (cmd == CanHeightCmd) theHolder->SetCanHeight(CanHeightCmd->GetNewDoubleValue(value) );
  if (cmd == CanThickCmd) theHolder->SetCanThickness(CanThickCmd->GetNewDoubleValue(value) );
  if (cmd == PuckHeightCmd) theHolder->SetSourcePuckHeight(PuckHeightCmd->GetNewDoubleValue(value) );
  if (cmd == LeadThicknessCmd) theHolder->SetLeadDiskThickness(LeadThicknessCmd->GetNewDoubleValue(value) );
  if (cmd == FoilThicknessCmd) theHolder->SetAlphaFoilThickness(FoilThicknessCmd->GetNewDoubleValue(value) );
  if (cmd == FoilMaterialCmd) theHolder->SetAlphaFoilMaterial(value);
  if (cmd == LeadHoleRCmd) theHolder->SetLeadHoleRadius(LeadHoleRCmd->GetNewDoubleValue(value) );
  if (cmd == PlateHoleRCmd) theHolder->SetPlateHoleRadius(PlateHoleRCmd->GetNewDoubleValue(value) );
  if (cmd == NumberOfCansCmd) theHolder->SetNumberOfCans(NumberOfCansCmd->GetNewIntValue(value) );
  if (cmd == ActivityCmd) theHolder->SetSourceActivity(ActivityCmd->GetNewListValue(value) );
  if (cmd == PuckRCmd) theHolder->SetSourcePuckRadius(PuckRCmd->GetNewListValue(value) );
  if (cmd == ActiveRCmd) theHolder->SetActiveRadius(ActiveRCmd->GetNewListValue(value) );
  if (cmd == ActiveLCmd) theHolder->SetActiveHeight(ActiveLCmd->GetNewListValue(value) );
  if (cmd == CanPosRCmd) theHolder->SetCanPositionR(CanPosRCmd->GetNewListValue(value) );
  if (cmd == CanPosPhiCmd) theHolder->SetCanPositionPhi(CanPosPhiCmd->GetNewListValue(value) );
  if (cmd == HolePosRCmd) theHolder->SetHolePositionR(HolePosRCmd->GetNewListValue(value) );
  if (cmd == HolePosPhiCmd) theHolder->SetHolePositionPhi(HolePosPhiCmd->GetNewListValue(value) );
}
