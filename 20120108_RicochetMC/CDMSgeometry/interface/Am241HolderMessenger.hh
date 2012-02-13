#ifndef Am241HolderMessenger_hh
#define Am241HolderMessenger_hh 1
////////////////////////////////////////////////////////////////////////
// $Id: Am241HolderMessenger.hh,v 1.3 2011/06/24 03:37:31 kelsey Exp $
//  File:        Am241HolderMessenger.hh                              //
//  Description: Configuration commands for Am241SourceHolder         //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        26 April 2011                                        //
//                                                                    //
//  20110427  M. Kelsey -- Add command to use hardwired gamma lines   //
//  20110623  M. Kelsey -- Add configurable foil material             //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSMessengerBase.hh"
#include "globals.hh"

class Am241SourceHolder;
class CDMS_UIcmdDoublesListAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;


class Am241HolderMessenger: public CDMSMessengerBase {
public:
  Am241HolderMessenger(Am241SourceHolder* holder);
  virtual ~Am241HolderMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  void MakeLocatorCommands();
  void MakePlateCommands();
  void MakeCanisterCommands();
  void MakeSourceCommands();
  void MakeShieldCommands();

  Am241SourceHolder* theHolder;

  G4UIcmdWith3VectorAndUnit*    SrcPosCmd;	// Source locator
  G4UIcmdWith3Vector*           SrcDirCmd;
  G4UIcmdWithABool*             UseZipCmd;	// Reference ZIP geometry
  G4UIcmdWithABool*             UseGammasCmd;	// Hardwired CDMS gamma lines
  G4UIcmdWithABool*             DrawSolidCmd;
  G4UIcmdWithAString*           MaterialCmd;
  G4UIcmdWithAnInteger*         SidesCmd;
  G4UIcmdWithADoubleAndUnit*    PlateRCmd;
  G4UIcmdWithADoubleAndUnit*    PlateLCmd;
  G4UIcmdWithADoubleAndUnit*    CanHeightCmd;
  G4UIcmdWithADoubleAndUnit*    CanThickCmd;
  G4UIcmdWithADoubleAndUnit*    PuckHeightCmd;
  G4UIcmdWithADoubleAndUnit*    LeadThicknessCmd;
  G4UIcmdWithADoubleAndUnit*    FoilThicknessCmd;
  G4UIcmdWithAString*           FoilMaterialCmd;
  G4UIcmdWithADoubleAndUnit*    LeadHoleRCmd;
  G4UIcmdWithADoubleAndUnit*    PlateHoleRCmd;
  G4UIcmdWithAnInteger*         NumberOfCansCmd;
  CDMS_UIcmdDoublesListAndUnit* ActivityCmd;
  CDMS_UIcmdDoublesListAndUnit* PuckRCmd;
  CDMS_UIcmdDoublesListAndUnit* ActiveRCmd;
  CDMS_UIcmdDoublesListAndUnit* ActiveLCmd;
  CDMS_UIcmdDoublesListAndUnit* CanPosRCmd;
  CDMS_UIcmdDoublesListAndUnit* CanPosPhiCmd;
  CDMS_UIcmdDoublesListAndUnit* HolePosRCmd;
  CDMS_UIcmdDoublesListAndUnit* HolePosPhiCmd;
};

#endif	/* Am241HolderMessenger_hh */
