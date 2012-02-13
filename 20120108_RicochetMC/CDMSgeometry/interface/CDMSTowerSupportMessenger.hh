#ifndef CDMSTowerSupportMessenger_h
#define CDMSTowerSupportMessenger_h 1
// $Id: CDMSTowerSupportMessenger.hh,v 1.2 2010/12/15 08:15:15 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSTowerSupportMessenger.hh                         //
//                                                                    //
//  Description: Messenger class to allow setting CDMS tower support  //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        8 December 2010                                      //
//                                                                    //
//  20101214  M. Kelsey -- Add parameters for special stage indices   //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class CDMSTowerSupportConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class CDMSTowerSupportMessenger: public G4UImessenger {
public:
  CDMSTowerSupportMessenger(CDMSTowerSupportConstruction* builder);
  virtual ~CDMSTowerSupportMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  CDMSTowerSupportConstruction* theBuilder;
  G4UIdirectory* cmdDir;

  G4UIcmdWithAnInteger*      NStagesCmd;
  G4UIcmdWithADoubleAndUnit* StageHeightCmd;
  G4UIcmdWithADoubleAndUnit* StageWallThickCmd;
  G4UIcmdWithADoubleAndUnit* StageGapCmd;
  G4UIcmdWithAnInteger*      VesselMountStageCmd;
  
  G4UIcmdWithADoubleAndUnit* SpoolLengthCmd; 
  G4UIcmdWithADoubleAndUnit* SpoolIRCmd;
  G4UIcmdWithADoubleAndUnit* SpoolThicknessCmd;
  G4UIcmdWithADoubleAndUnit* CTubeIRCmd;
  G4UIcmdWithADoubleAndUnit* CTubeThickCmd;

  G4UIcmdWithADoubleAndUnit* SquidLenCmd;
  G4UIcmdWithADoubleAndUnit* SquidThickCmd;
  G4UIcmdWithADoubleAndUnit* SquidWidthCmd;
  G4UIcmdWithAnInteger*      SquidStageCmd;
  G4UIcmdWithADoubleAndUnit* FETLenCmd;
  G4UIcmdWithAnInteger*      FETStageCmd;
};

#endif	/* CDMSTowerSupportMessenger_hh */
