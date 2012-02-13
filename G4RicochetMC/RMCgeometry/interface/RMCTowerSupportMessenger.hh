#ifndef RMCTowerSupportMessenger_h
#define RMCTowerSupportMessenger_h 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerSupportMessenger.hh                         //
//                                                                    //
//  Description: Messenger class to allow setting RMC tower support  //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        8 December 2010                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class RMCTowerSupportConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;


class RMCTowerSupportMessenger: public G4UImessenger {
public:
  RMCTowerSupportMessenger(RMCTowerSupportConstruction* builder);
  virtual ~RMCTowerSupportMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  RMCTowerSupportConstruction* theBuilder;
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

#endif	/* RMCTowerSupportMessenger_hh */
