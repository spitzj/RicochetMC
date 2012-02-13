#ifndef RMCVesselMessenger_h
#define RMCVesselMessenger_h 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVesselMessenger.hh                                //
//                                                                    //
//  Description: Messenger class to allow setting RMC cryostat        //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kesey (SLAC)                   //
//  Date:        14 January 2012                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class RMCVesselConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class RMC_UIcmdDoublesListAndUnit;


class RMCVesselMessenger: public G4UImessenger {
public:
  RMCVesselMessenger(RMCVesselConstruction* builder);
  virtual ~RMCVesselMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  RMCVesselConstruction* theBuilder;
  G4UIdirectory* cmdDir;

  G4UIcmdWithAnInteger* NStagesCmd;
  G4UIcmdWithAnInteger* NVacuumCmd;

  G4UIcmdWithADoubleAndUnit* Vessel0RadCmd;
  G4UIcmdWithADoubleAndUnit* VesselDeltaRadCmd;
  G4UIcmdWithADoubleAndUnit* VesselDeltaHeightCmd;
  G4UIcmdWithADoubleAndUnit* VesselExtraHeightCmd;
  G4UIcmdWithADoubleAndUnit* VesselGapCmd;

  G4UIcmdWithADoubleAndUnit* StemRadCmd;
  G4UIcmdWithADoubleAndUnit* VacRadCmd;
  G4UIcmdWithADoubleAndUnit* PipeThickCmd;
  G4UIcmdWithADoubleAndUnit* PipeBaseLenCmd;
};

#endif	/* RMCVesselMessenger_hh */
