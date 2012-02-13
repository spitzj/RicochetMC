#ifndef CDMSVesselMessenger_h
#define CDMSVesselMessenger_h 1
// $Id: CDMSVesselMessenger.hh,v 1.4 2010/12/23 23:41:15 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVesselMessenger.hh                               //
//                                                                    //
//  Description: Messenger class to allow setting CDMS cryostat       //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
//  20101208  M. Kelsey -- Move tower support structure to new class  //
//  20101215  M. Kelsey -- Add number of vacuum vessels               //
//  20101223  M. Kelsey -- Use new vector-input command (CDMS only)   //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class CDMSVesselConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class CDMS_UIcmdDoublesListAndUnit;


class CDMSVesselMessenger: public G4UImessenger {
public:
  CDMSVesselMessenger(CDMSVesselConstruction* builder);
  virtual ~CDMSVesselMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  CDMSVesselConstruction* theBuilder;
  G4UIdirectory* cmdDir;

  G4UIcmdWithAnInteger* NStagesCmd;
  G4UIcmdWithAnInteger* NVacuumCmd;

  CDMS_UIcmdDoublesListAndUnit* VesselThickCmd;
  G4UIcmdWithADoubleAndUnit* Vessel0RadCmd;
  G4UIcmdWithADoubleAndUnit* VesselDeltaRadCmd;
  G4UIcmdWithADoubleAndUnit* VesselDeltaHeightCmd;
  G4UIcmdWithADoubleAndUnit* VesselExtraHeightCmd;
  G4UIcmdWithADoubleAndUnit* VesselGapCmd;

  CDMS_UIcmdDoublesListAndUnit* LidThicknessCmd;
  G4UIcmdWithADoubleAndUnit* StemRadCmd;
  G4UIcmdWithADoubleAndUnit* VacRadCmd;
  G4UIcmdWithADoubleAndUnit* PipeThickCmd;
  G4UIcmdWithADoubleAndUnit* PipeBaseLenCmd;
};

#endif	/* CDMSVesselMessenger_hh */
