#ifndef RMCZipMessenger_h
#define RMCZipMessenger_h 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCZipMessenger.hh                                  //
//                                                                    //
//  Description: Messenger class to allow setting RMC single-ZIP     //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class RMCZipConstruction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;


class RMCZipMessenger: public G4UImessenger {
public:
  RMCZipMessenger(RMCZipConstruction* builder);
  virtual ~RMCZipMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  RMCZipConstruction* theBuilder;
  G4UIdirectory* cmdDir;
  G4UIcmdWithABool*          MakeHousingCmd;
  G4UIcmdWithAString*        ZipMaterialCmd;
  G4UIcmdWithAString*        HousingMaterialCmd;
  G4UIcmdWithADoubleAndUnit* ZipRadCmd;
  G4UIcmdWithADoubleAndUnit* ZipThickCmd;
  G4UIcmdWithADoubleAndUnit* ZipAxis1LenCmd;
  G4UIcmdWithADoubleAndUnit* ZipAxis2LenCmd;
  G4UIcmdWithAnInteger*      HousingSidesCmd;
  G4UIcmdWithADoubleAndUnit* HousingThicknessCmd;
  G4UIcmdWithADoubleAndUnit* ZipClearanceRCmd;
  G4UIcmdWithADoubleAndUnit* ZipClearanceZCmd;
  G4UIcmdWith3VectorAndUnit* ZipPositionCmd;
};

#endif	/* RMCZipMessenger_hh */
