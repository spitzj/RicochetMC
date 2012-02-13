#ifndef CDMSZipMessenger_h
#define CDMSZipMessenger_h 1
// $Id: CDMSZipMessenger.hh,v 1.3 2010/12/09 20:54:40 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSZipMessenger.hh                                  //
//                                                                    //
//  Description: Messenger class to allow setting CDMS single-ZIP     //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Michael Kesey (SLAC)                                 //
//  Date:        24 November 2010                                     //
//                                                                    //
//  20101203  M. Kelsey -- Add flag to suppress housing construction  //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class CDMSZipConstruction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;


class CDMSZipMessenger: public G4UImessenger {
public:
  CDMSZipMessenger(CDMSZipConstruction* builder);
  virtual ~CDMSZipMessenger();
  
  void SetNewValue(G4UIcommand* cmd, G4String value);

private:
  CDMSZipConstruction* theBuilder;
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

#endif	/* CDMSZipMessenger_hh */
