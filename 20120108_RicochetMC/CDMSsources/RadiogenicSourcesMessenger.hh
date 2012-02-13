#ifndef RadiogenicSourcesMessenger_hh
#define RadiogenicSourcesMessenger_hh
// $Id: RadiogenicSourcesMessenger.hh,v 1.1 2010/12/14 07:41:00 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RadiogenicSourcesMessenger.hh                        //
//  Description: User interface for all radioactive decays            //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        13 December 2010                                     //
//                                                                    //
//  20101213  M. Kelsey -- Copy CDMSmain/CDMS_ParticleSourceMessenger //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class RadiogenicSources;
class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;


class RadiogenicSourcesMessenger : public G4UImessenger {
public:
  RadiogenicSourcesMessenger(RadiogenicSources *fPtclGun);
  ~RadiogenicSourcesMessenger();
  
  void SetNewValue(G4UIcommand *command, G4String newValues);
  
private:
  RadiogenicSources          *fParticleGun;
  G4UIdirectory              *gunDirectory;

  G4UIcmdWithAString         *typeCmd;
  G4UIcmdWithAString         *shapeCmd;
  G4UIcmdWith3VectorAndUnit  *centreCmd;
  G4UIcmdWithADoubleAndUnit  *centre_rCmd;
  G4UIcmdWithADoubleAndUnit  *centre_phiCmd;
  G4UIcmdWithADoubleAndUnit  *centre_zCmd;
  G4UIcmdWithADoubleAndUnit  *halfzCmd;
  G4UIcmdWithADoubleAndUnit  *radiusCmd;
  G4UIcmdWithAString         *confineCmd;         
  G4UIcmdWithAString         *angtypeCmd;
  G4UIcmdWithAString         *energytypeCmd;
  G4UIcmdWithAnInteger       *verbosityCmd;
  G4UIcmdWithAnInteger       *CEflagCmd;
  G4UIcommand                *ionCmd;
  G4UIcmdWithAString         *particleCmd;
  G4UIcmdWith3VectorAndUnit  *positionCmd;
  G4UIcmdWith3Vector         *directionCmd;
  G4UIcmdWithADoubleAndUnit  *energyCmd;
  G4UIcmdWithoutParameter    *listCmd;
  
  G4ParticleTable *particleTable;
  G4bool   fShootIon; 			// Particle parameters
  G4int    fAtomicNumber;
  G4int    fAtomicMass;
  G4int    fIonCharge;
  G4double fIonExciteEnergy;
};

#endif	/* RadiogenicSourcesMessenger_hh */
