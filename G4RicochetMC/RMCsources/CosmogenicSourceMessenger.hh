#ifndef CosmogenicSourceMessenger_hh
#define CosmogenicSourceMessenger_hh
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSourceMessenger.hh                         //
//  Description: User interface for cosmic ray generators             //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Dennis Wright (SLAC)                   //
//  Date:        13 February 2012                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4UImessenger.hh"
#include "globals.hh"

class CosmogenicSource;
class G4ParticleTable;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;


class CosmogenicSourceMessenger : public G4UImessenger {
  public:
    CosmogenicSourceMessenger(CosmogenicSource* fPtclGun);
    ~CosmogenicSourceMessenger();
  
  //void SetNewValue(G4UIcommand *command, G4String newValues);
  
  private:
    CosmogenicSource* fParticleGun;
    G4UIdirectory* gunDirectory;
    G4UIcmdWithADoubleAndUnit* depthCmd;
    G4UIcmdWithAnInteger* verbosityCmd;
    G4UIcmdWithoutParameter* genTestCmd;
};

#endif
