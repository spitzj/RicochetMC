#ifndef CosmogenicSourceMessenger_hh
#define CosmogenicSourceMessenger_hh
// $Id: CosmogenicSourceMessenger.hh,v 1.2 2011/06/30 21:35:16 dwright Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CosmogenicSourceMessenger.hh                         //
//  Description: User interface for cosmic ray generators             //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        28 May 2011                                          //
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
