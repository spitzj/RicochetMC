//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSMiniGeomMessenger.hh                             //
//                                                                    //
//  Description: Messenger class to allow setting CDMS mini           //
//               geometry parameters                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        23 June 2010                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSMiniGeomMessenger_h
#define CDMSMiniGeomMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class CDMSMiniDetConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;


class CDMSMiniGeomMessenger: public G4UImessenger
{
  public:
    CDMSMiniGeomMessenger(CDMSMiniDetConstruction*);
   ~CDMSMiniGeomMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    CDMSMiniDetConstruction* miniSetup;
    
    G4UIdirectory* miniGeomDir;
    G4UIcmdWithADoubleAndUnit* ZipRadCmd;
    G4UIcmdWithADoubleAndUnit* ZipLenCmd;
    G4UIcmdWithADoubleAndUnit* DetBoxShimCmd;
};

#endif

