//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        NoLabMessenger.hh                                    //
//                                                                    //
//  Description: Messenger class to allow setting NoLab geometry      //
//               parameters                                           //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        23 June 2010                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef NoLabMessenger_h
#define NoLabMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NoLab;
class G4UIdirectory;


class NoLabMessenger: public G4UImessenger
{
  public:
    NoLabMessenger(NoLab*);
   ~NoLabMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    NoLab* labSetup;
    G4UIdirectory* labDir;
};

#endif

