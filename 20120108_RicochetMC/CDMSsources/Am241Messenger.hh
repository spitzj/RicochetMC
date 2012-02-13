//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        Am241Messenger.hh                                    //
//                                                                    //
//  Description: messenger class for Am241 source                     //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        2 September 2010                                     //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef Am241Messenger_h
#define Am241Messenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Am241Source;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;


class Am241Messenger: public G4UImessenger
{
  public:
    Am241Messenger(Am241Source*);
   ~Am241Messenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Am241Source* source;
    
    G4UIdirectory* Am241Dir;
    G4UIcmdWithADoubleAndUnit* SrcRadCmd;
    G4UIcmdWithADoubleAndUnit* SrcLenCmd;
    G4UIcmdWith3VectorAndUnit* SrcPosCmd;
    G4UIcmdWith3Vector* SrcDirCmd;
};

#endif

