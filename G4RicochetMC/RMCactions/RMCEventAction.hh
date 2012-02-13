////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCEventAction.cc                                    //
//  Description: Event action class                                   //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               (adapted from Dennis Wright (SLAC))                  //
//  Date:        8 January 2012                                       //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef RMCEventAction_h
#define RMCEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class RMCRunAction;


class RMCEventAction : public G4UserEventAction
{
  public:
    RMCEventAction(RMCRunAction*);
    virtual ~RMCEventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
  private:
    RMCRunAction* runAction;
};

#endif

    
