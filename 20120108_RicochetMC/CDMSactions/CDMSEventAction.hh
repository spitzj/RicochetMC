//
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSEventAction.hh                               //
//  Description: CDMSMini event action class                          //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        14 May 2010                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSEventAction_h
#define CDMSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class CDMSRunAction;


class CDMSEventAction : public G4UserEventAction
{
  public:
    CDMSEventAction(CDMSRunAction*);
    virtual ~CDMSEventAction();

    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
  private:

    CDMSRunAction* runAction;
   
  //   CDMSEventActionMessenger*  eventMessenger;
};

#endif

    
