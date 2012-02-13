#ifndef CDMSEventSelector_hh
#define CDMSEventSelector_hh 1
////////////////////////////////////////////////////////////////////////
// $Id: CDMSEventSelector.hh,v 1.1 2011/06/29 06:21:05 kelsey Exp $
//                                                                    //
//  File:        CDMSEventSelector.hh                                 //
//  Description: Accept/reject interface for use in CDMSRunAction     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        28 June 2011                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "globals.hh"

class G4Event;


class CDMSEventSelector {
public:
  CDMSEventSelector(const G4String& name, G4int verbose=0)
    : selectorName(name), verboseLevel(verbose) {}

  virtual ~CDMSEventSelector() {}

  // Subclasses must implement selection function
  virtual G4bool accept(const G4Event* evt) const = 0;
  virtual G4bool reject(const G4Event* evt) const { return !accept(evt); }

  virtual void SetVerboseLevel(G4int verbose) { verboseLevel = verbose; }

  virtual const G4String& GetName() const { return selectorName; }

protected:
  G4String selectorName;
  G4int verboseLevel;
};

#endif	/* CDMSEventSelector_hh */
