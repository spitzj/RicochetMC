#ifndef CDMS_SelectMISS_hh
#define CDMS_SelectMISS_hh 1
////////////////////////////////////////////////////////////////////////
// $Id: CDMS_SelectMISS.hh,v 1.1 2011/07/08 21:54:11 kelsey Exp $
//                                                                    //
//  File:        CDMS_SelectMISS.hh                                   //
//  Description: Find events containing multiple scatters per crystal //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        8 July 2011                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSactions/CDMSEventSelector.hh"
#include "CDMSgeometry/detectors/CDMSZipHit.hh"	/* CDMSZipHitCollection */
#include "globals.hh"

class G4Event;


class CDMS_SelectMISS : public CDMSEventSelector {
public:
  CDMS_SelectMISS(G4int verbose=0) : CDMSEventSelector("SelectMISS") {}
  virtual ~CDMS_SelectMISS() {}

  virtual G4bool accept(const G4Event* evt) const;

protected:
  CDMSZipHitsCollection* getHits(const G4Event* evt) const;
};

#endif	/* CDMS_SelectMISS_hh */
