////////////////////////////////////////////////////////////////////////
// $Id: Am241Lines.hh,v 1.4 2011/05/21 05:23:22 kelsey Exp $
//  File:        Am241Lines.hh                                        //
//  Description: Generate special Am-241 gammas for CDMS testing      //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        27 April 2011                                        //
//                                                                    //
//  20110429  M. Kelsey -- Add verbosity for diagnostics              //
//  20110520  M. Kelsey -- Move functional code to new base class     //
//////////////////////////////////////////////////////////////////////// 

#ifndef Am241Lines_hh
#define Am241Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Am241Lines : public CDMSGammaLines {
public:
  Am241Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Am241Lines();
};

#endif	/*  Am241Lines_hh */
