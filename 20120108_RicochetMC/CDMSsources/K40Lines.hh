////////////////////////////////////////////////////////////////////////
// $Id: K40Lines.hh,v 1.1 2011/05/24 21:50:31 kevmc Exp $
//  File:        K40Lines.hh                                          //
//  Description: Generate special K-40 gammas for CDMS testing        //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        24 May 2011                                          //
//                                                                    //
//  20110524  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef K40Lines_hh
#define K40Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class K40Lines : public CDMSGammaLines {
public:
  K40Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~K40Lines();
};

#endif	/*  K40Lines_hh */
