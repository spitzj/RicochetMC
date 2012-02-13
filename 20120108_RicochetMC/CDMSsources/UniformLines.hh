////////////////////////////////////////////////////////////////////////
// $Id: UniformLines.hh,v 1.1 2011/07/06 22:18:58 kevmc Exp $
//  File:        UniformLines.hh                                      //
//  Description: Generate uniform gammas 1-2500 keV for CDMS testing  //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        25 Jun 2011                                          //
//                                                                    //
//  20110625  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef UniformLines_hh
#define UniformLines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class UniformLines : public CDMSGammaLines {
public:
  UniformLines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~UniformLines();
};

#endif	/*  UniformLines_hh */
