////////////////////////////////////////////////////////////////////////
// $Id: Tl208Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Tl208Lines.hh                                         //
//  Description: Generate special Tl208 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Tl208Lines_hh
#define Tl208Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Tl208Lines : public CDMSGammaLines {
public:
  Tl208Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Tl208Lines();
};

#endif	/*  Tl208Lines_hh */
