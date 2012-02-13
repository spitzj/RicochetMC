////////////////////////////////////////////////////////////////////////
// $Id: Bi212Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Bi212Lines.hh                                         //
//  Description: Generate special Bi212 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Bi212Lines_hh
#define Bi212Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Bi212Lines : public CDMSGammaLines {
public:
  Bi212Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Bi212Lines();
};

#endif	/*  Bi212Lines_hh */
