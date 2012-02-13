////////////////////////////////////////////////////////////////////////
// $Id: Bi214Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Bi214Lines.hh                                         //
//  Description: Generate special Bi214 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Bi214Lines_hh
#define Bi214Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Bi214Lines : public CDMSGammaLines {
public:
  Bi214Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Bi214Lines();
};

#endif	/*  Bi214Lines_hh */
