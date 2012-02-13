////////////////////////////////////////////////////////////////////////
// $Id: Pb214Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Pb214Lines.hh                                         //
//  Description: Generate special Pb214 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Pb214Lines_hh
#define Pb214Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Pb214Lines : public CDMSGammaLines {
public:
  Pb214Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Pb214Lines();
};

#endif	/*  Pb214Lines_hh */
