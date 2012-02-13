////////////////////////////////////////////////////////////////////////
// $Id: Ac228Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Ac228Lines.hh                                         //
//  Description: Generate special Ac228 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Ac228Lines_hh
#define Ac228Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Ac228Lines : public CDMSGammaLines {
public:
  Ac228Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Ac228Lines();
};

#endif	/*  Ac228Lines_hh */
