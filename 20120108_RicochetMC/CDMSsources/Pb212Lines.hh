////////////////////////////////////////////////////////////////////////
// $Id: Pb212Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Pb212Lines.hh                                         //
//  Description: Generate special Pb212 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Pb212Lines_hh
#define Pb212Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Pb212Lines : public CDMSGammaLines {
public:
  Pb212Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Pb212Lines();
};

#endif	/*  Pb212Lines_hh */
