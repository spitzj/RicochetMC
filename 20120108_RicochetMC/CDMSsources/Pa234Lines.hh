////////////////////////////////////////////////////////////////////////
// $Id: Pa234Lines.hh,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Pa234Lines.hh                                         //
//  Description: Generate special Pa234 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Pa234Lines_hh
#define Pa234Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Pa234Lines : public CDMSGammaLines {
public:
  Pa234Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Pa234Lines();
};

#endif	/*  Pa234Lines_hh */
