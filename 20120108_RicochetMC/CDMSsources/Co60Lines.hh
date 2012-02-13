////////////////////////////////////////////////////////////////////////
// $Id: Co60Lines.hh,v 1.1 2011/05/24 21:50:31 kevmc Exp $
//  File:        Co60Lines.hh                                         //
//  Description: Generate special Co-60 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy                                       //
//  Date:        24 May 2011                                          //
//                                                                    //
//  20110524  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#ifndef Co60Lines_hh
#define Co60Lines_hh 1

#include "CDMSsources/CDMSGammaLines.hh"


class Co60Lines : public CDMSGammaLines {
public:
  Co60Lines(const G4ThreeVector& dir=origin, G4int verbose=0);    
  virtual ~Co60Lines();
};

#endif	/*  Co60Lines_hh */
