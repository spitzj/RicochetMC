////////////////////////////////////////////////////////////////////////
// $Id: K40Lines.cc,v 1.1 2011/05/24 21:50:31 kevmc Exp $
//  File:        K40Lines.hh                                          //
//  Description: Generate special K-40 gammas for CDMS testing        //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        24 May 2011                                          //
//                                                                    //
//  20110524  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/K40Lines.hh"


// Constructor fills spectrum

K40Lines::K40Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("K40Lines", dir, verbose) {
  AddLine(1460.83*keV, 0.11);    
  AddLine(2.955*keV,   0.0031); 
  AddLine(2.957*keV,   0.0062); 
  AddLine(3.190*keV,   0.00025); 
  AddLine(3.190*keV,   0.00049); 
  NormalizeLines();
}

K40Lines::~K40Lines() {}
