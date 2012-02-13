////////////////////////////////////////////////////////////////////////
// $Id: Co60Lines.cc,v 1.1 2011/05/24 21:50:31 kevmc Exp $
//  File:        Co60Lines.hh                                         //
//  Description: Generate special Co-60 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        24 May 2011                                          //
//                                                                    //
//  20110524  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Co60Lines.hh"


// Constructor fills spectrum

Co60Lines::Co60Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Co60Lines", dir, verbose) {
  AddLine(1173.24*keV, 0.999736);    
  AddLine(1332.50*keV, 0.999856); 
  AddLine(346.93*keV,  0.000076); 
  AddLine(826.06*keV,  0.000076); 
  AddLine(2158.57*keV, 0.0000111); 
  AddLine(2505*keV,    0.00000002); 
  AddLine(7.461*keV,   0.0000343); 
  AddLine(7.478*keV,   0.000067); 
  AddLine(8.265*keV,   0.00000413); 
  AddLine(8.265*keV,   0.0000081); 
  NormalizeLines();
}

Co60Lines::~Co60Lines() {}
