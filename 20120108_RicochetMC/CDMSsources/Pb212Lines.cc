////////////////////////////////////////////////////////////////////////
// $Id: Pb212Lines.cc,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Pb212Lines.hh                                         //
//  Description: Generate special Pb212 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Pb212Lines.hh"


// Constructor fills spectrum

Pb212Lines::Pb212Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Pb212Lines", dir, verbose) {
    	AddLine(115.183*keV, 0.592);
	AddLine(176.68*keV, 0.052);
	AddLine(238.632*keV, 43.3);
	AddLine(300.087*keV, 3.28);
	AddLine(415.2*keV, 0.143);
  NormalizeLines();
}

Pb212Lines::~Pb212Lines() {}
