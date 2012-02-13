////////////////////////////////////////////////////////////////////////
// $Id: Tl208Lines.cc,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Tl208Lines.hh                                         //
//  Description: Generate special Tl208 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Tl208Lines.hh"


// Constructor fills spectrum

Tl208Lines::Tl208Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Tl208Lines", dir, verbose) {
    	AddLine(211.40*keV, 0.178);
	AddLine(233.36*keV, 0.307);
	AddLine(252.61*keV, 0.69);
	AddLine(277.351*keV, 6.31);
	AddLine(485.95*keV, 0.050);
	AddLine(510.77*keV, 22.6);
	AddLine(583.191*keV, 84.5);
	AddLine(587.7*keV, 0.040);
	AddLine(650.1*keV, 0.036);
	AddLine(705.2*keV, 0.022);
	AddLine(722.04*keV, 0.201);
	AddLine(748.7*keV, 0.043);
	AddLine(763.13*keV, 1.81);
	AddLine(821.2*keV, 0.040);
	AddLine(860.564*keV, 12.42);
	AddLine(883.3*keV, 0.031);
	AddLine(927.6*keV, 0.131);
	AddLine(982.7*keV, 0.203);
	AddLine(1093.9*keV, 0.40);
	AddLine(1125.7*keV, 0.0050);
	AddLine(1160.8*keV, 0.011);
	AddLine(1185.1*keV, 0.017);
	AddLine(1282.8*keV, 0.052);
	AddLine(1381.1*keV, 0.007);
	AddLine(1647.5*keV, 0.0020);
	AddLine(1744.0*keV, 0.0020);
	AddLine(2614.533*keV, 99);
  NormalizeLines();
}

Tl208Lines::~Tl208Lines() {}
