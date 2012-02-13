////////////////////////////////////////////////////////////////////////
// $Id: Pb214Lines.cc,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Pb214Lines.hh                                         //
//  Description: Generate special Pb214 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Pb214Lines.hh"


// Constructor fills spectrum

Pb214Lines::Pb214Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Pb214Lines", dir, verbose) {
    	AddLine(53.2275*keV, 1.2);
	AddLine(107.22*keV, 0.015);
	AddLine(137.45*keV, 0.006);
	AddLine(141.3*keV, 0.004);
	AddLine(170.07*keV, 0.032);
	AddLine(196.20*keV, 0.069);
	AddLine(205.68*keV, 0.0115);
	AddLine(216.47*keV, 0.022);
    AddLine(238.4*keV, 0.015);	  
	AddLine(241.997*keV, 7.43);
	AddLine(258.87*keV, 0.524);
	AddLine(274.80*keV, 0.474);
	AddLine(295.224*keV, 19.3);
	AddLine(298.76*keV, 0.02);	  
	AddLine(305.26*keV, 0.031);
	AddLine(314.32*keV, 0.078);
	AddLine(323.83*keV, 0.028);
	AddLine(351.932*keV, 37.6);
	AddLine(462.00*keV, 0.221);
	AddLine(480.43*keV, 0.320);
	AddLine(487.09*keV, 0.422);
	AddLine(511.0*keV, 0.032);
	AddLine(533.66*keV, 0.186);
	AddLine(538.41*keV, 0.020);
	AddLine(543.81*keV, 0.069);
	AddLine(580.13*keV, 0.352);
	AddLine(765.96*keV, 0.078);
	AddLine(785.96*keV, 1.07);
	AddLine(839.04*keV, 0.587);

  NormalizeLines();
}

Pb214Lines::~Pb214Lines() {}
