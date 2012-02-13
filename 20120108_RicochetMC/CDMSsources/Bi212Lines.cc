////////////////////////////////////////////////////////////////////////
// $Id: Bi212Lines.cc,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Bi212Lines.hh                                         //
//  Description: Generate special Bi212 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Bi212Lines.hh"


// Constructor fills spectrum

Bi212Lines::Bi212Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Bi212Lines", dir, verbose) {
    	AddLine(39.858*keV, 1.091);
	AddLine(144.0*keV, 0.010);
	AddLine(164.0*keV, 0.005);
	AddLine(180.2*keV, 0.0032);
	AddLine(288.07*keV, 0.31);
	AddLine(295.1*keV, 0.024);
	AddLine(327.96*keV, 0.139);
	AddLine(433.6*keV, 0.012);
	AddLine(452.83*keV, 0.31);
	AddLine(473.6*keV, 0.046);
	AddLine(492.7*keV, 0.006);
	AddLine(576.0*keV, 0.0008);
	AddLine(620.4*keV, 0.0036);
	AddLine(727.330*keV, 6.58);
	AddLine(785.37*keV, 1.102);
	AddLine(893.408*keV, 0.378);
	AddLine(952.120*keV, 0.17);
	AddLine(1073.6*keV, 0.0160);
	AddLine(1078.62*keV, 0.564);
	AddLine(1512.7*keV, 0.29);
	AddLine(1620.50*keV, 1.49);
	AddLine(1679.7*keV, 0.058);
	AddLine(1806.0*keV, 0.090);

  NormalizeLines();
}

Bi212Lines::~Bi212Lines() {}
