////////////////////////////////////////////////////////////////////////
// $Id: Am241Lines.cc,v 1.4 2011/05/21 05:23:22 kelsey Exp $
//  File:        Am241Lines.hh                                        //
//  Description: Generate special Am-241 gammas for CDMS testing      //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        27 April 2011                                        //
//                                                                    //
//  20110429  M. Kelsey -- Add verbosity, bug fix direction generator //
//  20110520  M. Kelsey -- Move functional code to new base class     //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Am241Lines.hh"


// Constructor fills spectrum

Am241Lines::Am241Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Am241Lines", dir, verbose) {
  AddLine(26.3448*keV, 0.024);    
  AddLine(32.183*keV,  0.000174); 
  AddLine(33.196*keV,  0.00126);  
  AddLine(43.423*keV,  0.00073);  
  AddLine(55.56*keV,   0.000181); 
  AddLine(59.541*keV,  0.359);    
  AddLine(98.97*keV,   0.000203); 
  AddLine(102.98*keV,  0.000195); 
  AddLine(11.871*keV,  0.0066);   
  AddLine(13.761*keV,  0.0107);   
  AddLine(13.946*keV,  0.096);    
  AddLine(15.861*keV,  0.00153);  
  AddLine(16.109*keV,  0.00184);  
  AddLine(16.816*keV,  0.025);    
  AddLine(17.061*keV,  0.015);    
  AddLine(17.505*keV,  0.0065);   
  AddLine(17.751*keV,  0.057);    
  AddLine(17.992*keV,  0.0137);   
  AddLine(20.784*keV,  0.0139);   
  AddLine(21.099*keV,  0.0065);   
  AddLine(21.342*keV,  0.0059);   
  AddLine(21.491*keV,  0.0029);   
  AddLine(97.069*keV,  0.00008);  
  AddLine(101.059*keV, 0.00012);

  NormalizeLines();
}

Am241Lines::~Am241Lines() {}
