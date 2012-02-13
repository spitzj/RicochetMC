// $Id: Am241Source.hh,v 1.8 2011/07/22 21:08:36 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        Am241Source.hh                                       //
//  Description: specific source class for CDMS                       //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        2 September 2010                                     //
//                                                                    //
//  20101026  M. Kelsey -- Add GetMaximumSize()                       //
//  20101129  M. Kelsey -- Move particleGun to base class             //
//  20101209  M. Kelsey -- Reorder data members, inline SetXXX()      //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
//  20110427  M. Kelsey -- Move gamma lines to separate class.        //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#ifndef Am241Source_hh
#define Am241Source_hh 1

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class Am241Messenger;
class Am241Lines;
class G4LogicalVolume;


class Am241Source: public CDMSVSourceConstruction {
public:
  Am241Source();    
  virtual ~Am241Source();
  
  void GeneratePrimaries(G4Event*);

  virtual G4double GetRadius() const { return R_active; }
  virtual G4double GetLength() const { return L_active; }

  virtual G4LogicalVolume* BuildGeometry();
  virtual void PrintParameters(std::ostream& os) const;
  
  void SetRadius(G4double value) { R_active = value; }
  void SetLength(G4double value) { L_active = value; }
  
private:
  G4double R_active;       	// Radius of active region of source
  G4double L_active;       	// Length of active region of source

  Am241Lines* gammaLines;	// Gamma spectrum
  
  Am241Messenger* sourceMessenger;
};

#endif	/*  Am241Source_hh */
