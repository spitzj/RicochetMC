#ifndef TowerPlusShielding_hh
#define TowerPlusShielding_hh 1
// $Id: CDMSShieldConstruction.hh,v 1.19 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSShieldConstruction.hh                            //
//                                                                    //
//  Description: CDMS shielding and veto detector surrounding the     //
//               icebox                                               //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        3 December 2010                                      //
//                                                                    //
//  20110105  M. Kelsey -- Add Messenger.                             //
//  20110106  M. Kelsey -- Define struct for layer parameters, use in //
//		BuildShieldLayer.                                     //
//  20110107  M. Kelsey -- Move veto parameters to data members, set  //
//		in FillExtraParameters().  New function to make pipe  //
//		penetrations.                                         //
//  20110110  M. Kelsey -- Pipe hole function needs full three-vector //
// 		Add parameters for number of veto panels              //
//  20110111  M. Kelsey -- Set shield-layer parameters en masse.      //
//  20110112  M. Kelsey -- Allow for "default" LayerData constructor, //
//		add LayerData::Copy() function for input numerics.    //
//  20110113  M. Kelsey -- Check for even-sided polygon veto.         //
//  20110114  M. Kelsey -- Handle pipe holes in arbitrary veto panels //
//  20110328  M. Kelsey -- Access cryostat via Manager                //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "globals.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>

class CDMSShieldConstruction;
class CDMSTowerConstruction;
//class G4LogicalVolume;
//class G4VisAttributes;
//class G4VSolid;
//class G4VSensitiveDetector;

class TowerPlusShielding : public CDMSVDetectorGeometry {
public:
  TowerPlusShielding();
  virtual ~TowerPlusShielding();
  virtual G4double GetRadius() const { return theShield->GetRadius(); }
  virtual G4double GetLength() const { return theShield->GetLength(); }
  virtual void SetVerboseLevel(G4int verbose);
  virtual void FillExtraParameters();
  virtual G4LogicalVolume* BuildGeometry();
  virtual void PrintParameters(std::ostream& os);

private:
  CDMSShieldConstruction* theShield;
  CDMSTowerConstruction* theTower;
  

#endif
