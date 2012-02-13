// $Id: CDMSMiniDetConstruction.hh,v 1.6 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSMiniDetConstruction.hh                           //     
//  Description: single zip detector                                  //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        23 May 2010                                          //
//                                                                    //
//  20101026  M. Kelsey -- Add GetMaximumSize() function, follow name //
//		changes in base class                                 //
//  20101210  M. Kelsey -- Follow base class change to dimensions,    //
//		make SetX() functions inline, set location with Z     //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#ifndef CDMSMiniDetConstruction_h
#define CDMSMiniDetConstruction_h 1

#include "CDMSg4base/CDMSVDetectorGeometry.hh"

class CDMSMiniGeomMessenger;
class G4Region;
class G4LogicalVolume;


class CDMSMiniDetConstruction: public CDMSVDetectorGeometry {
public:
  CDMSMiniDetConstruction();
  virtual ~CDMSMiniDetConstruction();
  
  virtual G4double GetRadius() const { return Zip_Rout; }
  virtual G4double GetLength() const { return Zip_z; }
  virtual G4LogicalVolume* BuildGeometry();
  
  void SetZipRad(G4double val) { Zip_Rout = val; }
  void SetZipLen(G4double val) { Zip_z = val; position.setZ(val/2.); }
  void SetDetBoxShimThick(G4double val) { DetBoxShim = val; }
  
  G4String DetCollName;
  
private:
  G4double Zip_Rout;
  G4double Zip_z;
  G4double DetBoxShim;
  
  G4Region* detectorRegion;
  CDMSMiniGeomMessenger* geomMessenger;
};

#endif
