// $Id: CKGDetConstruction.hh,v 1.7 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CKGDetConstruction.hh                                //     
//  Description: Full multiple tower detector with icebox             //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        27 October 2010                                      //
//                                                                    //
//  20101028  M. Kelsey -- Follow CDMSVDetectorGeometry changes       //
//  20101210  M. Kelsey -- Follow CDMSVDetectorGeometry changes       //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#ifndef CKGDetConstruction_hh
#define CKGDetConstruction_hh 1

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "G4RotationMatrix.hh"
#include <vector>

class G4Region;
class G4LogicalVolume;
class G4Material;


class CKGDetConstruction: public CDMSVDetectorGeometry
{
  public:
  
    CKGDetConstruction();
    virtual ~CKGDetConstruction();

    G4double GetRadius() const { return 2*m; }	// Copy from "DetBox" in .cc
    G4double GetLength() const { return 4*m; }
    G4LogicalVolume* BuildGeometry();

    void SetZipRad(G4double);
    void SetZipLen(G4double);
    void SetDetBoxShimThick(G4double);

    G4String DetCollName;

  private:
    G4LogicalVolume* BuildZip(G4Material*);
    G4LogicalVolume* BuildZipHousing(G4double length);
    G4LogicalVolume* BuildSupportTower(G4double radius, G4double length);

    G4LogicalVolume* BuildTowerLid(G4bool invert);
    G4LogicalVolume* BuildStage();
    G4LogicalVolume* BuildSideStrip(G4double stripLength);
 
    // NOTE:  Vector passed by value because function changes contents
    G4LogicalVolume* BuildDiskWithHoles(G4double radius, G4double thick,
                                        std::vector<G4ThreeVector>,
                                        G4RotationMatrix*);

    G4LogicalVolume* BuildVessel(G4double tubeRadius, G4double tubeHeight,
                                 G4double tubeThick, G4double hole1Rad,
                                 G4double zHole1, G4double hole2Rad = 0.0,
                                 G4double zHole2 = 0.0);
 
    G4int NTowers;
    G4int NTowerSides;
    G4int NStages;
    
    G4double ZipRad;
    G4double ZipThick;
    G4double ZipGap;

    G4double ZipAxis1Len;
    G4double ZipAxis2Len;
    G4double Zip_z;
    G4double ZipClearance;
    G4double LidClearance;
    G4double HousingThickness;
    G4double ExtraSpace;

    G4double SpoolLen; 
    G4double SpoolIR;
    G4double SpoolThickness;
    G4double CTubeIR;
    G4double CTubeThick;

    G4double StageHeight;
    G4double StageWallThick;
    G4double StageGap;

    G4double SquidLen;
    G4double SquidThick;
    G4double SquidWidth;
    G4double FETLen;

    G4double VesselThick[6];
    G4double Vessel0Rad;
    G4double VesselDeltaRad;
    G4double VesselDeltaHeight;
    G4double VesselExtraHeight;
    G4double VesselGap;
    G4double LidThickness[3];

    G4double StemRad;
    G4double VacRad;

    G4double PipeThick;
    G4double PipeBaseLen;

    G4double StripThick;
    G4double StripWidth;

    G4double eps;

    G4Region* detectorRegion;
  //    CDMSMiniGeomMessenger* geomMessenger;
};

#endif
