#ifndef CDMSVesselConstruction_hh
#define CDMSVesselConstruction_hh 1
// $Id: CDMSVesselConstruction.hh,v 1.20 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVesselConstruction.hh                            //
//                                                                    //
//  Description: Construction of vacuum cryostat with supports        //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        16 November 2010                                     //
//                                                                    //
//  20101130  M. Kelsey -- Add reporting.                             //
//  20101208  M. Kelsey -- Move tower support structure to new class  //
//  20101209  M. Kelsey -- Pass verbosity to contained objects, add   //
//		support functions to build mother volume              //
//  20101210  M. Kelsey -- Follow base class change to dimensions     //
//  20101211  M. Kelsey -- Drop unused parameters                     //
//  20101215  M. Kelsey -- Add number of vacuum vessels, multi pipes  //
//  20101218  M. Kelsey -- Convert fixed length arrays to vectors.    //
//  20101223  M. Kelsey -- Add functions to return complete vectors.  //
//  20110103  M. Kelsey -- Add access to pipe parameter vectors.      //
//  20110105  M. Kelsey -- Rename pipe variables, add pipe len list;  //
//		add offsets from ends for cryo and vacuum pipes.      //
//  20110325  M. Kelsey -- Drop UseTower() functions; Manager instead //
//  20110722  M. Kelsey -- Remove materials table; use singleon in .cc//
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMSVDetectorGeometry.hh"
#include "globals.hh"
#include "G4Transform3D.hh"
#include <vector>

class G4LogicalVolume;
class G4VSolid;
class G4VisAttributes;
class CDMSGeometryManager;
class CDMSTowerConstruction;
class CDMSTowerSupportConstruction;
class CDMSTowerPattern;
class CDMSVesselMessenger;


class CDMSVesselConstruction : public CDMSVDetectorGeometry {
public:
  CDMSVesselConstruction();
  virtual ~CDMSVesselConstruction();

  virtual G4double GetRadius() const;	  // Radius of maximum extent
  virtual G4double GetLength() const;	  // Z-length of maximum extent

  virtual G4LogicalVolume* BuildGeometry();	// Construct for positioning
  virtual void FillExtraParameters();
  virtual void PrintParameters(std::ostream& os) const;

  void SetVerboseLevel(G4int verbose=0);

  void SetNStages(G4int value)                  { NStages = value; }
  void SetNVacuumStages(G4int value)            { NVacuum = value; }
  void SetVesselThick(G4int i, G4double value)  { VesselThick[i] = value; }
  void SetVessel0Rad(G4double value)            { Vessel0Rad = value; }
  void SetVesselDeltaRad(G4double value)        { VesselDeltaRad = value; }
  void SetVesselDeltaHeight(G4double value)     { VesselDeltaHeight = value; }
  void SetVesselExtraHeight(G4double value)     { VesselExtraHeight = value; }
  void SetVesselGap(G4double value)             { VesselGap = value; }
  void SetLidThickness(G4int i, G4double value) { LidThickness[i] = value; }
  void SetStemInnerRad(G4double value)          { StemInnerRad = value; }
  void SetVacInnerRad(G4double value)           { VacInnerRad = value; }
  void SetPipeThick(G4double value)             { PipeThick = value; }
  void SetPipeOutsideLen(G4double value)        { PipeOutsideLen = value; }

  void SetVesselThick(const std::vector<G4double>& value);
  void SetLidThickness(const std::vector<G4double>& value);

  G4int    GetNStages() const             { return NStages; }
  G4int    GetNVacuumStages() const       { return NVacuum; }
  G4double GetVesselThick(G4int i) const  { return VesselThick[i]; }
  G4double GetVessel0Rad() const          { return Vessel0Rad; }
  G4double GetVesselDeltaRad() const      { return VesselDeltaRad; }
  G4double GetVesselDeltaHeight() const   { return VesselDeltaHeight; }
  G4double GetVesselExtraHeight() const   { return VesselExtraHeight; }
  G4double GetVesselGap() const           { return VesselGap; }
  G4double GetLidThickness(G4int i) const { return LidThickness[i]; }
  G4double GetStemInnerRad() const        { return StemInnerRad; }
  G4double GetVacInnerRad() const         { return VacInnerRad; }
  G4double GetPipeThick() const           { return PipeThick; }
  G4double GetPipeOutsideLen() const      { return PipeOutsideLen; }

  G4double GetPipeLen(G4int i) const      { return pipeLen[i]; }
  G4double GetStemRad(G4int i) const      { return stemRad[i]; }
  G4double GetZStemHole(G4int i) const    { return zStemHole[i]; }
  G4double GetVacRad(G4int i) const       { return vacRad[i]; }
  G4double GetZVacHole(G4int i) const     { return zVacHole[i]; }
  
  const std::vector<G4double>& GetVesselHeight() const { return vesselHeight; }
  const std::vector<G4double>& GetVesselZPos() const   { return zVessel; }
  const std::vector<G4double>& GetVesselThick() const  { return VesselThick; }
  const std::vector<G4double>& GetLidThickness() const { return LidThickness; }

  const std::vector<G4double>& GetPipeLenList() const   { return pipeLen; }
  const std::vector<G4double>& GetStemRadList() const   { return stemRad; }
  const std::vector<G4double>& GetZStemHoleList() const { return zStemHole; }
  const std::vector<G4double>& GetVacRadList() const    { return vacRad; }
  const std::vector<G4double>& GetZVacHoleList() const  { return zVacHole; }

private:
  G4LogicalVolume* BuildEnvelope();	// Creates outermost "mother volume"

  // These functions place major physical volumes in the top-most mother
  void PlaceTowerSupports(G4LogicalVolume* world);
  void PlaceVesselStage(G4int stage, G4LogicalVolume* world);

  // These functions create nested volumes for placement
  G4LogicalVolume* BuildVesselTop(G4int stage);
  G4LogicalVolume* BuildVesselSide(G4int stage);
  G4LogicalVolume* BuildVesselBottom(G4int stage);
  G4LogicalVolume* BuildCryoPipe(G4int stage);
  G4LogicalVolume* BuildVacuumPipe(G4int stage);

  // These functions return the orientation matrix for the pipe stems
  G4Transform3D GetCryoPipePosition(G4int stage);
  G4Transform3D GetVacuumPipePosition(G4int stage);

  // Remove cylindrical spaces for each tower from mother's physical volume
  G4VSolid* SubtractTowerPositions(G4VSolid* icebox);

  void UpdateTowerPattern();

  CDMSGeometryManager* theManager;
  CDMSTowerPattern* towerPattern;
  CDMSTowerConstruction* towerBuilder;
  CDMSTowerSupportConstruction* supportBuilder;

  G4int    NStages;		// Total number of nested cryostats
  G4int    NVacuum;		// Number of vacuum-pumped (outer) cryostats
  G4int    NTowers;		// Number of towers (copied from towerPattern)

  std::vector<G4double> VesselThick;
  G4double Vessel0Rad;
  G4double VesselDeltaRad;
  G4double VesselDeltaHeight;
  G4double VesselExtraHeight;
  G4double VesselGap;		// Vertical distance between vessels 4 and 5

  std::vector<G4double> LidThickness;
  G4double StemInnerRad;
  G4double VacInnerRad;
  G4double StemZInset;
  G4double VacZInset;
  G4double PipeThick;
  G4double PipeOutsideLen;

  void CopyTowerParameters();	// Copy values from TowerConstruction to here
  G4int    NTowerSides;
  G4double ExtraSpace;
  G4double housingRadius;
  G4double housingThickness;
  G4double towerHeight;
  G4double towerRadius;

  void CopySupportParameters();
  G4int    supportNStages;
  G4int    supportMountStage;
  G4int    supportSquidStage;
  G4int    supportFETStage;
  G4double supportRadius;
  G4double supportBottom;
  G4double spoolLength;
  G4double stageHeight;
  std::vector<G4double> stageGaps;

  void FillVesselParameters();
  void ComputeVesselTopCoordinates();
  void ComputeVesselDimensions();
  void ComputePipeCoordinates();
  std::vector<G4double> vesselRadius;
  std::vector<G4double> vesselHeight;
  std::vector<G4double> zVesselTop;
  std::vector<G4double> zVessel;
  std::vector<G4double> zVesselBottom;
  std::vector<G4double> pipeLen;
  std::vector<G4double> stemRad;
  std::vector<G4double> zStemHole;
  std::vector<G4double> vacRad;
  std::vector<G4double> zVacHole;
  std::vector<G4VisAttributes*> vesselVisAtt;

  CDMSVesselMessenger* messenger;
};

#endif	/* CDMSVesselConstruction_hh */
