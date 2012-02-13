#ifndef RMCTowerSupportConstruction_hh
#define RMCTowerSupportConstruction_hh 1
// $Id: RMCTowerSupportConstruction.hh,v 1.10 2011/07/22 21:09:17 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCTowerSupportConstruction.hh                       //
//                                                                    //
//  Description: Construction of tower support structures in cryostat //
//		 Internal coordinates put Z=0 at spool/CF joint               //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "RMCg4base/RMCVDetectorGeometry.hh"
#include "globals.hh"
#include <vector>

class G4LogicalVolume;
class G4VSolid;
class G4VisAttributes;
class RMCTowerConstruction;
class RMCTowerSupportMessenger;


class RMCTowerSupportConstruction : public RMCVDetectorGeometry {
public:
  RMCTowerSupportConstruction();
  virtual ~RMCTowerSupportConstruction() {}

  virtual G4double GetRadius() const { return supportRadius; }
  virtual G4double GetLength() const { return supportHeight; }

  virtual G4LogicalVolume* BuildGeometry();	// Construct for positioning
  virtual void FillExtraParameters();
  virtual void PrintParameters(std::ostream& os) const;

  void SetNStages(G4int value)           { NStages = value; }
  void SetStageHeight(G4double value)    { StageHeight = value; }
  void SetStageWallThick(G4double value) { StageWallThick = value; }
  void SetStageGap(G4double value)       { StageGap = value; }
  void SetVesselMountStage(G4int value)	 { VesselMountStage = value; }
  void SetSpoolLength(G4double value)    { SpoolLength = value; } 
  void SetSpoolIR(G4double value)        { SpoolIR = value; }
  void SetSpoolThickness(G4double value) { SpoolThickness = value; }
  void SetCTubeIR(G4double value)        { CTubeIR = value; }
  void SetCTubeThick(G4double value)     { CTubeThick = value; }
  void SetSquidLength(G4double value)    { SquidLength = value; }
  void SetSquidThick(G4double value)     { SquidThick = value; }
  void SetSquidWidth(G4double value)     { SquidWidth = value; }
  void SetSquidStage(G4int value)	 { SquidStage = value; }
  void SetFETLength(G4double value)      { FETLength = value; }
  void SetFETStage(G4int value)		 { FETStage = value; }

  G4int    GetNStages() const          { return NStages; }
  G4double GetStageHeight() const      { return StageHeight; }
  G4double GetStageWallThick() const   { return StageWallThick; }
  G4double GetStageGap() const         { return StageGap; }
  G4int    GetVesselMountStage() const { return VesselMountStage; }
  G4double GetSpoolLength() const      { return SpoolLength; } 
  G4double GetSpoolIR() const          { return SpoolIR; }
  G4double GetSpoolThickness() const   { return SpoolThickness; }
  G4double GetCTubeIR() const          { return CTubeIR; }
  G4double GetCTubeThick() const       { return CTubeThick; }
  G4double GetCTubeLength() const      { return CTubeLength; }
  G4double GetSquidLength() const      { return SquidLength; }
  G4double GetSquidThick() const       { return SquidThick; }
  G4double GetSquidWidth() const       { return SquidWidth; }
  G4int    GetSquidStage() const       { return SquidStageIndex; }
  G4double GetFETLength() const        { return FETLength; }
  G4int    GetFETStage() const         { return FETStageIndex; }

  G4double GetSupportHeight() const    { return supportHeight; }
  G4double GetSupportRadius() const    { return supportRadius; }
  G4double GetSupportTop() const       { return supportTop; }
  G4double GetSupportBottom() const    { return supportBottom; }
  G4double GetStageGapSum() const      { return StageGapSum; }
  const std::vector<G4double>& GetStageGaps() const { return StageGaps; }

private:
  // These functions position the components of the tower support
  void AttachCopperSpool(G4LogicalVolume* theSupport);
  void AttachSupportTube(G4LogicalVolume* theSupport);
  void AttachSupportStage(G4int stage, G4LogicalVolume* theSupport);
  void AttachReadout(G4int side, G4LogicalVolume* theSupport);

  // These functions create nested volumes for placement
  G4LogicalVolume* BuildSquidCard();
  G4LogicalVolume* BuildFETCard();
  G4LogicalVolume* BuildSideStrip(G4double length);

  RMCTowerConstruction* towerBuilder;

  G4int    NStages;		// Will override by RMCVesselConstruction
  G4double StageHeight;
  G4double StageWallThick;
  G4double StageGap;
  G4int    VesselMountStage;	// Index of stage mounted to inner cryostat

  G4double SpoolLength; 
  G4double SpoolIR;
  G4double SpoolThickness;
  G4double CTubeIR;
  G4double CTubeThick;

  G4double SquidLength;
  G4double SquidThick;		// Also used for FET thickness
  G4double SquidWidth;		// Also used for FET width
  G4int    SquidStage;		// Index of stage gap where SQUID is mounted
  G4int    SquidStageIndex;
  G4double FETLength;
  G4int    FETStage;		// Index of stage gap where FET card is mounted
  G4int    FETStageIndex;

  void CopyTowerParameters();	// Copy values from TowerConstruction to here
  G4int    NTowerSides;
  G4double housingRadius;
  G4double stripThick;
  G4double stripWidth;

  void FillSupportParameters();
  G4double supportHeight;
  G4double supportRadius;
  G4double CTubeLength;
  std::vector<G4double> StageGaps;
  G4double StageGapSum;
  G4double supportOffset;	// Reference Z coordinates for construction
  G4double supportBottom;
  G4double supportTop;

  void FillReadoutParameters();
  G4double cardRadius;
  G4double stripRadius;
  G4double squidHeight; 
  G4double fetHeight;
  G4double stripLoLen;
  G4double stripLoHeight;
  G4double stripMidLen;
  G4double stripMidHeight;
  G4double stripHiLen;
  G4double stripHiHeight;

  RMCTowerSupportMessenger* messenger;
};

#endif	/* RMCTowerSupportConstruction_hh */
