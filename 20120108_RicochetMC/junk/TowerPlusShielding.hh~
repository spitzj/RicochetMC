#ifndef CDMSShieldConstruction_hh
#define CDMSShieldConstruction_hh 1
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

class CDMSVesselConstruction;
class CDMSShieldMessenger;
class G4LogicalVolume;
class G4VisAttributes;
class G4VSolid;
class G4VSensitiveDetector;

class CDMSShieldConstruction : public CDMSVDetectorGeometry {
public:
  // Data structure to carry around parameters for each sheilding layer
  struct LayerData {
    G4double Gap, SideThick, TopThick, BottomThick;	// Preset
    G4double Radius, Length, Mid, Z;			// Computed
    G4String Name;					// For volumes
    G4String matName;					// For material
    G4VisAttributes* visAtt;				// Owned

    LayerData(const char* name="LayerData", G4double gap=0., G4double side=0.,
	      G4double top=0., G4double bottom=0., const char* mat="");
    ~LayerData();
    void FillOuterShield(G4double& radius, G4double& length, G4double& z);
    void FillInnerShield(G4double& radius, G4double& length, G4double& z);
    void Set(G4double g, G4double s, G4double t, G4double b) {
      Gap = g; SideThick = s; TopThick = t; BottomThick = b;
    }
    void Copy(const LayerData& layer) {
      Set(layer.Gap, layer.SideThick, layer.TopThick, layer.BottomThick);
    }
  };

public:
  CDMSShieldConstruction();
  virtual ~CDMSShieldConstruction();
  virtual G4double GetRadius() const { return mother.Radius; }
  virtual G4double GetLength() const { return mother.Length; }
  virtual G4LogicalVolume* BuildGeometry();
  virtual void FillExtraParameters();
  virtual void PrintParameters(std::ostream& os) const;

  virtual void SetVerboseLevel(G4int verbose=0);

  void UseIcebox();

  void SetMuMetalParams(const LayerData& buf)   {muMetal.Copy(buf);}
  void SetInnerPolyParams(const LayerData& buf) {innerPoly.Copy(buf);}
  void SetInnerLeadParams(const LayerData& buf) {innerLead.Copy(buf);}
  void SetOuterLeadParams(const LayerData& buf) {outerLead.Copy(buf);}
  void SetOuterPolyParams(const LayerData& buf) {outerPoly.Copy(buf);}

  void SetScintThick(G4double value) {scintThick = value;}

  void SetNTopPanelsX(G4int value) { NTopPanelsX = value; }
  void SetNTopPanelsY(G4int value) { NTopPanelsY = value; }
  void SetOverhang(G4double value) { overhang = value; }
  void SetOverlapX(G4double value) { overlapX = value; }
  void SetOverlapY(G4double value) { overlapY = value; }
  void SetTopVetoSupportHeight(G4double value) { topVetoSupportHeight = value; }
  void SetBottomVetoSupportHeight(G4double value) { bottomVetoSupportHeight = value; }

  void SetNSidePanels(G4int value);	// Will check for even-sided polygon
  void SetSidePanelClearance(G4double value) { sidePanelClearance = value; }
  void SetSidePanelOverlap(G4double value)   { sidePanelOverlap = value; }
  void SetSideCornerOverlap(G4double value)  { sideCornerOverlap = value; }

  const LayerData& GetMuMetalParams() const   {return muMetal;}
  const LayerData& GetInnerPolyParams() const {return innerPoly;}
  const LayerData& GetInnerLeadParams() const {return innerLead;}
  const LayerData& GetOuterLeadParams() const {return outerLead;}
  const LayerData& GetOuterPolyParams() const {return outerPoly;}

  G4double GetScintThick() const {return scintThick;}

  G4int    GetNTopPanelsX() const { return NTopPanelsX; }
  G4int    GetNTopPanelsY() const { return NTopPanelsY; }
  G4double GetOverhang() const { return overhang; }
  G4double GetOverlapX() const { return overlapX; }
  G4double GetOverlapY() const { return overlapY; }
  G4double GetTopVetoSupportHeight() const { return topVetoSupportHeight; }
  G4double GetBottomVetoSupportHeight() const { return bottomVetoSupportHeight; }

  G4int    GetNSidePanels() const { return NSidePanels; }
  G4double GetSidePanelClearance() const { return sidePanelClearance; }
  G4double GetSidePanelOverlap() const   { return sidePanelOverlap; }
  G4double GetSideCornerOverlap() const  { return sideCornerOverlap; }

private:
  G4LogicalVolume* CreateMotherVolume();

  void PlaceShielding(G4LogicalVolume* mother);
  void PlaceShieldLayer(const LayerData& layer, G4LogicalVolume* mother, G4bool innerType);
  G4LogicalVolume* BuildShieldLayer(const LayerData& layer, G4bool innerType);

  // Add pipe penetration (along Y axis) through solid at specified location
  G4VSolid* 
  MakePipeHoleInSolid(G4VSolid* solid, G4double radius, G4double thick,
		      G4double holeX, G4double holeY, G4double holeZ) {
    return MakePipeHoleInSolid(solid, radius, thick,
			       G4ThreeVector(holeX,holeY,holeZ));
  }

  G4VSolid* MakePipeHoleInSolid(G4VSolid* solid,
				G4double radius, G4double thick,
				const G4ThreeVector& hole);

  void PlaceVeto(G4LogicalVolume* mother);
  void PlaceEndcapVeto(G4LogicalVolume* mother);
  void PlaceBarrelVeto(G4LogicalVolume* mother);
  void PlaceTopSideVeto(G4LogicalVolume* mother);
  void PlaceMiddleSideVeto(G4LogicalVolume* mother);
  void PlaceBottomSideVeto(G4LogicalVolume* mother);

  // Identify and install pipe penetrations in side panels
  G4bool PanelIntersectsVacStem(G4int index, G4double& panelY);
  G4bool PanelIntersectsFridgeStem(G4int index, G4double& panelY);
  G4LogicalVolume* MakeVacHoleInPanel(G4VSolid* panel, G4int index,
				      G4double holeY);
  G4LogicalVolume* MakeFridgeHoleInPanel(G4VSolid* panel, G4int index,
				      G4double holeY);
  G4LogicalVolume* MakePipeHoleInPanel(G4VSolid* panel, G4int index,
				       G4double radius, G4double holeY,
				       G4double holeZ);

  void GetCornersOfPanel(G4int index, G4Point3D& right, G4Point3D& left);
  G4double DistanceToPipe(const G4Point3D& right, const G4Point3D& left);

  G4VSensitiveDetector* BuildSensitiveDetector();

  CDMSVesselConstruction* icebox;

  G4String DetectorName;                // Sensitive detector name

  void FillIceboxParameters();		// Copy parameters from cryostat
  G4double iceboxRad;
  G4double iceboxLength;
  G4double iceboxZ;
  G4double fridgeStemRad;
  G4double vacStemRad;
  G4double zFridgeStem;
  G4double zVacStem;

  void FillMotherParameters();		// Compute from outer veto positions
  void FillShieldParameters();
  LayerData mother;
  LayerData muMetal;
  LayerData innerPoly;
  LayerData innerLead;
  LayerData outerLead;
  LayerData outerPoly;

  G4double outerShieldRad;
  G4double outerShieldLength;
  G4double outerShieldZ;

  G4double scintThick;

  void FillVetoParameters();
  void FillEndcapVetoParameters();
  G4int    NTopPanelsX;
  G4int    NTopPanelsY;
  G4int    NTopPanels;	// == NTopPanelsX + NTopPanelsY
  G4double overhang;
  G4double overlapX;
  G4double overlapY;
  G4double topVetoSupportHeight;
  G4double bottomVetoSupportHeight;
  G4double topVetoLength;
  G4double topPanelLength;
  G4double topPanelWidth;
  std::vector<G4double> panelX;
  std::vector<G4double> panelY;
  std::vector<G4double> panelZtop;
  std::vector<G4double> panelZbot;

  void FillBarrelVetoParameters();
  G4int    NSidePanels;
  G4double sidePanelClearance;
  G4double sidePanelOverlap;
  G4double sideCornerOverlap;
  G4double sideRingRadius1;
  G4double sideRingRadius2;
  G4double sidePanelWidth;
  G4double sidePanelShift;
  G4double midPanelShift;
  G4double loPanelBottom;
  G4double loPanelTop;
  G4double loPanelHeight;
  G4double loPanelZ;
  G4double midPanelHeight;
  G4double midPanelZ;
  G4double hiPanelBottom;
  G4double hiPanelTop;
  G4double hiPanelHeight;
  G4double hiPanelZ;
  G4int    iVacPanel;			// Indices of panels with penetrations
  G4int    iFridgePanel;
  std::vector<G4Transform3D> topTransform;	// Placement for side panels
  std::vector<G4Transform3D> middleTransform;
  std::vector<G4Transform3D> bottomTransform;

  CDMSShieldMessenger* messenger;
};


// Write out shielding layer parameters

std::ostream& operator<<(std::ostream& os,
			 const CDMSShieldConstruction::LayerData& layer);

#endif
