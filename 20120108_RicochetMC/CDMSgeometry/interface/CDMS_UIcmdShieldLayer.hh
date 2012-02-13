#ifndef CDMS_UIcmdShieldLayer_hh
#define CDMS_UIcmdShieldLayer_hh 1
// $Id: CDMS_UIcmdShieldLayer.hh,v 1.1 2011/01/12 19:50:42 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMS_UIcmdShieldLayer.hh                             //     
//  Description: Messenger interface to fill ShieldConstr::LayerData  //
//                                                                    //
//  This command sets the parameter names automatically, and provides //
//  built-in guidance.  Client code does not need to SetGuidance().   //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        12 January 2011                                      //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "G4UIcommand.hh"
#include "globals.hh"
#include "CDMSgeometry/detectors/CDMSShieldConstruction.hh"


class CDMS_UIcmdShieldLayer : public G4UIcommand {
public:
  CDMS_UIcmdShieldLayer(const char* theCommandPath,
			G4UImessenger* theMessenger);
  virtual ~CDMS_UIcmdShieldLayer() {}

  // Convert parameter list into form with default units
  virtual G4int DoIt(G4String parameterList);

  // Convert input string into shield-layer parameters using input units
  const CDMSShieldConstruction::LayerData&
  GetNewLayerData(const char* paramString);

  // Convert input string into shield-layer parameters in user's buffer
  void GetNewLayerData(const char* paramString,
		       CDMSShieldConstruction::LayerData& lyr);

  // Extract just the units portion of the input string (fifth argument)
  G4double GetNewUnitValue(const char* paramString);

  // Construct input data to string (e.g., for printing)
  G4String
  ConvertToString(const CDMSShieldConstruction::LayerData& lyr,
		  const char* unit);

  G4String
  ConvertToStringWithBestUnit(const CDMSShieldConstruction::LayerData& lyr);

  G4String 
  ConvertToStringWithDefaultUnit(const CDMSShieldConstruction::LayerData& lyr);

  // Assign name to data buffer, not to individual parameters
  void SetParameterName(const char* theName,
			G4bool omittable, G4bool currentAsDefault=false);

  void SetDefaultValue(const CDMSShieldConstruction::LayerData& defVal);

  // Set units configuration -- DefaultUnit is exclusive w.r.t. other two
  void SetUnitCategory(const char* unitCategory);
  void SetUnitCandidates(const char* candidateList);
  void SetDefaultUnit(const char* defUnit);

private:
  size_t Tokenize(const G4String& parameterList);	// Copy to buffer
  std::vector<G4String> tokens;

  CDMSShieldConstruction::LayerData layerBuffer;	// Parameter values

private:
  G4UIparameter* gapParam;	// Parameters for each input value
  G4UIparameter* sideParam;	// NOTE:  These are owned by base class!
  G4UIparameter* topParam;
  G4UIparameter* bottomParam;
  G4UIparameter* unitParam;
};

#endif	/* CDMS_UIcmdShieldLayer_hh */
