// $Id: CDMS_UIcmdShieldLayer.cc,v 1.3 2011/01/22 17:26:34 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMS_UIcmdShieldLayer.cc                             //     
//  Description: Messenger interface to fill ShieldConstr::LayerData  //
//                                                                    //
//  NOTE:  Some of the code is adapted from G4UIcmdWith3VectorAndUnit //
//  with simplifications.  Guidance is predefined in constructor.     //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        12 January 2011                                      //
//                                                                    //
//  20110122 M. Kelsey -- Bug fix; "i" was uninitialized at line 139  //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSgeometry/interface/CDMS_UIcmdShieldLayer.hh"
#include "G4UImessenger.hh"
#include "G4UIparameter.hh"
#include "G4Tokenizer.hh"


// Constructor and destructor

CDMS_UIcmdShieldLayer::CDMS_UIcmdShieldLayer(const char* theCommandPath,
					     G4UImessenger* theMessenger)
  : G4UIcommand(theCommandPath, theMessenger),
    gapParam(new G4UIparameter("Gap",'d',false)),
    sideParam(new G4UIparameter("Side",'d',false)),
    topParam(new G4UIparameter("Top",'d',false)),
    bottomParam(new G4UIparameter("Bottom",'d',false)),
    unitParam(new G4UIparameter("Unit",'s',false))
{
  SetParameter(gapParam);
  SetParameter(sideParam);
  SetParameter(topParam);
  SetParameter(bottomParam);
  SetParameter(unitParam);

  SetGuidance("Dimensions for CDMS shielding layer (Gap Side Top Bottom)");
  SetGuidance("Gap    : separation between interior object and this layer");
  SetGuidance("Side   : thickness of side wall for this layer");
  SetGuidance("Top    : thickness of upper wall for this layer");
  SetGuidance("Bottom : thickness of lower wall for this layer");
  SetGuidance("Units may be specified at end of list");
}


// Convert parameter list into form with default units

G4int CDMS_UIcmdShieldLayer::DoIt(G4String parameterList) {
  // convert a value from input to default units, overwriting input string
  if (Tokenize(parameterList)) {
    G4String default_unit = unitParam->GetDefaultValue();
    if (default_unit != "" && tokens.size() >= 5) {
      G4double scale = ValueOf(tokens[4]) / ValueOf(default_unit);
      layerBuffer.Gap         = ConvertToDouble(tokens[0]) * scale;
      layerBuffer.SideThick   = ConvertToDouble(tokens[1]) * scale;
      layerBuffer.TopThick    = ConvertToDouble(tokens[2]) * scale;
      layerBuffer.BottomThick = ConvertToDouble(tokens[3]) * scale;
      
      // reconstruct parameter list, including any excess tokens
      parameterList = ConvertToStringWithDefaultUnit(layerBuffer);
      for (size_t i=4; i<tokens.size(); i++) {
	parameterList += " " + tokens[i];
      }
    }
  }

  return G4UIcommand::DoIt(parameterList);	// Reprocess through base
}


// Convert input string into shield-layer parameters using input units

const CDMSShieldConstruction::LayerData&
CDMS_UIcmdShieldLayer::GetNewLayerData(const char* paramString) {
  layerBuffer.Set(0.,0.,0.,0.);

  if (Tokenize(paramString) >= 4) {		// Only complete list allowed
    G4double scale = 1.;
    if (tokens.size() >= 5) scale = ValueOf(tokens[4]);
    else {
      G4String default_unit = unitParam->GetDefaultValue();
      scale = ValueOf(default_unit);
    }

    layerBuffer.Gap         = ConvertToDouble(tokens[0]) * scale;
    layerBuffer.SideThick   = ConvertToDouble(tokens[1]) * scale;
    layerBuffer.TopThick    = ConvertToDouble(tokens[2]) * scale;
    layerBuffer.BottomThick = ConvertToDouble(tokens[3]) * scale;
  }    

  return layerBuffer;
}

void CDMS_UIcmdShieldLayer::
GetNewLayerData(const char* paramString,
		CDMSShieldConstruction::LayerData& lyr) {
  GetNewLayerData(paramString);			// Ignore return value
  lyr.Copy(layerBuffer);
}


// Extract just the units portion of the input string (fifth argument)

G4double CDMS_UIcmdShieldLayer::GetNewUnitValue(const char* paramString) {
  return (Tokenize(paramString) >= 5) ? ValueOf(tokens[4]) : 0.;
}


// Construct input data to string (e.g., for printing)

G4String CDMS_UIcmdShieldLayer::
ConvertToString(const CDMSShieldConstruction::LayerData& lyr,
		const char* unit) {
  G4String line = G4UIcommand::ConvertToString(lyr.Gap) + " "
    + G4UIcommand::ConvertToString(lyr.SideThick) + " "
    + G4UIcommand::ConvertToString(lyr.TopThick) + " "
    + G4UIcommand::ConvertToString(lyr.BottomThick) + " ";
  line += unit;

  return line;
}

G4String CDMS_UIcmdShieldLayer::
ConvertToStringWithBestUnit(const CDMSShieldConstruction::LayerData& lyr) {
  // FIXME:  Don't know how to dynamically get the unit alone from G4BestUnit
  return ConvertToStringWithDefaultUnit(lyr);
}

G4String CDMS_UIcmdShieldLayer::
ConvertToStringWithDefaultUnit(const CDMSShieldConstruction::LayerData& lyr) {
  return ConvertToString(lyr, unitParam->GetDefaultValue());
}


// Assign name to data buffer, not to individual parameters

void CDMS_UIcmdShieldLayer::
SetParameterName(const char* theName, G4bool omittable,
		 G4bool currentAsDefault) {
  layerBuffer.Name = theName;

  // Parameter flags apply "globally" to everything except units
  for (G4int i=0; i<GetParameterEntries()-1; i++) {
    GetParameter(i)->SetOmittable(omittable);
    GetParameter(i)->SetCurrentAsDefault(currentAsDefault);
  }
}

void CDMS_UIcmdShieldLayer::
SetDefaultValue(const CDMSShieldConstruction::LayerData& defVal) {
  layerBuffer = defVal;
}


// Set units configuration -- DefaultUnit is exclusive w.r.t. other two

void CDMS_UIcmdShieldLayer::SetDefaultUnit(const char* defUnit) {
  SetUnitCandidates(CategoryOf(defUnit));

  unitParam->SetDefaultValue(defUnit);
  unitParam->SetOmittable(true);
}

void CDMS_UIcmdShieldLayer::SetUnitCategory(const char* unitCategory) {
  SetUnitCandidates(UnitsList(unitCategory));
}

void CDMS_UIcmdShieldLayer::SetUnitCandidates(const char* candidateList) {
  if (!unitParam->GetDefaultValue().empty()) return;	// User set default

  unitParam->SetParameterCandidates(candidateList);
}


// Copy input string components to buffer for analysis

size_t CDMS_UIcmdShieldLayer::Tokenize(const G4String& parameterList) {
  tokens.clear();
  G4Tokenizer tokenList(parameterList);
  G4String str;
  while( (str = tokenList()) != "" ) {
    tokens.push_back(str);
  }

  return tokens.size();		// Zero can be used as "error condition"
}
