// $Id: CDMS_UIcmdDoublesListAndUnit.cc,v 1.2 2010/12/23 23:39:39 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMS_UIcmdDoublesListAndUnit.hh                      //     
//  Description: Subclass of DoublesList to support (maybe optional)  //
//		 unit string following input list of values.          //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        22 December 2010                                     //
//                                                                    //
//  20101223  M. Kelsey -- Add local data members to handle units,    //
//		drop use of G4UIparameter.                            //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMS_UIcmdDoublesListAndUnit.hh"
#include "G4UIparameter.hh"
#include "globals.hh"
#include <vector>


// Constructor and destructor

CDMS_UIcmdDoublesListAndUnit::
CDMS_UIcmdDoublesListAndUnit(const char* theCommandPath,
			     G4UImessenger* theMessenger)
  : CDMS_UIcmdDoublesList(theCommandPath, theMessenger) {}

CDMS_UIcmdDoublesListAndUnit::~CDMS_UIcmdDoublesListAndUnit() {}


// Convert user input to raw list of string values (to pick off units)

const std::vector<G4String>& 
CDMS_UIcmdDoublesListAndUnit::GetTokenizedList(const char* st) {
  paramList.clear();

  std::istringstream inputString(st);
  G4String value;
  while (inputString.good() && !inputString.eof()) {
    inputString >> value;
    paramList.push_back(value);
  }

  return paramList;
}

// Convert list of values with unit to std::vector

const std::vector<G4double>& 
CDMS_UIcmdDoublesListAndUnit::GetNewListValue(const char* st) {
  GetNewListRawValue(st);		  // Fill currentList with raw numbers

  G4double unitConv = GetCurrentUnitValue();

  // Rescale input values with unit
  for (size_t i=0; i<currentList.size(); i++) currentList[i] *= unitConv;

  return currentList;
}


// Convert list of values to std::vector, ignoring units at end

const std::vector<G4double>& 
CDMS_UIcmdDoublesListAndUnit::GetNewListRawValue(const char* st) {
  const std::vector<G4String>& tokens = GetTokenizedList(st);

  if (!IsAValidUnit(tokens.back())) {	// No units specified, use base class
    CDMS_UIcmdDoublesList::GetNewListValue(st);
  } else {				// Process list excluding final element
    currentList.clear();
    for (size_t i=0; i<tokens.size()-1; i++)
      currentList.push_back(ConvertToDouble(tokens[i]));
  }

  return currentList;
}


// Get scaling factor for unit at end of string, ignoring list of values

G4double CDMS_UIcmdDoublesListAndUnit::GetNewUnitValue(const char* st) {
  GetTokenizedList(st);
  return GetCurrentUnitValue();
}


// Use last filled set of parameters, or defaults, to return unit scale

G4double CDMS_UIcmdDoublesListAndUnit::GetCurrentUnitValue() const {
  return (IsAValidUnit(paramList.back()) ? ValueOf(paramList.back())
	  : ValueOf(unitDefault));
}


// Convert list plus unit string to a space-delimited string

G4String CDMS_UIcmdDoublesListAndUnit::
ConvertToString(const std::vector<G4double>& list, const char* unit) {
  return CDMS_UIcmdDoublesList::ConvertToString(list)+" "+G4String(unit);
}


// Convert list to a string of digits; Best unit from default or set category

G4String CDMS_UIcmdDoublesListAndUnit::
ConvertToStringWithBestUnit(const std::vector<G4double>& list) {
  // FIXME:  Don't know how to get the unit by itself from G4BestUnit
  return ConvertToStringWithDefaultUnit(list);
}

// Convert list to a string of digits and unit. Best unit from set category

G4String CDMS_UIcmdDoublesListAndUnit::
ConvertToStringWithDefaultUnit(const std::vector<G4double>& list) {
  G4double unitConv = ValueOf(unitDefault);

  std::ostringstream output;
  for (size_t i=0; i<list.size(); i++)
    output << (i>0?" ":"") << list[i]/unitConv;	// Rescale to given units

  output << " " << unitDefault;

  return output.str();
}


// Assign list of candidate unit names, and choose the internal unit as default

void 
CDMS_UIcmdDoublesListAndUnit::SetUnitCandidates(const char* unitCands) {
  unitCandidates = unitCands;
  if (!unitDefault.empty()) return;	// User already set a default

  // Parse all the candidates, and set default to the internal (or nearest)
  // FIXME:  Assumes all candidate strings are valid! (non-zero ValueOf())
  const std::vector<G4String>& cands = GetTokenizedList(unitCands);
  int best = 0;
  G4double bestLog = std::abs(std::log10(ValueOf(cands[0])));
  for (size_t i=1; i<cands.size(); i++) {
    G4double ulog = std::abs(std::log10(ValueOf(cands[i])));
    if (ulog < bestLog) {
      best = i;
      bestLog = ulog;
    }
  }

  unitDefault = cands[best];	// Save internal unit for use as default
}


// Compare input string to registered list of valid units

G4bool CDMS_UIcmdDoublesListAndUnit::IsAValidUnit(const G4String& aUnit) const {
  return (!aUnit.empty() && unitCandidates.find(aUnit) != std::string::npos);
}
  
