#ifndef CDMS_UIcmdDoublesListAndUnit_hh
#define CDMS_UIcmdDoublesListAndUnit_hh 1
// $Id: CDMS_UIcmdDoublesListAndUnit.hh,v 1.2 2010/12/23 23:39:39 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMS_UIcmdDoublesListAndUnit.hh                      //     
//  Description: Subclass of DoublesList to support (maybe optional)  //
//		 unit string following input list of values.          //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        22 December 2010                                     //
//                                                                    //
//  20101223  M. Kelsey -- Add local data members to handle units     //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMS_UIcmdDoublesList.hh"
#include "globals.hh"
#include <vector>

class G4UImessenger;


class CDMS_UIcmdDoublesListAndUnit : public CDMS_UIcmdDoublesList {
public:
  CDMS_UIcmdDoublesListAndUnit(const char* theCommandPath,
			       G4UImessenger* theMessenger);

  virtual ~CDMS_UIcmdDoublesListAndUnit();

  // Convert list of values with unit to std::vector
  const std::vector<G4double>& GetNewListValue(const char* st);

  // Convert list of values to std::vector, ignoring units at end
  const std::vector<G4double>& GetNewListRawValue(const char* st);

  // Get scaling factor for unit at end of string, ignoring list of values
  G4double GetNewUnitValue(const char* st);

  // Get scaling factor for last input string, or from default
  G4double GetCurrentUnitValue() const;

  // Convert list plus unit string to a string
  G4String ConvertToString(const std::vector<G4double>& list, const char* unit);

  // Convert list to a string of digits and unit. Best unit is
  // chosen from the unit category of default unit (in case SetDefaultUnit()
  // is defined) or category defined by SetUnitCategory().
  G4String ConvertToStringWithBestUnit(const std::vector<G4double>& list);

  // Convert list to a string of digits and unit. Best unit is
  // chosen from the category defined by SetUnitCategory() in case default
  // unit is not defined.
  G4String ConvertToStringWithDefaultUnit(const std::vector<G4double>& list);

  // These three methods must be used alternatively.
  // The user cannot ommit the unit as the last parameter of the command if
  // SetUnitCategory() or SetUnitCandidates() is used.

  // Available categories can be found in G4SystemOfUnits.hh.
  // Only the units categorized in the given category will be accepted.
  void SetUnitCategory(const char* unitCat) {
    SetUnitCandidates(UnitsList(unitCat));
  }

  // Units listed in the argument of this method must be separated by space(s).
  // Only the units listed in the candidate list will be accepted.
  void SetUnitCandidates(const char* unitCands);

  // Define both default unit and category.  User may omit unit string.
  void SetDefaultUnit(const char* defUnit) { unitDefault = defUnit; }

protected:
  G4bool IsAValidUnit(const G4String& aUnit) const;

  const std::vector<G4String>& GetTokenizedList(const char* st);

  std::vector<G4String> paramList;	// Buffer to store all user arguments

  G4String unitCategory;		// Local buffers instead of G4UIparam
  G4String unitCandidates;
  G4String unitDefault;
};

#endif	/* CDMS_UIcmdDoublesListAndUnit_hh */
