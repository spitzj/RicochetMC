#ifndef CDMS_TUIcmdWithList_hh
#define CDMS_TUIcmdWithList_hh 1
// $Id: CDMS_TUIcmdWithList.hh,v 1.3 2011/05/24 21:46:57 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMS_TUIcmdWithList.hh                               //     
//  Description: Templated base class to support macro command with   //
//		 fixed (N) or variable-length list of arguments (all  //
//		 of type T)		 		 	      //
//                                                                    //
//  NOTE:        Only double, int, bool or G4String will work with    //
//               GEANT4 code.                                         //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        14 December 2010                                     //
//                                                                    //
//  20110524  M. Kelsey -- Bug fix: ConvertToString not templated!    //
//		Add optional second template argument to set fixed    //
//		length list                                           //
//////////////////////////////////////////////////////////////////////// 

#include "G4UIcommand.hh"
#include "globals.hh"
#include <vector>

class G4UImessenger;


template <typename T, size_t N=0>
class CDMS_TUIcmdWithList : public G4UIcommand {
public:
  CDMS_TUIcmdWithList(const char* theCommandPath, G4UImessenger* theMessenger);
  virtual ~CDMS_TUIcmdWithList();

  // Convert string which represents list of values to std::vector
  const std::vector<T>& GetNewListValue(const char* st);

  // Write input list of values to space-delimited string
  G4String ConvertToString(const std::vector<T>& list);

  // Set the parameter name for the whole list, treated as a string.
  //  If "omittable" is set as true, the user of this command can ommit
  // the value(s) when he/she applies the command. If "omittable" is false,
  // the user must supply all three values.
  //  "currentAsDefault" flag is valid only if "omittable" is true. If this
  // flag is true, the current values are used as the default values when the 
  // user ommit some of the parameters. If this flag is false, the values
  // given by the SetDefaultValue() method are used. 
  void SetParameterName(const char* theName, G4bool omittable,
			G4bool currentAsDefault=false);

  // Set the default values of the parameters. These default values are used
  // when the user of this command omits the parameter values, and
  // "ommitable" is true and "currentAsDefault" is false.
  void SetDefaultValue(const std::vector<T>& defVal);

protected:
  std::vector<T> defaultList;
  std::vector<T> currentList;
};

// Templated functions must be implemented in each compilation unit
#include "CDMSg4base/CDMS_TUIcmdWithList.icc"

#endif	/* CDMS_TUIcmdWithList_hh */
