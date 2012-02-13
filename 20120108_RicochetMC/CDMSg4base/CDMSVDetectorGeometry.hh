#ifndef CDMSVDetectorGeometry_hh
#define CDMSVDetectorGeometry_hh 1
// $Id: CDMSVDetectorGeometry.hh,v 1.4 2011/06/29 22:24:54 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVDetectorGeometry.hh                             //
//                                                                    //
//  Description: virtual base class for CDMS detector geometries      //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        19 May 2010                                          //
//                                                                    //
//  20101026  M. Kelsey -- Make generic for use with labs and sources //
//  20101122  M. Kelsey -- Add copy-number system for building multis //
//  20101128  M. Kelsey -- Add function to force parameter calcs.     //
//  20101130  M. Kelsey -- Add reporting functions for verbosity.     //
//  20101210  M. Kelsey -- Replace GetMaximumSize() with radius and   //
//		length; Make GetLocation() location, add SetLocation  //
//  20110421  M. Kelsey -- Make list-printing functions templates.    //
//  20110629  M. Kelsey -- Fix constness of ctor arg                  //
//////////////////////////////////////////////////////////////////////// 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>
#include <vector>

class G4LogicalVolume;


class CDMSVDetectorGeometry {
public:
  CDMSVDetectorGeometry(const G4String& nameString = "none") :
    verboseLevel(0), copyNumber(0), tolerance(0.01*mm),
    name(nameString) {}

  virtual ~CDMSVDetectorGeometry() {}

  // Subclasses MUST implement these functions
  virtual G4double GetRadius() const = 0;	  // Radius of maximum extent
  virtual G4double GetLength() const = 0;	  // Z-length of maximum extent

  virtual G4LogicalVolume* BuildGeometry() = 0;	  // Construct for positioning
  virtual G4LogicalVolume* BuildGeometryCopy(G4int copy);

  virtual const G4ThreeVector& GetPosition() const { return position; }

  // Subclasses should override these functions if required
  virtual void FillExtraParameters() {}
  virtual void PrintParameters(std::ostream& os) const;

  virtual void SetVerboseLevel(G4int verbose) { verboseLevel = verbose; }
  virtual void SetPosition(const G4ThreeVector& val) { position = val; }

  const G4String GetName() const {return name;}

protected:
  template <class T>
  void PrintArrayParameter(std::ostream& os, const T par[], G4int N) const;

  template <class T>
  void PrintVectorParameter(std::ostream& os, const std::vector<T>& par) const;

  G4int verboseLevel;
  G4int copyNumber;
  G4ThreeVector position;	// Will default to origin

  G4double tolerance;
  G4String name;
};

// Global reporting classes (inlined) as call-throughs

inline std::ostream& 
operator<<(std::ostream& os, const CDMSVDetectorGeometry& theDet) {
  theDet.PrintParameters(os);
  return os;
}

inline std::ostream& 
operator<<(std::ostream& os, const CDMSVDetectorGeometry* theDet) {
  if (theDet) theDet->PrintParameters(os);
  return os;
}

// Templated functions for global reporting

// NOTE:  No bounds checking is possible here, user MUST provide correct N
template <class T>
inline void 
CDMSVDetectorGeometry::PrintArrayParameter(std::ostream& os, 
					   const T par[],
					   G4int N) const {
  for (G4int i=0; i<N; i++) os << " " << par[i];
}

template <class T>
inline void 
CDMSVDetectorGeometry::PrintVectorParameter(std::ostream& os,
					    const std::vector<T>& par) const {
  for (size_t i=0; i<par.size(); i++) os << " " << par[i];
}

#endif
