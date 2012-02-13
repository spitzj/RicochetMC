#ifndef RMCVDetectorGeometry_hh
#define RMCVDetectorGeometry_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVDetectorGeometry.hh                              //
//                                                                    //
//  Description: virtual base class for  detector geometries          //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        8 January 2012                                       //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "globals.hh"
#include "G4ThreeVector.hh"
#include <iosfwd>
#include <vector>

class G4LogicalVolume;


class RMCVDetectorGeometry {
public:
    RMCVDetectorGeometry(const G4String& nameString = "none") :
    verboseLevel(0), copyNumber(0), tolerance(0.01*mm),
    name(nameString) {}
    
    virtual ~RMCVDetectorGeometry() {}
    
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
operator<<(std::ostream& os, const RMCVDetectorGeometry& theDet) {
    theDet.PrintParameters(os);
    return os;
}

inline std::ostream& 
operator<<(std::ostream& os, const RMCVDetectorGeometry* theDet) {
    if (theDet) theDet->PrintParameters(os);
    return os;
}

// Templated functions for global reporting

// NOTE:  No bounds checking is possible here, user MUST provide correct N
template <class T>
inline void 
RMCVDetectorGeometry::PrintArrayParameter(std::ostream& os, 
                                          const T par[],
                                          G4int N) const {
    for (G4int i=0; i<N; i++) os << " " << par[i];
}

template <class T>
inline void 
RMCVDetectorGeometry::PrintVectorParameter(std::ostream& os,
                                           const std::vector<T>& par) const {
    for (size_t i=0; i<par.size(); i++) os << " " << par[i];
}

#endif
