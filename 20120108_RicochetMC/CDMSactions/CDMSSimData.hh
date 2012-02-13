// $Id: CDMSSimData.hh,v 1.2 2011/07/01 03:39:30 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSSimData.hh                                       //
//  Description: class to output G4 hit collections to text file      //
//		 DMC units are eV and meters.                         //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        20 April 2011                                        //
//                                                                    //
//  20110420  M. Kelsey -- Adapted from CDMSSimRoot                   //
//  20110630  M. Kelsey -- Move filename construction here from RunAc //
//////////////////////////////////////////////////////////////////////// 

#ifndef CDMSSimData_hh
#define CDMSSimData_hh 1

#include "CDMSgeometry/detectors/CDMSZipHit.hh"
#include "CDMSgeometry/detectors/CDMSVetoHit.hh"
#include <iosfwd>

class G4Event;

class CDMSSimData
{
public:
  CDMSSimData(G4int verbose=0);
  virtual ~CDMSSimData();
  
  void SetupFile(const G4String& prefix, G4int runNumber, G4int labCode=0);
  void ProcessHits(const G4Event*, CDMSZipHitsCollection*);
  void ProcessHits(const G4Event*, CDMSVetoHitsCollection*);
  void WriteAndClose();
  
  void SetVerboseLevel(G4int verbose=1) { verboseLevel = verbose; }
  G4bool GetVerboseLevel() const { return verboseLevel; }
  
private:
  G4String FileName(const G4String& prefix, G4int runNumber, G4int labCode=0) const;

  void LoadHitData(const G4Event* evt, const CDMSZipHit* ZipHit);
  void LoadHitData(const G4Event* evt, const CDMSVetoHit* VetoHit);
  void WriteHitData();

  G4int verboseLevel;

  std::ostream* outputFile;

  // Leaf variables for hit data (zip and veto)
  G4double ev;
  G4double dettype;
  G4double detnum;
  G4double ts;
  G4double p;
  G4double type;
  G4double e1;
  G4double d3;
  G4double px3;
  G4double py3;
  G4double pz3;
  G4double x3;
  G4double y3;
  G4double z3;
  G4double t3;
  G4double px1;
  G4double py1;
  G4double pz1;
  G4double x1;
  G4double y1;
  G4double z1;
  G4double t1;
  G4double empty;

  // Leaf variables for event data
  G4double evnum;
  G4double htpevt;
  G4double zppevt;
};

#endif
