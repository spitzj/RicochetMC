////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCSimRoot.hh                                        //
//  Description: class to output G4 hit collections to Root trees     //
//		 DMC units are eV and meters.                                 //
//                                                                    //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        15 July 2010                                         //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#ifndef RMCSimRoot_hh
#define RMCSimRoot_hh 1

#include "RMCgeometry/detectors/RMCZipHit.hh"
#include "RMCgeometry/detectors/RMCVetoHit.hh"

class G4Event;
class TFile;
class TTree;


class RMCSimRoot
{
public:
  RMCSimRoot(G4int verbose=0);
  virtual ~RMCSimRoot();
  
  void SetupFile(const G4String& prefix, G4int runNumber, G4int labCode);
  void ProcessHits(const G4Event*, RMCZipHitsCollection*);
  void ProcessHits(const G4Event*, RMCVetoHitsCollection*);
  void ProcessEventInfo(const G4Event*);
  void WriteAndClose();
  
  void SetVerboseLevel(G4int verbose=1) { verboseLevel = verbose; }
  G4bool GetVerboseLevel() const { return verboseLevel; }
  
private:
  G4String FileName(const G4String& prefix, G4int runNumber, G4int labCode) const;

  void CreateEventTree();
  void CreateVetoTree();
  void CreateAZipTree(G4int i);

  void LoadHitData(const G4Event* evt, const RMCZipHit* ZipHit);
  void LoadHitData(const G4Event* evt, const RMCVetoHit* VetoHit);
  void LoadEventData(const G4Event* evt);

  G4int verboseLevel;
  
  // Leaf variables for zip trees 
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

  // Leaf variables for veto tree 
  G4double vev;
  G4double vdettype;
  G4double vdetnum;
  G4double vts;
  G4double vp;
  G4double vtype;
  G4double ve1;
  G4double vd3;
  G4double vpx3;
  G4double vpy3;
  G4double vpz3;
  G4double vx3;
  G4double vy3;
  G4double vz3;
  G4double vt3;
  G4double vpx1;
  G4double vpy1;
  G4double vpz1;
  G4double vx1;
  G4double vy1;
  G4double vz1;
  G4double vt1;
  G4double vempty;

  // Leaf variables for event tree
  G4double evnum;
  G4double htpevt;
  G4double zppevt;
  G4double primaryPID;
  G4double primaryX;
  G4double primaryY;
  G4double primaryZ;
  G4double primaryPx;
  G4double primaryPy;
  G4double primaryPz;
  G4double primaryT0;
  
  enum {NZIPS = 110};
  
  // TFile and TTrees
  TFile* tfile;
  TTree* eventTree;
  TTree* vetoTree;
  TTree* zipTree[RMCSimRoot::NZIPS];
};

#endif
