////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVetoHit.hh                                        //
//  Description: hit definition class for veto scintillator,          //
//               modified from original RMCmain version               //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        25 January 2011                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#ifndef RMCVetoHit_h
#define RMCVetoHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"


class RMCVetoHit : public G4VHit {
public:
  RMCVetoHit();
  RMCVetoHit(const RMCVetoHit&);
  
  ~RMCVetoHit();
  
  const RMCVetoHit& operator=(const RMCVetoHit&);
  int operator == (const RMCVetoHit&) const;
  
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  
  void Draw();
  void Print();
  
  void SetTransform(const G4RotationMatrix* rot, const G4ThreeVector& trans);
  
  void SetReplicaNum(G4int Nrep) {repNum = Nrep;}
  void SetTrackID(G4int tID) {trackID = tID;}
  void SetStepNum(G4int nStep) {stepNum = nStep;}
  void SetParentID(G4int parID) {parentID = parID;}
  void SetPID(G4double pid) {pID = pid;}
  void SetPreStepKE(G4double ke) {kineticEnergy = ke;}
  void SetEdep(G4double de) {edep = de;}
  
  void SetPostStepMomentum(const G4ThreeVector& postMom) {postStepMomentum = postMom;}
  void SetPostStepPosition(const G4ThreeVector& postPos) {postStepPosition = postPos;}
  void SetPostStepTime(G4double postTime) {postStepTime = postTime;}
  void SetPreStepMomentum(const G4ThreeVector& preMom) {preStepMomentum = preMom;}
  void SetPreStepPosition(const G4ThreeVector& prePos) {preStepPosition = prePos;}
  void SetPreStepTime(G4double preTime) {preStepTime = preTime;}
  
  G4int GetReplicaNum() const {return repNum;}
  G4int GetTrackID() const {return trackID;}
  G4int GetStepNum() const {return stepNum;}
  G4int GetParentID() const {return parentID;}
  G4double GetPID() const {return pID;}
  
  G4double GetPreStepKE() const {return kineticEnergy;}
  G4double GetEdep() const {return edep;}
  const G4ThreeVector& GetPostStepMomentum() const {return postStepMomentum;}
  G4double GetPostStepTime() const {return postStepTime;}
  const G4ThreeVector& GetPreStepMomentum() const {return preStepMomentum;}
  G4double GetPreStepTime() const {return preStepTime;}
  
  // Return hit locations in local coordinates for this panel
  G4ThreeVector GetPostStepPosition() const {
    return globalToLocal.TransformPoint(postStepPosition);
  }

  G4ThreeVector GetPreStepPosition() const {
    return globalToLocal.TransformPoint(preStepPosition);
  }
  
private:
  G4int repNum;
  G4int trackID;
  G4int stepNum;
  G4int parentID;
  G4double pID;
  
  G4double kineticEnergy;
  G4double edep;
  G4ThreeVector postStepMomentum;
  G4ThreeVector postStepPosition;
  G4double postStepTime;
  G4ThreeVector preStepMomentum;
  G4ThreeVector preStepPosition;
  G4double preStepTime;
  
  G4AffineTransform globalToLocal;
};


// Data and memory management

typedef G4THitsCollection<RMCVetoHit> RMCVetoHitsCollection;
extern G4Allocator<RMCVetoHit> RMC_VetoHitAllocator;

inline void* RMCVetoHit::operator new(size_t) {
  void* aHit;
  aHit = (void*) RMC_VetoHitAllocator.MallocSingle();
  return aHit;
}

inline void RMCVetoHit::operator delete(void* aHit) {
  RMC_VetoHitAllocator.FreeSingle((RMCVetoHit*) aHit);
}

#endif	/* RMCVetoHit_h */
