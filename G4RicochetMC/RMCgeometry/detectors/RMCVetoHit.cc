////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVetoHit.cc                                      //
//  Description: definition of hit information for veto scintillator  //
//               modified from original RMCmaim version              //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        25 January 2011                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/detectors/RMCVetoHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Square.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

G4Allocator<RMCVetoHit> RMC_VetoHitAllocator;


RMCVetoHit::RMCVetoHit() {}

RMCVetoHit::RMCVetoHit(const RMCVetoHit& right) : G4VHit() {
  *this = right;
}

RMCVetoHit::~RMCVetoHit() {}

const RMCVetoHit& RMCVetoHit::operator=(const RMCVetoHit& right)
{
  if (this != &right) {
    repNum = right.repNum;
    trackID = right.trackID;
    stepNum = right.stepNum;
    parentID = right.parentID;
    pID = right.pID;;
    
    kineticEnergy = right.kineticEnergy;
    edep = right.edep;
    postStepMomentum = right.postStepMomentum;
    postStepPosition = right.postStepPosition;
    postStepTime = right.postStepTime;
    preStepMomentum = right.preStepMomentum;
    preStepPosition = right.preStepPosition;
    preStepTime = right.preStepTime;

    globalToLocal = right.globalToLocal;
  }

  return *this;
}

int RMCVetoHit::operator==(const RMCVetoHit& /*right*/) const { return 0; }


void RMCVetoHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    G4Square square(postStepPosition);
    square.SetScreenSize(10.04);
    square.SetFillStyle(G4Square::filled);
    G4Color color(1.0,0.0,1.0);
    G4VisAttributes attribs(color);
    square.SetVisAttributes(attribs);
    pVVisManager->Draw(square);
  }
}


void RMCVetoHit::Print() {
  G4cout << "  trackID: " << trackID << "  Veto" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(postStepPosition,"Length") 
         << "  time: " << G4BestUnit(postStepTime,"Time") << G4endl;
}


void RMCVetoHit::SetTransform(const G4RotationMatrix* rot,
			       const G4ThreeVector& trans) {
  globalToLocal.SetNetRotation(*rot);
  globalToLocal.SetNetTranslation(trans);
  globalToLocal.Invert();
}
