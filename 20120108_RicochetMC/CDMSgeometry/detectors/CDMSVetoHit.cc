// $Id: CDMSVetoHit.cc,v 1.4 2011/05/04 20:17:22 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSVetoHit.cc                                      //
//  Description: definition of hit information for veto scintillator  //
//               modified from original CDMSmaim version              //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        25 January 2011                                      //
//                                                                    //
//  20110427  M. Kelsey -- Add registration of ZIP coordinates, make  //
//		three-vector accessors return const-refs.             //
//  20110502  M. Kelsey -- Need to invert transform for local coords  //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/detectors/CDMSVetoHit.hh"

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Square.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

G4Allocator<CDMSVetoHit> CDMS_VetoHitAllocator;


CDMSVetoHit::CDMSVetoHit() {}

CDMSVetoHit::CDMSVetoHit(const CDMSVetoHit& right) : G4VHit() {
  *this = right;
}

CDMSVetoHit::~CDMSVetoHit() {}

const CDMSVetoHit& CDMSVetoHit::operator=(const CDMSVetoHit& right)
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

int CDMSVetoHit::operator==(const CDMSVetoHit& /*right*/) const { return 0; }


void CDMSVetoHit::Draw() {
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


void CDMSVetoHit::Print() {
  G4cout << "  trackID: " << trackID << "  Veto" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(postStepPosition,"Length") 
         << "  time: " << G4BestUnit(postStepTime,"Time") << G4endl;
}


void CDMSVetoHit::SetTransform(const G4RotationMatrix* rot,
			       const G4ThreeVector& trans) {
  globalToLocal.SetNetRotation(*rot);
  globalToLocal.SetNetTranslation(trans);
  globalToLocal.Invert();
}
