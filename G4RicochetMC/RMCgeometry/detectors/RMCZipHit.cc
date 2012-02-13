////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCZipHit.cc                                        //
//  Description: definition of hit information for zip detector,      //
//               modified from original RMCmini version              //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        7 July 2010                                          //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/detectors/RMCZipHit.hh"

#include "RMCgeometry/interface/RMCGeometryManager.hh"
#include "RMCgeometry/detectors/RMCZipConstruction.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

G4Allocator<RMCZipHit> RMC_ZipHitAllocator;


RMCZipHit::RMCZipHit() {}

RMCZipHit::RMCZipHit(const RMCZipHit& right) : G4VHit() {
  *this = right;
}

RMCZipHit::~RMCZipHit() {}

const RMCZipHit& RMCZipHit::operator=(const RMCZipHit& right) {
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


int RMCZipHit::operator==(const RMCZipHit& /*right*/) const { return 0; }


void RMCZipHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    G4Circle circle(postStepPosition);
    circle.SetScreenSize(10.04);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0/255.,100/255.,128/255.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}


void RMCZipHit::Print() {
  G4cout << "  trackID: " << trackID << "  Zip" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(postStepPosition,"Length") 
         << "  time: " << G4BestUnit(postStepTime,"Time") << G4endl;
}


void RMCZipHit::SetTransform(const G4RotationMatrix* rot,
			      const G4ThreeVector& trans) {
  globalToLocal.SetNetRotation(*rot);
  globalToLocal.SetNetTranslation(trans);
  globalToLocal.Invert();

  // Need to query ZIP directly, in order to apply additional Z shift
  RMCZipConstruction* theZip = RMCGeometryManager::Instance()->GetZip();
  if (!theZip) theZip = new RMCZipConstruction();	// Failsafe, but wrong

  G4ThreeVector zipZoffset(0., 0., theZip->GetZipThick()/2.);
  globalToLocal += zipZoffset;
}
