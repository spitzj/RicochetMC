// $Id: CDMSZipHit.cc,v 1.6 2011/05/04 20:17:22 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSZipHit.cc                                        //
//  Description: definition of hit information for zip detector,      //
//               modified from original CDMSmini version              //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        7 July 2010                                          //
//                                                                    //
//  20110427  M. Kelsey -- Add coordinate transformation, modify      //
//		position accessors to return local coordinates.       //
//  20110502  M. Kelsey -- Add hardwired Z-coordinate offset.         //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/detectors/CDMSZipHit.hh"

#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSgeometry/detectors/CDMSZipConstruction.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"

G4Allocator<CDMSZipHit> CDMS_ZipHitAllocator;


CDMSZipHit::CDMSZipHit() {}

CDMSZipHit::CDMSZipHit(const CDMSZipHit& right) : G4VHit() {
  *this = right;
}

CDMSZipHit::~CDMSZipHit() {}

const CDMSZipHit& CDMSZipHit::operator=(const CDMSZipHit& right) {
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


int CDMSZipHit::operator==(const CDMSZipHit& /*right*/) const { return 0; }


void CDMSZipHit::Draw() {
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


void CDMSZipHit::Print() {
  G4cout << "  trackID: " << trackID << "  Zip" 
	 << "  energy deposit: " << G4BestUnit(edep,"Energy")
	 << "  position: " << G4BestUnit(postStepPosition,"Length") 
         << "  time: " << G4BestUnit(postStepTime,"Time") << G4endl;
}


void CDMSZipHit::SetTransform(const G4RotationMatrix* rot,
			      const G4ThreeVector& trans) {
  globalToLocal.SetNetRotation(*rot);
  globalToLocal.SetNetTranslation(trans);
  globalToLocal.Invert();

  // Need to query ZIP directly, in order to apply additional Z shift
  CDMSZipConstruction* theZip = CDMSGeometryManager::Instance()->GetZip();
  if (!theZip) theZip = new CDMSZipConstruction();	// Failsafe, but wrong

  G4ThreeVector zipZoffset(0., 0., theZip->GetZipThick()/2.);
  globalToLocal += zipZoffset;
}
