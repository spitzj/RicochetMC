////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCVetoSD.cc                                        //
//  Description: sensitive detector class for veto shield             //
//               scintillator                                         //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        25 January 2011                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/detectors/RMCVetoSD.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VTouchable.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


RMCVetoSD::RMCVetoSD(const G4String& SDname)
 :G4VSensitiveDetector(SDname)
{
  G4String HCname = "vetoHitsCollection";
  collectionName.insert(HCname);
}

RMCVetoSD::~RMCVetoSD() {}

void RMCVetoSD::Initialize(G4HCofThisEvent* HCE)
{
  vetoHitsCollection =
    new RMCVetoHitsCollection(SensitiveDetectorName,collectionName[0]);
  G4int HCID = -1;
  if(HCID < 0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection(HCID, vetoHitsCollection);
}


G4bool RMCVetoSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.) return false;
	
  RMCVetoHit* vetoHit = new RMCVetoHit();

  // Transformation from global to local coordinates, needed by DMC
  /***** FOR NOW, DON'T DO THE TRANSFORMATION
  vetoHit->SetTransform(postPoint->GetTouchable()->GetRotation(),
		       postPoint->GetTouchable()->GetTranslation());
  *****/

  vetoHit->SetReplicaNum(prePoint->GetTouchable()->GetReplicaNumber() + 1);
  vetoHit->SetTrackID(aStep->GetTrack()->GetTrackID() ); 
  vetoHit->SetStepNum(aStep->GetTrack()->GetCurrentStepNumber() );
  vetoHit->SetParentID(aStep->GetTrack()->GetParentID() );

  G4double pid;
  if (aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus") {
    pid = (aStep->GetTrack()->GetDefinition()->GetPDGCharge()) +
          (1000 * aStep->GetTrack()->GetDefinition()->GetBaryonNumber());
  } else {
    pid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  }
  vetoHit->SetPID(pid);

  vetoHit->SetPreStepKE(prePoint->GetKineticEnergy() );
  vetoHit->SetEdep(edep);

  vetoHit->SetPostStepMomentum(postPoint->GetMomentum() );
  vetoHit->SetPostStepPosition(postPoint->GetPosition() );
  vetoHit->SetPostStepTime(postPoint->GetGlobalTime() );

  vetoHit->SetPreStepMomentum(prePoint->GetMomentum() );
  vetoHit->SetPreStepPosition(prePoint->GetPosition() );
  vetoHit->SetPreStepTime(prePoint->GetGlobalTime() );

  /*
   1. EV (starts with 1)
   2. DT (detector #)
   3. TS (track and step: ttttsssss)
   4. PA (parent track  : tttt00000)
   5. TY (type. -ve =>  : -zzzaaa)
   6. E1 (KE of this track@1)
   7. D3 (energy deposition@3)
   8. PX3 (momentum@3)
   9. PY3
  10. PZ3
  11. X3 (position@3)
  12. Y3
  13. Z3 
  14. PX1 (momentum@1)
  15. PY1
  16. PZ1
  17. X1 (position@1)
  18. Y1
  19. Z1 
 */
	
  vetoHitsCollection->insert(vetoHit);
  //  vetoHit->Print();
  vetoHit->Draw();
	
  return true;
}


void RMCVetoSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
  if (verboseLevel>0) {
    G4int NbHits = vetoHitsCollection->entries();     
    G4cout << "\n--------> Hits Collection: in this event there are " << NbHits
	   << " hits in the Veto: " << G4endl;

    G4double TotEnergyDep = 0;
	
    for (G4int ii=0; ii<NbHits; ii++) {
			//(*vetoHitsCollection)[ii]->Print();
			TotEnergyDep += (*vetoHitsCollection)[ii]->GetEdep();
    }

    G4cout << "--------> Total deposited energy in " 
	   << SensitiveDetectorName << " : " 
	   << std::setw(5) << G4BestUnit(TotEnergyDep,"Energy")
	   << "\n" << G4endl;
  }
}
