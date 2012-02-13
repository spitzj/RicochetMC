////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCZipSD.cc                                         //
//  Description: sensitive detector class for RMCMini                //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        15 July 2010                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/detectors/RMCZipSD.hh"
#include "G4UnitsTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4VTouchable.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


RMCZipSD::RMCZipSD(const G4String& SDname)
 :G4VSensitiveDetector(SDname)
{
  G4String HCname = "zipHitsCollection";
  collectionName.insert(HCname);
}

RMCZipSD::~RMCZipSD() {}

void RMCZipSD::Initialize(G4HCofThisEvent* HCE)
{
  zipHitsCollection = new RMCZipHitsCollection(SensitiveDetectorName,collectionName[0]);
  G4int HCID = -1;
  if(HCID < 0) HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection(HCID, zipHitsCollection);
}


G4bool RMCZipSD::ProcessHits(G4Step* aStep, G4TouchableHistory* /*ROhist*/)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* postPoint = aStep->GetPostStepPoint();

  // Make sure step actually deposited energy
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false;
	
  RMCZipHit* zipHit = new RMCZipHit();

  // Transformation from global to local coordinates, needed by DMC
  zipHit->SetTransform(prePoint->GetTouchable()->GetRotation(),
		       prePoint->GetTouchable()->GetTranslation());

  zipHit->SetReplicaNum(prePoint->GetTouchable()->GetReplicaNumber() + 1);
  zipHit->SetTrackID(aStep->GetTrack()->GetTrackID() ); 
  zipHit->SetStepNum(aStep->GetTrack()->GetCurrentStepNumber() );
  zipHit->SetParentID(aStep->GetTrack()->GetParentID() );

  G4double pid;
  if (aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus") {
    pid = (aStep->GetTrack()->GetDefinition()->GetPDGCharge()) +
          (1000 * aStep->GetTrack()->GetDefinition()->GetBaryonNumber());
  } else {
    pid = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();
  }
  zipHit->SetPID(pid);

  zipHit->SetPreStepKE(prePoint->GetKineticEnergy() );
  zipHit->SetEdep(edep);

  zipHit->SetPostStepMomentum(postPoint->GetMomentum() );
  zipHit->SetPostStepPosition(postPoint->GetPosition() );
  zipHit->SetPostStepTime(postPoint->GetGlobalTime() );

  zipHit->SetPreStepMomentum(prePoint->GetMomentum() );
  zipHit->SetPreStepPosition(prePoint->GetPosition() );
  zipHit->SetPreStepTime(prePoint->GetGlobalTime() );

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
	
  zipHitsCollection->insert(zipHit);
  //  zipHit->Print();
  zipHit->Draw();
	
  return true;
}


void RMCZipSD::EndOfEvent(G4HCofThisEvent* /*HCE*/)
{
  if (verboseLevel>0) {
    G4int NbHits = zipHitsCollection->entries();     
    G4cout << "\n--------> Hits Collection: in this event there are " << NbHits
	   << " hits in the Zip: " << G4endl;

    G4double TotEnergyDep = 0;
	
    for (G4int ii=0; ii<NbHits; ii++) {
			//(*zipHitsCollection)[ii]->Print();
			TotEnergyDep += (*zipHitsCollection)[ii]->GetEdep();
    }

    G4cout << "--------> Total deposited energy in " 
	   << SensitiveDetectorName << " : " 
	   << std::setw(5) << G4BestUnit(TotEnergyDep,"Energy")
	   << "\n" << G4endl;
  }
}
