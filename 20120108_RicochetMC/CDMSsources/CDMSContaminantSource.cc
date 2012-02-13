////////////////////////////////////////////////////////////////////////
// $Id: CDMSContaminantSource.cc,v 1.4 2011/07/08 23:26:34 kelsey Exp $
//  File:        CDMSContaminantSource.cc                             //
//  Description: Base class (fully functional) for sources which are  //
//		 contaminants on/in detector components               //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        5 July 2011                                          //
//                                                                    //
//  20110707  Drop list of contaminants; use MultiGenerator instead.  //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/CDMSContaminantSource.hh"
/* #include "CDMSsources/CDMSContaminantMessenger.hh" */
#include "G4Event.hh"
#include "G4GenericIon.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4Point3D.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RadioactiveDecay.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Vector3D.hh"
#include "Randomize.hh"


// Constructors and destructor

CDMSContaminantSource::CDMSContaminantSource(const G4String& name,
					     const G4String& volume)
  : CDMSVSourceConstruction(name), volumeName(volume), sources(name),
    numberOfParticles(1), halfAngle(90.*deg), messenger(0) {}

CDMSContaminantSource::~CDMSContaminantSource() {
  //*** delete messenger;
}


// Create event (single or multiple particles)

void CDMSContaminantSource::GeneratePrimaries(G4Event* evt) {
  if (verboseLevel>1) G4cout << GetName() << "::GeneratePrimaries" << G4endl;

  if (numberOfParticles<1) return;		// Avoid unnecessary work

  if (evt->GetEventID() < 1) FindVolumes();	// Get updated geometry
  if (lVolumes.empty() && pVolumes.empty()) return;  // Avoid unnecessary work

  for (G4int i=0; i<numberOfParticles; i++) {
    GeneratePrimary(evt);
  }
}

void CDMSContaminantSource::GeneratePrimary(G4Event* evt) {
  if (verboseLevel>1) G4cout << GetName() << "::GeneratePrimary" << G4endl;

  G4VPhysicalVolume* detvol = ChooseVolume();
  if (!detvol) return;				// No volumes, no work to do

  G4ThreeVector srcLoc, srcDir;			// Buffers for selection
  G4double ekin;
  G4ParticleDefinition* srcType;

  ChoosePointAndDir(detvol, srcLoc, srcDir);
  sources.SetDirection(srcDir, halfAngle);
  sources.shoot(ekin, srcDir, srcType);		// Get generated particle

#ifdef CDMS_RDM_COLLIMATION
  // Configure collimation for radioactive decay products
  if (dynamic_cast<G4GenericIon*>(srcType)) {
    static G4RadioactiveDecay* rdmProcess = findRadioactiveDecayProcess();
    if (rdmProcess) {
      rdmProcess->SetDecayCollimation(srcDir, halfAngle);
      if (verboseLevel>2)
	G4cout << " decay daughters collimated to "
	       << generator.GetDirection() << G4endl;
    }
  }
#endif	/* CDMS_RDM_COLLIMATION */

  particleGun->SetParticlePosition(srcLoc);
  particleGun->SetParticleMomentumDirection(srcDir);
  particleGun->SetParticleEnergy(ekin);

  // Add gamma production or radioactive decay to event
  particleGun->GeneratePrimaryVertex(evt);  
}


// Select random point on surface for source, directed locally outward

void 
CDMSContaminantSource::ChoosePointAndDir(G4VPhysicalVolume* detvol,
					 G4ThreeVector& location,
					 G4ThreeVector& direction) const {
  if (verboseLevel>1) G4cout << GetName() << "::ChoosePointAndDir" << G4endl;

  G4VSolid* shape = detvol->GetLogicalVolume()->GetSolid();
  G4ThreeVector localSource = shape->GetPointOnSurface();
  G4ThreeVector localOut = shape->SurfaceNormal(localSource);

  // Convert local to global coordinates -- points and vectors are different!
  G4Transform3D toGlobal = GetWorldFrame(detvol);
  location = toGlobal * G4Point3D(localSource);
  direction = toGlobal * G4Vector3D(localOut);

  if (verboseLevel>2)
    G4cout << " source at " << location << " toward " << direction << G4endl;
}


// Select physical volume from list to use with event

G4VPhysicalVolume* CDMSContaminantSource::ChooseVolume() const {
  if (pVolumes.empty()) return 0;		// Empty list, do nothing

  if (verboseLevel>1) G4cout << GetName() << "::ChooseVolume" << G4endl;

  G4int ivol = pVolumes.size() * G4UniformRand();
  if (verboseLevel>2) 
    G4cout << " selected " << ivol << " of " << pVolumes.size() << G4endl;

  return pVolumes[ivol];
}


// Collect named volumes from currently active geometry

void CDMSContaminantSource::FindVolumes() {
  FindLogicalVolumes();
  FindPhysicalVolumes();
}

void CDMSContaminantSource::FindLogicalVolumes() {
  if (verboseLevel) G4cout << GetName() << "::FindLogicalVolumes" << G4endl;

  lVolumes.clear();

  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  if (!lvStore) return;

  G4LogicalVolume* vol = 0;
  for (size_t iv=0; iv<lvStore->size(); iv++) {
    vol = (*lvStore)[iv];
    if (vol && vol->GetName() == volumeName) lVolumes.push_back(vol);
  }
}

void CDMSContaminantSource::FindPhysicalVolumes() {
  if (verboseLevel) G4cout << GetName() << "::FindPhysicalVolumes" << G4endl;

  pVolumes.clear();

  G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  if (!pvStore) return;

  G4VPhysicalVolume* vol = 0;
  for (size_t iv=0; iv<pvStore->size(); iv++) {
    vol = (*pvStore)[iv];
    if (vol && vol->GetName() == volumeName) pVolumes.push_back(vol);
  }
}


// Determine transformation from local to global coordinates

G4Transform3D 
CDMSContaminantSource::GetWorldFrame(G4VPhysicalVolume* vol) const {
  if (verboseLevel>1) G4cout << GetName() << "::GetWorldFrame" << G4endl;

  // Get local positioning of volume in mother
  G4Transform3D localToMother(*(vol->GetFrameRotation()),
			      vol->GetFrameTranslation());

  // FIXME:  No way to traverse up the hierarchy; mother is a Logical!
  G4LogicalVolume* mother = vol->GetMotherLogical();

  return localToMother;
}


// Search process list for RadioactiveDecay (there should be only one)

G4RadioactiveDecay* 
CDMSContaminantSource::findRadioactiveDecayProcess() const {
  if (verboseLevel > 1)
    G4cout << "CDMSGammaSphere::findRadioactiveDecayProcess" << G4endl;

  G4RadioactiveDecay* rdm = 0;		// Buffer for return value

  static G4ParticleDefinition* allIons = G4GenericIon::Definition();

  G4ProcessManager* procMan = allIons->GetProcessManager();
  if (procMan) {
    G4ProcessVector* ionProcs = procMan->GetProcessList();
    if (ionProcs) {
      for (G4int ip=0; ip<ionProcs->size(); ip++) {
	rdm = dynamic_cast<G4RadioactiveDecay*>((*ionProcs)[ip]);
	if (rdm) break;			// Found the right process
      }
    }
  }

  if (rdm && verboseLevel > 1)
    G4cout << " found RadioactiveDecay associated with GenericIon" << G4endl;

  return rdm;		// Will be null if conditionals above failed
}
