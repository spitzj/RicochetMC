// $Id: CDMSGeomMessenger.cc,v 1.17 2011/06/10 06:50:43 dwright Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSGeomMessenger.cc                                 //
//                                                                    //
//  Description: Messenger class to construct complete CDMS models    //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        26 October 2010                                      //
//                                                                    //
//  20101203  M. Kelsey -- Discard "icebox" as separate activity.     //
//  20101208  M. Kelsey -- Add free-standing tower-support detector.  //
//  20110105  M. Kelsey -- Drop "shield" -- move to MultiTower det;   //
//		update geometry verbosity before command processing;  //
//		add support for big "RadiogenicSources" class.        //
//  20110211  M. Kelsey -- Add empty shielding as detector option.    //
//  20110215  M. Kelsey -- Use geometry manager to get components.    //
//  20110421  M. Kelsey -- Add Am-241 source plate as detector option.//
//  20110426  M. Kelsey -- Use Am-241 plate as source, leave old "toy"//
//  20110524  M. Kelsey -- Use new CreateCommand interface.           //
//  20110525  M. Kelsey -- Add GammaSphere source.                    //
////////////////////////////////////////////////////////////////////////

#include "CDMSgeometry/interface/CDMSGeomMessenger.hh"

#include "CDMSgeometry/interface/CDMSGeomConstructor.hh"
#include "CDMSgeometry/interface/CDMSGeometryManager.hh"
#include "CDMSg4base/CDMSVLabConstruction.hh"
#include "CDMSgeometry/labs/NoLab.hh"
#include "CDMSsources/Am241Source.hh"
#include "CDMSsources/CDMSDemoSource.hh"
#include "CDMSsources/CDMSNeutronWall.hh"
#include "CDMSsources/CDMSGammaSphere.hh"
#include "CDMSsources/RadiogenicSources.hh"
#include "CDMSsources/CosmogenicSource.hh"
#include "CDMSsources/SurfaceCosmogenicSource.hh"

#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"


// Constructor and destructor

CDMSGeomMessenger::CDMSGeomMessenger(CDMSGeomConstructor* theSetup)
  : CDMSMessengerBase("/CDMS/", "UI commands to configure CDMS geometry"),
    setup(theSetup), geomManager(CDMSGeometryManager::Instance()) {
  UpdateGeomCmd = CreateCommand<G4UIcmdWithoutParameter>("updateGeom",
	"Update CDMS geometry");
  UpdateGeomCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  UpdateGeomCmd->SetGuidance("if you changed geometrical value(s).");
  UpdateGeomCmd->AvailableForStates(G4State_Idle);

  SelectLabCmd = CreateCommand<G4UIcmdWithAString>("Lab",
	"Define laboratory/cavern for CDMS detector");
  SelectLabCmd->AvailableForStates(G4State_Idle);
  SelectLabCmd->SetCandidates("NoLab SUF SnoLab Soudan DUSEL MITReactor");

  SelectDetCmd = CreateCommand<G4UIcmdWithAString>("Detector",
	"Define active element for CDMS detector");
  SelectDetCmd->SetGuidance(" mini    : Old CDMSMini detector model");
  SelectDetCmd->SetGuidance(" zip     : Single iZip or mZip with housing");
  SelectDetCmd->SetGuidance(" tower   : Single instrumented tower");
  SelectDetCmd->SetGuidance(" cryo    : Empty cryostat (for testing)");
  SelectDetCmd->SetGuidance(" support : Tower support strucutre (for testing)");
  SelectDetCmd->SetGuidance(" shield  : Bare veto shielding (for testing)");
  SelectDetCmd->SetGuidance(" holder  : Holder plate for Am-241 (for testing)");
  SelectDetCmd->SetGuidance(" KGDet   : Strawman multi-tower design");
  SelectDetCmd->SetGuidance(" 100kg   : Configurable multiple-tower detector");
  SelectDetCmd->AvailableForStates(G4State_Idle);
  SelectDetCmd->SetCandidates("mini zip tower cryo support shield holder KGDet 100kg");

  SelectSrcCmd = CreateCommand<G4UIcmdWithAString>("Source",
	"Add background source for CDMS detector");
  SelectSrcCmd->SetGuidance(" beam, demo : Simple electron gun");
  SelectSrcCmd->SetGuidance(" Am241      : Holder plate for Am-241 on ZIP");
  SelectSrcCmd->SetGuidance(" Am241toy   : Dummy disk with fixed gammas");
  SelectSrcCmd->SetGuidance(" sphere     : Spherical radiogenic gammas");
  SelectSrcCmd->SetGuidance(" bkg        : Consolidated radiogenic sources");
  SelectSrcCmd->SetGuidance(" cosmu      : cosmic ray muon source");
  SelectSrcCmd->SetGuidance(" SurfCosmu  : cosmic ray muon source with surface spectrum");
  SelectSrcCmd->SetGuidance("neutronwall : a wall of neutrons");
  SelectSrcCmd->AvailableForStates(G4State_Idle);
  SelectSrcCmd->SetCandidates("Am241 Am241toy sphere beam demo bkg cosmu SurfCosmu neutronwall");
}

CDMSGeomMessenger::~CDMSGeomMessenger() {
  delete UpdateGeomCmd;
  delete SelectLabCmd;
  delete SelectDetCmd;
  delete SelectSrcCmd;
}


// Apply commands

void CDMSGeomMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  CDMSMessengerBase::SetNewValue(cmd,value);	// Check base class first!

  setup->SetVerboseLevel(verboseLevel);		// Update detector verbosity
  
  if (cmd == UpdateGeomCmd) setup->UpdateGeometry();
  if (cmd == SelectLabCmd)  assignLab(value);
  if (cmd == SelectDetCmd)  assignDetector(value);
  if (cmd == SelectSrcCmd)  assignSource(value);
}

void CDMSGeomMessenger::ActionAfterSetVerbose() {
  if (verboseLevel)
    G4cout << " CDMSGeomMessenger::ActionAfterSetVerbose" << G4endl;

  setup->SetVerboseLevel(verboseLevel);
}

void CDMSGeomMessenger::assignLab(const G4String& name) {
  if (verboseLevel) G4cout << " assignLab(" << name << ")" << G4endl;

  CDMSVLabConstruction* theLab = 0;
  if (name == "NoLab") {
    theLab = geomManager->GetNoLab();
  } else if (name == "SnoLab") {
    theLab = geomManager->GetSnoLab();
  } else if (name == "MITReactor") {
    theLab = geomManager->GetMITReactor();
  }

  if (theLab) setup->SetLab(theLab);
  else G4cerr << "Sorry, no " << name << " laboratory yet" << G4endl;
}


void CDMSGeomMessenger::assignDetector(const G4String& name) {
  if (verboseLevel) G4cout << " assignDetector(" << name << ")" << G4endl;

  CDMSVDetectorGeometry* theDet = 0;
  if (name == "mini")    theDet = geomManager->GetMiniDet();
  if (name == "zip")     theDet = geomManager->GetZip();
  if (name == "tower")   theDet = geomManager->GetTower();
  if (name == "cryo")    theDet = geomManager->GetVessel();
  if (name == "support") theDet = geomManager->GetTowerSupport();
  if (name == "shield")  theDet = geomManager->GetShield();
  if (name == "holder")  theDet = geomManager->GetAm241Holder();
  if (name == "KGDet")   theDet = geomManager->GetCKGDet();
  if (name == "100kg")   theDet = geomManager->GetMultiTower();

  if (theDet) setup->SetDetector(theDet);
  else G4cerr << " Sorry, no " << name << " detector yet" << G4endl;
}

void CDMSGeomMessenger::assignSource(const G4String& name) {
  if (verboseLevel) G4cout << " assignSource(" << name << ")" << G4endl;

  CDMSVSourceConstruction* theSrc = 0;
  if (name == "Am241")     theSrc = geomManager->GetAm241Holder();
  if (name == "Am241toy")  theSrc = new Am241Source;
  if (name == "beam")      theSrc = new CDMSDemoSource;
  if (name == "demo")      theSrc = new CDMSDemoSource;
  if (name == "neutronwall") theSrc = new CDMSNeutronWall;
  if (name == "sphere")    theSrc = new CDMSGammaSphere;
  if (name == "bkg")       theSrc = new RadiogenicSources;
  if (name == "cosmu")     theSrc = new CosmogenicSource;
  if (name == "SurfCosmu") theSrc = new SurfaceCosmogenicSource;

  if (theSrc) {
    setup->AddSource(theSrc);				      // For geometry
    G4RunManager::GetRunManager()->SetUserAction(theSrc);     // For generation
  } else G4cerr << " Sorry, no " << name << " source yet" << G4endl;
}
