////////////////////////////////////////////////////////////////////////
//  File:        RMCGeomMessenger.cc                                  //
//                                                                    //
//  Description: Messenger class to construct complete RMC models     //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               Adapted from: Michael Kelsey (SLAC)                  //
//  Date:        14 January 2012                                      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "RMCgeometry/interface/RMCGeomMessenger.hh"

#include "RMCgeometry/interface/RMCGeomConstructor.hh"
#include "RMCgeometry/interface/RMCGeometryManager.hh"
#include "RMCg4base/RMCVLabConstruction.hh"
#include "RMCgeometry/labs/NoLab.hh"
//#include "RMCsources/RMCDemoSource.hh"
#include "RMCsources/RMCNeutronWall.hh"
//#include "RMCsources/RMCGammaSphere.hh"
//#include "RMCsources/RadiogenicSources.hh"
#include "RMCsources/CosmogenicSource.hh"
//#include "RMCsources/SurfaceCosmogenicSource.hh"

#include "G4RunManager.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"


// Constructor and destructor

RMCGeomMessenger::RMCGeomMessenger(RMCGeomConstructor* theSetup)
: RMCMessengerBase("/RMC/", "UI commands to configure RMC geometry"),
setup(theSetup), geomManager(RMCGeometryManager::Instance()) {
    UpdateGeomCmd = CreateCommand<G4UIcmdWithoutParameter>("updateGeom",
                                                           "Update RMC geometry");
    UpdateGeomCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
    UpdateGeomCmd->SetGuidance("if you changed geometrical value(s).");
    UpdateGeomCmd->AvailableForStates(G4State_Idle);
    
    SelectLabCmd = CreateCommand<G4UIcmdWithAString>("Lab",
                                                     "Define laboratory/cavern for RMC detector");
    SelectLabCmd->AvailableForStates(G4State_Idle);
    SelectLabCmd->SetCandidates("NoLab MITReactor");
    
    SelectDetCmd = CreateCommand<G4UIcmdWithAString>("Detector",
                                                     "Define active element for RMC detector");
    SelectDetCmd->SetGuidance(" zip     : Single iZip or mZip with housing");
    SelectDetCmd->SetGuidance(" tower   : Single instrumented tower");
    SelectDetCmd->SetGuidance(" cryo    : Empty cryostat (for testing)");
    SelectDetCmd->SetGuidance(" support : Tower support strucutre (for testing)");
    SelectDetCmd->SetGuidance(" shield  : Bare veto shielding (for testing)");
    SelectDetCmd->SetGuidance(" experiment  : Modular multi-component detector");
    SelectDetCmd->AvailableForStates(G4State_Idle);
    SelectDetCmd->SetCandidates("zip tower cryo support shield experiment");
    
    SelectSrcCmd = CreateCommand<G4UIcmdWithAString>("Source",
                                                     "Add background source for RMC detector");
    //SelectSrcCmd->SetGuidance(" beam, demo : Simple electron gun");
    //SelectSrcCmd->SetGuidance(" Am241      : Holder plate for Am-241 on ZIP");
    //SelectSrcCmd->SetGuidance(" Am241toy   : Dummy disk with fixed gammas");
    //SelectSrcCmd->SetGuidance(" sphere     : Spherical radiogenic gammas");
    //SelectSrcCmd->SetGuidance(" bkg        : Consolidated radiogenic sources");
    SelectSrcCmd->SetGuidance(" cosmu      : cosmic ray muon source");
    //SelectSrcCmd->SetGuidance(" SurfCosmu  : cosmic ray muon source with surface spectrum");
    SelectSrcCmd->SetGuidance("neutronwall : a wall of neutrons");
    SelectSrcCmd->AvailableForStates(G4State_Idle);
    SelectSrcCmd->SetCandidates("Am241 Am241toy sphere beam demo bkg cosmu SurfCosmu neutronwall");
}

RMCGeomMessenger::~RMCGeomMessenger() {
    delete UpdateGeomCmd;
    delete SelectLabCmd;
    delete SelectDetCmd;
    //delete SelectSrcCmd;
}


// Apply commands

void RMCGeomMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
    RMCMessengerBase::SetNewValue(cmd,value);	// Check base class first!
    
    setup->SetVerboseLevel(verboseLevel);		// Update detector verbosity
    
    if (cmd == UpdateGeomCmd) setup->UpdateGeometry();
    if (cmd == SelectLabCmd)  assignLab(value);
    if (cmd == SelectDetCmd)  assignDetector(value);
}

void RMCGeomMessenger::ActionAfterSetVerbose() {
    if (verboseLevel)
        G4cout << " RMCGeomMessenger::ActionAfterSetVerbose" << G4endl;
    
    setup->SetVerboseLevel(verboseLevel);
}

void RMCGeomMessenger::assignLab(const G4String& name) {
    if (verboseLevel) G4cout << " assignLab(" << name << ")" << G4endl;
    
    RMCVLabConstruction* theLab = 0;
    if (name == "NoLab") {
        theLab = geomManager->GetNoLab();
    } else if (name == "MITReactor") {
        theLab = geomManager->GetMITReactor();
    }
    
    if (theLab) setup->SetLab(theLab);
    else G4cerr << "Sorry, no " << name << " laboratory yet" << G4endl;
}


void RMCGeomMessenger::assignDetector(const G4String& name) {
    if (verboseLevel) G4cout << " assignDetector(" << name << ")" << G4endl;
    
    RMCVDetectorGeometry* theDet = 0;
    if (name == "zip")          theDet = geomManager->GetZip();
    if (name == "tower")        theDet = geomManager->GetTower();
    if (name == "cryo")         theDet = geomManager->GetVessel();
    if (name == "support")      theDet = geomManager->GetTowerSupport();
    if (name == "shield")       theDet = geomManager->GetShield();
    if (name == "experiment")   theDet = geomManager->GetExperiment();
    
    if (theDet) setup->SetDetector(theDet);

    else G4cerr << " Sorry, no " << name << " detector yet" << G4endl;
}

void RMCGeomMessenger::assignSource(const G4String& name) {
 if (verboseLevel) G4cout << " assignSource(" << name << ")" << G4endl;
 
 RMCVSourceConstruction* theSrc = 0;
 //if (name == "Am241")     theSrc = geomManager->GetAm241Holder();
 //if (name == "Am241toy")  theSrc = new Am241Source;
 //if (name == "beam")      theSrc = new RMCDemoSource;
 //if (name == "demo")      theSrc = new RMCDemoSource;
 if (name == "neutronwall") theSrc = new RMCNeutronWall;
 //if (name == "sphere")    theSrc = new RMCGammaSphere;
 //if (name == "bkg")       theSrc = new RadiogenicSources;
 if (name == "cosmu")     theSrc = new CosmogenicSource;
 //if (name == "SurfCosmu") theSrc = new SurfaceCosmogenicSource;
 
 if (theSrc) {
 setup->AddSource(theSrc);				      // For geometry
 G4RunManager::GetRunManager()->SetUserAction(theSrc);     // For generation
 } else G4cerr << " Sorry, no " << name << " source yet" << G4endl;
 }
