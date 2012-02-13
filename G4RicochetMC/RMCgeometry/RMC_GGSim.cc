////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMC_GGSim.cc                                         //
//  Description: main() for RMC general geometry simulation           //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//  Date:        8 January 2012                                       //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "RMCactions/RMCRunAction.hh"
#include "RMCactions/RMCEventAction.hh"
#include "RMCgeometry/interface/RMCGeomConstructor.hh"
#include "RMCsources/RMCPhysicsList.hh"
/*#include "RMCRunAction.hh"
#include "RMCEventAction.hh"
#include "RMCGeomConstructor.hh"
#include "RMCPhysicsList.hh"*/
#include "QGSP_BERT.hh"
#include "QGSP_BERT_HP.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif


int main(int argc,char** argv)
{
    // Choose the Random engine
    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    
    // Construct the default run manager
    G4RunManager* runManager = new G4RunManager;
    
    // Geometry and sources
    RMCGeomConstructor* detector = new RMCGeomConstructor;
    runManager->SetUserInitialization(detector);
    
    // Physics
    G4VUserPhysicsList* physics = new RMCPhysicsList(0);	// verbosity
    runManager->SetUserInitialization(physics);
    
    // Run action
    G4cout << "Setting the RMCRunAction..." << G4endl;
    RMCRunAction* run_action = new RMCRunAction;
    runManager->SetUserAction(run_action);
    
    // Event action
    G4cout << "Setting the RMCEventAction..." << G4endl;
    RMCEventAction* event_action = new RMCEventAction(run_action);
    runManager->SetUserAction(event_action);
    
    //
    // Initialize G4 kernel
    //
    G4cout << "Initializing the G4 kernel..." << G4endl;
    runManager->Initialize();
    
#ifdef G4VIS_USE
    // Initialize visualization
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
#endif
    
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if (argc==2)   // batch mode w macro only
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);    
    }
    if (argc==3) //batch mode w macro and Run ID
    {
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        G4String setRunID = "/run/setRunID ";
        G4String RunID = argv[2];
        UImanager->ApplyCommand(setRunID+RunID);
        UImanager->ApplyCommand(command+fileName);
    }
    else
    {  // interactive mode : define UI session
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
        UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
        if (ui->IsGUI())
            UImanager->ApplyCommand("/control/execute visTutor/gui.mac");     
        ui->SessionStart();
        delete ui;
#endif
    }
    
    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    //                 owned and deleted by the run manager, so they should not
    //                 be deleted in the main() program !
#ifdef G4VIS_USE
    delete visManager;
#endif                
    delete runManager;
    
    return 0;
}
