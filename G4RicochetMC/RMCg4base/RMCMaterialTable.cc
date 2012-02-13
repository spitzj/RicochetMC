////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMaterialTable.cc                                  //     
//  Description: Singleton manager for all material construction      //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               (adapted from Michael Kelsey (SLAC))                 //
//  Date:        22 July 2011                                         //
//                                                                    //
//  All materials used in RMC (regardless of whether they are used    //
//  in any particlar job) must be instantiated at the beginning of    //
//  the job.  NeutronHP depends on the active material database to    //
//  initialize its cross-section tables, and will crash if materials  //
//  are added after that initialization is done (e.g., at geometry    //
//  construction time).                                               //
//////////////////////////////////////////////////////////////////////// 

//#include "RMCg4base/RMCMaterialTable.hh"
#include "RMCMaterialTable.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"


// Force call to constructor, which ensures all materials are created

const RMCMaterialTable RMCMaterialTable::Instance;


// Constructor -- make sure you add new materials here!

RMCMaterialTable::RMCMaterialTable() {
  // Pre-load all NIST materials
  GetMaterial("G4_AIR");
  GetMaterial("G4_Al");
  GetMaterial("G4_C");
  GetMaterial("G4_Cu");
  GetMaterial("G4_Galactic");
  GetMaterial("G4_Ge");
  GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  GetMaterial("G4_POLYCARBONATE");
  GetMaterial("G4_Pb");
  GetMaterial("G4_Si");
  GetMaterial("G4_POLYETHYLENE");
  GetMaterial("G4_Pb");

  // Create materials
  CreateGreenstone();
  CreateNorite();
  CreateMuMetal();
  CreateScintillator();
}

// Greenstone rock is primary component of Soudan cavern
// => Composition taken from M. Tarka, Neutrons from (alpha,n) reactions,
//    Background RMC ebook note (2007)
void RMCMaterialTable::CreateGreenstone() {
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* rock = new G4Material("Greenstone", 2.85*g/cm3, 14);
  rock->AddElement(nist->FindOrBuildElement("H"), 0.018);
  rock->AddElement(nist->FindOrBuildElement("C"), 0.001);
  rock->AddElement(nist->FindOrBuildElement("O"), 0.59311);
  rock->AddElement(nist->FindOrBuildElement("Na"),0.01666);
  rock->AddElement(nist->FindOrBuildElement("Mg"),0.0325);
  rock->AddElement(nist->FindOrBuildElement("Al"),0.060);
  rock->AddElement(nist->FindOrBuildElement("Si"),0.16867);
  rock->AddElement(nist->FindOrBuildElement("P"), 0.0029);
  rock->AddElement(nist->FindOrBuildElement("K"), 0.0026);
  rock->AddElement(nist->FindOrBuildElement("Ca"),0.045);
  rock->AddElement(nist->FindOrBuildElement("Mn"),0.001);
  rock->AddElement(nist->FindOrBuildElement("Fe"),0.0534);
  rock->AddElement(nist->FindOrBuildElement("Ti"),0.00366);
}

// Norite rock is primary component of SNOLAB cavern
void RMCMaterialTable::CreateNorite() {
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* rock = new G4Material("Norite", 2.87*g/cm3, 14);
  rock->AddElement(nist->FindOrBuildElement("H"), 0.0015);
  rock->AddElement(nist->FindOrBuildElement("C"), 0.0004);
  rock->AddElement(nist->FindOrBuildElement("O"), 0.4582);
  rock->AddElement(nist->FindOrBuildElement("Na"),0.0222);
  rock->AddElement(nist->FindOrBuildElement("Mg"),0.0328);
  rock->AddElement(nist->FindOrBuildElement("Al"),0.0892);
  rock->AddElement(nist->FindOrBuildElement("Si"),0.2610);
  rock->AddElement(nist->FindOrBuildElement("P"), 0.0012);
  rock->AddElement(nist->FindOrBuildElement("S"), 0.0020);
  rock->AddElement(nist->FindOrBuildElement("K"), 0.0115);
  rock->AddElement(nist->FindOrBuildElement("Ca"),0.0520);
  rock->AddElement(nist->FindOrBuildElement("Mn"),0.0013);
  rock->AddElement(nist->FindOrBuildElement("Fe"),0.0619);
  rock->AddElement(nist->FindOrBuildElement("Ti"),0.0050);
}

// Mu-metal magnetic shield

void RMCMaterialTable::CreateMuMetal() {
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* mumetal = new G4Material("MuMetal", 8.75*g/cm3, 2);
  mumetal->AddElement(nist->FindOrBuildElement("Fe"), 0.19);
  mumetal->AddElement(nist->FindOrBuildElement("Ni"), 0.81);
}

// Plastic scintillator (C10H11) used in cosmic-ray veto

void RMCMaterialTable::CreateScintillator() {
  G4NistManager* nist = G4NistManager::Instance();

  G4Material* scint = new G4Material("Scintillator", 1.032*g/cm3, 2);
  scint->AddElement(nist->FindOrBuildElement("H"), 11);
  scint->AddElement(nist->FindOrBuildElement("C"), 10);
}


// Load material either from NIST database or from local defintions

G4Material* RMCMaterialTable::GetMaterial(const G4String& mat) {
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(mat);

  // If material not found in NIST database, try for custom definition
  if (!material) material = G4Material::GetMaterial(mat);

  if (!material) 
    G4cerr << " Material " << mat << " not defined" << G4endl;
    
  return material;
}
