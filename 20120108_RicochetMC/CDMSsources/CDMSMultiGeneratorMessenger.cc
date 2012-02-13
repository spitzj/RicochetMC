// $Id: CDMSMultiGeneratorMessenger.cc,v 1.3 2011/07/21 21:22:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//  File:        CDMSMultiGeneratorMessenger.cc                       //
//  Description: messenger class for multiple-particle generator      //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        6 July 2011                                          //
//                                                                    //
//  20110707  M. Kelsey -- Extract source additions from GammaSphere  //
//  20110721  M. Kelsey -- Remove "SetCandidates" from AddSpectrumCmd //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/CDMSMultiGeneratorMessenger.hh"
#include "CDMSsources/CDMSMultiGenerator.hh"
#include "CDMSsources/Ac228Lines.hh"
#include "CDMSsources/Am241Lines.hh"
#include "CDMSsources/Bi212Lines.hh"
#include "CDMSsources/Bi214Lines.hh"
#include "CDMSsources/Co60Lines.hh"
#include "CDMSsources/K40Lines.hh"
#include "CDMSsources/Pa234Lines.hh"
#include "CDMSsources/Pb212Lines.hh"
#include "CDMSsources/Pb214Lines.hh"
#include "CDMSsources/Tl208Lines.hh"
#include "CDMSsources/UniformLines.hh"
#include "globals.hh"
#include "G4Tokenizer.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"


// Constructor and destructor -- note dynamic creation of path

CDMSMultiGeneratorMessenger::
CDMSMultiGeneratorMessenger(CDMSMultiGenerator* theGenerator)
  : CDMSMessengerBase(("/CDMS/"+theGenerator->GetName()+"/").data(),
		      "Weighted-composition particle generator"),
    generator(theGenerator), AddSpectrumCmd(0), AddIsotopeCmd(0), PrintCmd(0) {
  AddSpectrumCmd = CreateCommand<G4UIcmdWithAString>("AddSpectrum",
	"Add hardcoded gamma spectrum to source");
  AddSpectrumCmd->SetParameterName("source",false);
  AddSpectrumCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  AddSpectrumCmd->SetGuidance("<src> [weight]: Optional weight/fraction for each spectrum");
  AddSpectrumCmd->SetGuidance(" Ac228 : Gammas from actinium-228");
  AddSpectrumCmd->SetGuidance(" Am241 : Gammas and X-rays from Am-241 chain");
  AddSpectrumCmd->SetGuidance(" Bi212 : Gammas from bismuth-212");
  AddSpectrumCmd->SetGuidance(" Bi214 : Gammas from bismuth-214");
  AddSpectrumCmd->SetGuidance(" Co60  : Gammas and X-rays from cobalt-60");
  AddSpectrumCmd->SetGuidance(" K40   : Gammas and X-rays from potassium-40");
  AddSpectrumCmd->SetGuidance(" Pa234 : Gammas from protactinium-234");
  AddSpectrumCmd->SetGuidance(" Pb212 : Gammas from lead-212");
  AddSpectrumCmd->SetGuidance(" Pb214 : Gammas from lead-214");
  AddSpectrumCmd->SetGuidance(" Tl208 : Gammas from thallium-208");
  AddSpectrumCmd->SetGuidance(" Uniform : Flat spectrum from 1-2500 keV");

  AddIsotopeCmd = CreateCommand<G4UIcmdWithAString>("AddIsotope",
	"Add (Z,A) radionuclide to source");
  AddSpectrumCmd->SetGuidance("<Z> <A> [weight]: Optional weight/fraction for isotope");
  AddIsotopeCmd->SetParameterName("source",false);
  AddIsotopeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  PrintCmd = CreateCommand<G4UIcmdWithoutParameter>("PrintSources",
	"Report current composition");
  PrintCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

CDMSMultiGeneratorMessenger::~CDMSMultiGeneratorMessenger() {
  delete AddSpectrumCmd;
  delete AddIsotopeCmd;
  delete PrintCmd;
}


// Apply user command inputs

void CDMSMultiGeneratorMessenger::ActionAfterSetVerbose() {
  if (generator) generator->SetVerboseLevel(verboseLevel);
}

void 
CDMSMultiGeneratorMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
  if (verboseLevel > 2) {
    G4cout << "CDMSMultiGeneratorMessenger::SetNewValue "
	   << cmd->GetCommandPath() << " " << value << G4endl;
  }

  CDMSMessengerBase::SetNewValue(cmd, value);	// To process global commands

  if (cmd == AddSpectrumCmd) AddSpectrum(value);
  if (cmd == AddIsotopeCmd)  AddIsotope(value);
  if (cmd == PrintCmd) generator->PrintSources(G4cout);
}


// Add hardcoded gamma spectrum using name, splitting off optional weight

void CDMSMultiGeneratorMessenger::AddSpectrum(const G4String& value) {
  if (verboseLevel>1) 
    G4cout << "CDMSMultiGeneratorMessenger::AddSpectrum " << value << G4endl;

  // Split string into two pieces: name and weight
  G4Tokenizer split(value);
  G4String name    = split();	// First whitespace delimited token
  G4String sWeight = split();	// Second whitespace delimited token
  G4double weight  = sWeight.empty() ? 1. : StoD(sWeight);

  if (name.contains("Ac228"))   generator->AddSpectrum(new Ac228Lines, weight);
  if (name.contains("Am241"))   generator->AddSpectrum(new Am241Lines, weight);
  if (name.contains("Bi212"))   generator->AddSpectrum(new Bi212Lines, weight);
  if (name.contains("Bi214"))   generator->AddSpectrum(new Bi214Lines, weight);
  if (name.contains("Co60"))    generator->AddSpectrum(new Co60Lines, weight);
  if (name.contains("K40"))     generator->AddSpectrum(new K40Lines, weight);
  if (name.contains("Pa234"))   generator->AddSpectrum(new Pa234Lines, weight);
  if (name.contains("Pb212"))   generator->AddSpectrum(new Pb212Lines, weight);
  if (name.contains("Pb214"))   generator->AddSpectrum(new Pb214Lines, weight);
  if (name.contains("Tl208"))   generator->AddSpectrum(new Tl208Lines, weight);
  if (name.contains("Uniform")) generator->AddSpectrum(new UniformLines, weight);
}


// Add radioactive nucleus to be decayed automatically

void CDMSMultiGeneratorMessenger::AddIsotope(const G4String& value) {
  if (verboseLevel>1) 
    G4cout << "CDMSMultiGeneratorMessenger::AddIsotope " << value << G4endl;

  // Split string into three pieces: Z, A and weight
  G4Tokenizer split(value);
  G4int Zval = StoI(split());	// First whitespace delimited token
  G4int Aval = StoI(split());	// Second whitespace delimited token
  G4String sWeight = split();	// Third whitespace delimited token
  G4double weight  = sWeight.empty() ? 1. : StoD(sWeight);

  if (Zval>0 && Aval>0) generator->AddNucleus(Zval, Aval, weight);
  else G4cerr << "Invalid isotope: Z=" << Zval << " A=" << Aval << G4endl;
}
