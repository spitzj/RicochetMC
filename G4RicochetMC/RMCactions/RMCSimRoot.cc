////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCSimRoot.cc                                       //
//  Description: class to fill Root trees from G4 hit collections.    //
//		 DMC units are eV and meters.                         //
//                                                                    //
//  Author:      Dennis Wright (SLAC)                                 //
//  Date:        15 July 2010                                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////
 
#include "RMCactions/RMCSimRoot.hh"
#include "G4Event.hh"

#include "TFile.h"
#include "TTree.h"

#include <sstream>
#include <iomanip>
#include <time.h>

 
RMCSimRoot::RMCSimRoot(G4int verbose)
  : verboseLevel(verbose), tfile(0), eventTree(0), vetoTree(0) {
  for (G4int i = 0; i < RMCSimRoot::NZIPS; i++) zipTree[i] = 0;
}

RMCSimRoot::~RMCSimRoot()
{
  WriteAndClose();	// Just in case the user forgot
}


void RMCSimRoot::ProcessEventInfo(const G4Event* evt)
{
  if (!evt) return;	// Sanity check

  evnum = evt->GetEventID();
  htpevt = 0.0;
  zppevt = 0.0;

  if (verboseLevel>1)
    G4cout << ">>> ProcessEventInfo for Event #" << evnum << G4endl;

  LoadEventData(evt);
  if(!eventTree) CreateEventTree();
  eventTree->Fill();
}


void RMCSimRoot::ProcessHits(const G4Event* evt, RMCZipHitsCollection* ZHC)
{
  if (!evt || !ZHC) return;	// Sanity check

  G4int nZipHits = ZHC->entries();

  if (verboseLevel>1)
    G4cout << ">>> ProcessHits with " << nZipHits << " ZIP hits" << G4endl;

  // Fill the zip hits
  for (G4int i = 0; i < nZipHits; i++) {
    LoadHitData(evt, (*ZHC)[i]);

    G4int idt = (*ZHC)[i]->GetReplicaNum() - 1;
    if (!zipTree[idt]) CreateAZipTree(idt);

    zipTree[idt]->Fill();
  }
}

void RMCSimRoot::ProcessHits(const G4Event* evt, RMCVetoHitsCollection* VHC)
{
  if (!evt || !VHC) return;	// Sanity check

  G4int nVetoHits = VHC->entries();

  if (verboseLevel>1)
    G4cout << ">>> ProcessHits with " << nVetoHits << " veto hits" << G4endl;

  for (G4int i = 0; i < nVetoHits; i++) {
    LoadHitData(evt, (*VHC)[i]);

    if (!vetoTree) CreateVetoTree();
    vetoTree->Fill();
  }
}

void RMCSimRoot::LoadHitData(const G4Event* evt, const RMCZipHit* ZipHit)
{
  if (!evt || !ZipHit) return;	// Sanity check

  ev = evt->GetEventID();
  dettype = 1.0;
  detnum = ZipHit->GetReplicaNum();
  ts = ZipHit->GetTrackID()*100000 + ZipHit->GetStepNum();
  p = ZipHit->GetParentID();
  type = ZipHit->GetPID();
  e1 = ZipHit->GetPreStepKE()/eV;
  d3 = ZipHit->GetEdep()/eV;
  
  px3 = ZipHit->GetPostStepMomentum().x()/eV;
  py3 = ZipHit->GetPostStepMomentum().y()/eV;
  pz3 = ZipHit->GetPostStepMomentum().z()/eV;
  x3 = ZipHit->GetPostStepPosition().x()/meter;
  y3 = ZipHit->GetPostStepPosition().y()/meter;
  z3 = ZipHit->GetPostStepPosition().z()/meter;
  t3 = ZipHit->GetPostStepTime()/ns;
  
  px1 = ZipHit->GetPreStepMomentum().x()/eV;
  py1 = ZipHit->GetPreStepMomentum().y()/eV;
  pz1 = ZipHit->GetPreStepMomentum().z()/eV;
  x1 = ZipHit->GetPreStepPosition().x()/meter;
  y1 = ZipHit->GetPreStepPosition().y()/meter;
  z1 = ZipHit->GetPreStepPosition().z()/meter;
  t1 = ZipHit->GetPreStepTime()/ns;
  empty = 0.0;
}

void RMCSimRoot::LoadHitData(const G4Event* evt, const RMCVetoHit* VetoHit) {
  if (!evt || !VetoHit) return;	// Sanity check

  vev = evt->GetEventID();
  vdettype = 1.0;
  vdetnum = VetoHit->GetReplicaNum();
  vts = VetoHit->GetTrackID()*100000 + VetoHit->GetStepNum();
  vp = VetoHit->GetParentID();
  vtype = VetoHit->GetPID();
  ve1 = VetoHit->GetPreStepKE()/eV;
  vd3 = VetoHit->GetEdep()/eV;
  
  vpx3 = VetoHit->GetPostStepMomentum().x()/eV;
  vpy3 = VetoHit->GetPostStepMomentum().y()/eV;
  vpz3 = VetoHit->GetPostStepMomentum().z()/eV;
  vx3 = VetoHit->GetPostStepPosition().x()/meter;
  vy3 = VetoHit->GetPostStepPosition().y()/meter;
  vz3 = VetoHit->GetPostStepPosition().z()/meter;
  vt3 = VetoHit->GetPostStepTime()/ns;
  
  vpx1 = VetoHit->GetPreStepMomentum().x()/eV;
  vpy1 = VetoHit->GetPreStepMomentum().y()/eV;
  vpz1 = VetoHit->GetPreStepMomentum().z()/eV;
  vx1 = VetoHit->GetPreStepPosition().x()/meter;
  vy1 = VetoHit->GetPreStepPosition().y()/meter;
  vz1 = VetoHit->GetPreStepPosition().z()/meter;
  vt1 = VetoHit->GetPreStepTime()/ns;
  vempty = 0.0;
}

void RMCSimRoot::LoadEventData(const G4Event* evt)
{
  if (!evt) return;

  G4PrimaryVertex* thisPrimaryVertex = evt->GetPrimaryVertex();
  G4PrimaryParticle* thisPrimaryParticle = thisPrimaryVertex->GetPrimary();

  primaryPID = thisPrimaryParticle->GetPDGcode();
  primaryX = thisPrimaryVertex->GetX0();
  primaryY = thisPrimaryVertex->GetY0();
  primaryZ = thisPrimaryVertex->GetZ0();
  primaryPx = thisPrimaryParticle->GetPx();
  primaryPy = thisPrimaryParticle->GetPy();
  primaryPz = thisPrimaryParticle->GetPz();
  primaryT0 = sqrt(pow(primaryPx,2) + pow(primaryPy,2) + pow(primaryPz,2) + pow(thisPrimaryParticle->GetMass(),2)) - thisPrimaryParticle->GetMass();
}


void RMCSimRoot::WriteAndClose()
{
  if (verboseLevel) G4cout << ">>> RMCSimRoot::WriteAndClose" << G4endl;

  for (G4int i = 0; i < RMCSimRoot::NZIPS; i++) {
    if (zipTree[i]) {
      zipTree[i]->Write();
      delete zipTree[i];
      zipTree[i] = 0;
    }
  }

  if (vetoTree) {
    vetoTree->Write();
    delete vetoTree;
    vetoTree = 0;
  }

  if (eventTree) {
    eventTree->Write();
    delete eventTree;
    eventTree = 0;
  }

  if (tfile) {
    tfile->Write("", TObject::kOverwrite);
    tfile->Close();
    delete tfile;
    tfile = 0;
  }
}


void RMCSimRoot::SetupFile(const G4String& prefix, G4int runNumber,
			    G4int labCode) {
  if (tfile) WriteAndClose();			// Replace previous data file

  G4String fileName = FileName(prefix, runNumber, labCode);

  if (verboseLevel)
    G4cout << " Creating RMC simulation ROOT file " << fileName << G4endl;

  tfile = new TFile(fileName,"RECREATE");	// Overwrite pervious version
  tfile->mkdir("G4SimDir");
  tfile->mkdir("G4SettingsInfoDir");
  tfile->cd("G4SimDir");
  CreateEventTree();
}


// Filename format is PREFIX_LLYYMMDD_RRRR.root (LL is 10 + lab code)

G4String RMCSimRoot::FileName(const G4String& prefix, G4int runNumber,
			       G4int labCode) const {
  time_t now;		  // Get <time.h> structure with local calendar data
  time(&now);
  struct tm* sNow = localtime(&now);

  G4int yymmdd = 100*(100*(sNow->tm_year%100)+(sNow->tm_mon+1))+sNow->tm_mday;

  std::ostringstream fn;

  // NOTE: The lab code must be one digit, and gets a prefixed "1"
  fn << prefix << "_" << std::setfill('0') << std::setw(2) << 10+(labCode%10)
     << std::setw(6) << yymmdd << "_" << std::setw(4) << runNumber << ".root";

  G4String fileName(fn.str());
  return fileName;
}


void RMCSimRoot::CreateEventTree()
{
  if (eventTree) return;	// Avoid leaks

  if (verboseLevel)
    G4cout << " eventName = mvevent"
	   << " eventDescription = G4 simulated event info" << G4endl;

  eventTree = new TTree("mvevent", "G4 simulated event info");
  eventTree->Branch("Event#", &evnum, "evnum/D");
  eventTree->Branch("HitsPerEvent", &htpevt, "htpevt/D");
  eventTree->Branch("ZipsPerEvent", &zppevt, "zppevt/D");
  eventTree->Branch("PrimaryPID", &primaryPID, "primaryPID/D");
  eventTree->Branch("PrimaryX", &primaryX, "primaryX/D");
  eventTree->Branch("PrimaryY", &primaryY, "primaryY/D");
  eventTree->Branch("PrimaryZ", &primaryZ, "primaryZ/D");
  eventTree->Branch("PrimaryPx", &primaryPx, "primaryPx/D");
  eventTree->Branch("PrimaryPy", &primaryPy, "primaryPy/D");
  eventTree->Branch("PrimaryPz", &primaryPz, "primaryPz/D");
  eventTree->Branch("PrimaryT0", &primaryT0, "primaryT0/D");
}


void RMCSimRoot::CreateVetoTree()
{ 
  if (vetoTree) return;		// Avoid leaks

  G4String vetoTreeName = "mcveto";
  G4String vetoTreeDescription = "G4 data for tree mcveto";

  if (verboseLevel)
    G4cout << " vetoName = " << vetoTreeName 
	   << " vetoDescription = " << vetoTreeDescription << G4endl; 

  vetoTree = new TTree(vetoTreeName, vetoTreeDescription);

  vetoTree->Branch("Empty", &vempty, "empty/D");
  vetoTree->Branch("EventNum", &vev, "ev/D");
  vetoTree->Branch("DetType", &vdettype, "dt/D");
  vetoTree->Branch("DetNum", &vdetnum, "dt/D");
  vetoTree->Branch("TrkStep", &vts, "ts/D");
  vetoTree->Branch("Parent", &vp, "p/D");
  vetoTree->Branch("PType", &vtype, "type/D");
  vetoTree->Branch("KE", &ve1, "e1/D");
  vetoTree->Branch("Edep", &vd3, "d3/D");
  vetoTree->Branch("Xmom3", &vpx3, "px3/D");
  vetoTree->Branch("Ymom3", &vpy3, "py3/D");
  vetoTree->Branch("Zmom3", &vpz3, "pz3/D");
  vetoTree->Branch("X3", &vx3, "x3/D");
  vetoTree->Branch("Y3", &vy3, "y3/D");
  vetoTree->Branch("Z3", &vz3, "z3/D");
  vetoTree->Branch("Time3", &vt3, "t3/D");
  vetoTree->Branch("Xmom1", &vpx1, "px1/D");
  vetoTree->Branch("Ymom1", &vpy1, "py1/D");
  vetoTree->Branch("Zmom1", &vpz1, "pz1/D");
  vetoTree->Branch("X1", &vx1, "x1/D");
  vetoTree->Branch("Y1", &vy1, "y1/D");
  vetoTree->Branch("Z1", &vz1, "z1/D");
  vetoTree->Branch("Time1", &vt1, "t1/D");
}


void RMCSimRoot::CreateAZipTree(G4int i)
{
  if (zipTree[i]) return;	// Avoid leaks

  std::ostringstream os;
  os << i+1;
  G4String zipTreeName = "mczip" + os.str();
  G4String zipTreeDescription = "G4 data for tree " + zipTreeName;

  if (verboseLevel)
    G4cout << " zipName = " << zipTreeName 
	   << " zipDescription = " << zipTreeDescription << G4endl; 

  zipTree[i] = new TTree(zipTreeName, zipTreeDescription);

  zipTree[i]->Branch("Empty", &empty, "empty/D");
  zipTree[i]->Branch("EventNum", &ev, "ev/D");
  zipTree[i]->Branch("DetType", &dettype, "dt/D");
  zipTree[i]->Branch("DetNum", &detnum, "dt/D");
  zipTree[i]->Branch("TrkStep", &ts, "ts/D");
  zipTree[i]->Branch("Parent", &p, "p/D");
  zipTree[i]->Branch("PType", &type, "type/D");
  zipTree[i]->Branch("KE", &e1, "e1/D");
  zipTree[i]->Branch("Edep", &d3, "d3/D");
  zipTree[i]->Branch("Xmom3", &px3, "px3/D");
  zipTree[i]->Branch("Ymom3", &py3, "py3/D");
  zipTree[i]->Branch("Zmom3", &pz3, "pz3/D");
  zipTree[i]->Branch("X3", &x3, "x3/D");
  zipTree[i]->Branch("Y3", &y3, "y3/D");
  zipTree[i]->Branch("Z3", &z3, "z3/D");
  zipTree[i]->Branch("Time3", &t3, "t3/D");
  zipTree[i]->Branch("Xmom1", &px1, "px1/D");
  zipTree[i]->Branch("Ymom1", &py1, "py1/D");
  zipTree[i]->Branch("Zmom1", &pz1, "pz1/D");
  zipTree[i]->Branch("X1", &x1, "x1/D");
  zipTree[i]->Branch("Y1", &y1, "y1/D");
  zipTree[i]->Branch("Z1", &z1, "z1/D");
  zipTree[i]->Branch("Time1", &t1, "t1/D");
}
