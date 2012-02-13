////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCSimData.cc                                       //
//  Description: class to output G4 hit collections to text file      //
//		 DMC units are eV and meters.                         //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        20 April 2011                                        //
//                                                                    //
////////////////////////////////////////////////////////////////////////
 
#include "RMCactions/RMCSimData.hh"
#include "G4Event.hh"
#include <fstream>
#include <sstream>
#include <iomanip>

 
RMCSimData::RMCSimData(G4int verbose)
  : verboseLevel(verbose), outputFile(0) {}

RMCSimData::~RMCSimData()
{
  WriteAndClose();	// Just in case the user forgot
}


void RMCSimData::SetupFile(const G4String& prefix, G4int runNumber,
			    G4int labCode) {
  if (outputFile) WriteAndClose();		// Replace previous data file

  G4String fileName = FileName(prefix, runNumber, labCode);

  if (verboseLevel)
    G4cout << " Creating RMC simulation data file " << fileName << G4endl;

  outputFile = new std::ofstream(fileName, std::ios_base::trunc);

  // Write header and column names to top of file
  *outputFile << "\n Run " << runNumber
	      << "\nEV" << "\tDT" << "\tTS" << "\tP" << "\tType" << "\tE1"
	      << "\tD3" << "\tPX3" << "\tPY3" << "\tPZ3" << "\tX3"
	      << "\tY3" << "\tZ3" << "\tT3" << "\tPX1" << "\tPY1"
	      << "\tPZ1" << "\tX1" << "\tY1" << "\tZ1" << "\tT1" << std::endl;
  
  outputFile->precision(10);
}


// Filename format is PREFIX_LLYYMMDD_RRRR.root (LL is 10 + lab code)

G4String RMCSimData::FileName(const G4String& prefix, G4int runNumber,
			       G4int labCode) const {
  time_t now;		  // Get <time.h> structure with local calendar data
  time(&now);
  struct tm* sNow = localtime(&now);

  G4int yymmdd = 100*(100*(sNow->tm_year%100)+(sNow->tm_mon+1))+sNow->tm_mday;

  std::ostringstream fn;

  // NOTE: The lab code must be one digit, and gets a prefixed "1"
  fn << prefix << "_" << std::setfill('0') << std::setw(2) << 10+(labCode%10)
     << std::setw(6) << yymmdd << "_" << std::setw(4) << runNumber << ".txt";

  G4String fileName(fn.str());
  return fileName;
}


void RMCSimData::ProcessHits(const G4Event* evt, RMCZipHitsCollection* ZHC)
{
  if (!evt || !ZHC) return;	// Sanity check

  G4int nZipHits = ZHC->entries();

  // Fill the zip hits
  for (G4int i = 0; i < nZipHits; i++) {
    LoadHitData(evt, (*ZHC)[i]);
    WriteHitData();
  }
}

void RMCSimData::ProcessHits(const G4Event* evt, RMCVetoHitsCollection* VHC)
{
  if (!evt || !VHC) return;	// Sanity check

  G4int nVetoHits = VHC->entries();

  for (G4int i = 0; i < nVetoHits; i++) {
    LoadHitData(evt, (*VHC)[i]);
    WriteHitData();
  }
}

void RMCSimData::LoadHitData(const G4Event* evt, const RMCZipHit* ZipHit)
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

void RMCSimData::LoadHitData(const G4Event* evt, const RMCVetoHit* VetoHit) {
  if (!evt || !VetoHit) return;	// Sanity check

  ev = evt->GetEventID();
  dettype = 1.0;
  detnum = VetoHit->GetReplicaNum();
  ts = VetoHit->GetTrackID()*100000 + VetoHit->GetStepNum();
  p = VetoHit->GetParentID();
  type = VetoHit->GetPID();
  e1 = VetoHit->GetPreStepKE()/eV;
  d3 = VetoHit->GetEdep()/eV;
  
  px3 = VetoHit->GetPostStepMomentum().x()/eV;
  py3 = VetoHit->GetPostStepMomentum().y()/eV;
  pz3 = VetoHit->GetPostStepMomentum().z()/eV;
  x3 = VetoHit->GetPostStepPosition().x()/meter;
  y3 = VetoHit->GetPostStepPosition().y()/meter;
  z3 = VetoHit->GetPostStepPosition().z()/meter;
  t3 = VetoHit->GetPostStepTime()/ns;
  
  px1 = VetoHit->GetPreStepMomentum().x()/eV;
  py1 = VetoHit->GetPreStepMomentum().y()/eV;
  pz1 = VetoHit->GetPreStepMomentum().z()/eV;
  x1 = VetoHit->GetPreStepPosition().x()/meter;
  y1 = VetoHit->GetPreStepPosition().y()/meter;
  z1 = VetoHit->GetPreStepPosition().z()/meter;
  t1 = VetoHit->GetPreStepTime()/ns;
  empty = 0.0;
}


void RMCSimData::WriteHitData() {
  if (!outputFile) return;

  *outputFile << ev << "\t" << detnum << "\t" << ts << "\t" << p << "\t"
	      << type << "\t" << e1 << "\t" << d3 << "\t"
	      << px3 << "\t" << py3 << "\t" << pz3 << "\t"
	      << x3 << "\t" << y3 << "\t" << z3 << "\t" << t3 << "\t"
	      << px1 << "\t" << py1 << "\t" << pz1 << "\t"
	      << x1 << "\t" << y1 << "\t" << z1 << "\t" << t1 << std::endl;
}


void RMCSimData::WriteAndClose() {
  // FIXME:  GCC 4.0.1 on MacOSX claims there's no "close()" function!
  //  if (outputFile) outputFile->close();
  delete outputFile;
  outputFile = 0;
}
