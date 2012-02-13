// $Id: RadiogenicSources.cc,v 1.2 2011/01/05 19:16:21 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RadiogenicSources.hh                                 //
//  Description: Consolidated class for all radioactive decays        //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        12 December 2010                                     //
//                                                                    //
//  20101212  M. Kelsey -- Copy of CDMSmain/CDMS_ParticleSource.      //
////////////////////////////////////////////////////////////////////////

#include "CDMSsources/RadiogenicSources.hh"
#include "CDMSsources/RadiogenicSourcesMessenger.hh"

#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PrimaryParticle.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "Randomize.hh"
#include <cmath>
#include <fstream>
using std::ifstream;


RadiogenicSources::RadiogenicSources()
  : CDMSVSourceConstruction(), SourcePosType("Point"), Shape("NULL"),
    centre_z(0.), centre_r(0.), centre_phi(0.), halfz(25.*cm),
    Radius(0.), Confine(false), VolName("NULL"), AngDistType("iso"),
    MinTheta(0.), MaxTheta(pi*rad), MinPhi(0.), MaxPhi(twopi*rad),
    Theta(0.), Phi(0.), EnergyDisType("Mono"), ceFlag(1), MonoEnergy(1*MeV),
    NumberOfParticlesToBeGenerated(1), particle_definition(0),
    particle_momentum_direction(1.,0.,0.),
    particle_energy(1.*MeV), particle_charge(0.), particle_time(0.),
    gNavigator(0), bot(0), zyl(0), top(0), ass(0), count(0),
    Bi210rank(0), C14rank(0), Cf252rank(0), Cs137rank(0), Ga68rank(0),
    K40rank1(0), K40rank2(0), Pb210rank1(0), Pb210rank2(0), Sb125rank(0),
    Tl204rank(0), particleTable(0),
    theMessenger(new RadiogenicSourcesMessenger(this)) {
  particleTable = G4ParticleTable::GetParticleTable();
  SetParticleDefinition(particleTable->FindParticle("geantino"));

  gNavigator =
    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
}

RadiogenicSources::~RadiogenicSources() {
  delete theMessenger;
}


void RadiogenicSources::SetPosDisType(const G4String& PosType) 
{
  SourcePosType = PosType;

  if (SourcePosType == "PointRandomCopperTower")
    myMaterial = "CopperTower";
  else if (SourcePosType == "PointRandomInnerShieldAir")
    myMaterial = "InnerShieldAir";
  else  if (SourcePosType == "PointRandomCopperCans")
    myMaterial = "CopperCans";
  else  if (SourcePosType == "PointRandomZIPs")
    myMaterial = "ZIPs";
  else myMaterial = "vacuum";

  if (verboseLevel) G4cout << "MyMaterial: " << myMaterial << G4endl;
}

void RadiogenicSources::ConfineSourceToVolume(const G4String& Vname)
{
  VolName = Vname;
  if(verboseLevel == 2) G4cout << VolName << G4endl;

  // checks if selected volume exists
  G4VPhysicalVolume *tempPV      = NULL;
  G4PhysicalVolumeStore *PVStore = 0;
  G4String theRequiredVolumeName = VolName;
  PVStore = G4PhysicalVolumeStore::GetInstance();
  G4int i = 0;
  G4bool found = false;
  if(verboseLevel == 2) G4cout << PVStore->size() << G4endl;

  // recasting required since PVStore->size() is actually a signed int...
  while (!found && i<(G4int)PVStore->size())
    {
      tempPV = (*PVStore)[i];
      found  = tempPV->GetName() == theRequiredVolumeName;
      if(verboseLevel == 2)
	G4cout << i << " " << " " << tempPV->GetName() 
	       << " " << theRequiredVolumeName << " " << found << G4endl;
      if (!found)
	{i++;}
    }

  // found = true then the volume exists else it doesnt.
  if(found == true) 
  {
    if(verboseLevel >= 1)
      G4cout << "Volume " << VolName << " exists" << G4endl;
    Confine = true;
  }
  else if(VolName=="NULL")
    Confine = false;
  else 
  {
    G4cout << " **** Error: Volume does not exist **** " << G4endl;
    G4cout << " Ignoring confine condition" << G4endl;
    VolName = "NULL";
    Confine = false;
  }
}

void RadiogenicSources::SetAngDistType(const G4String& atype)
{
  AngDistType = atype;
}


void RadiogenicSources::GeneratePointSource()
{
  // Generates Points given the point source.
  if(SourcePosType == "Point")
    {
      //    particle_position = CentreCoords;
      particle_position.setX(centre_r * cos(centre_phi));
      particle_position.setY(centre_r * sin(centre_phi));
      particle_position.setZ(centre_z);

      particle_position += CentreCoords;
      //      G4cout << "Position: " << particle_position << G4endl;
    }

  else
    if(verboseLevel >= 1)
      G4cout << "Error SourcePosType is not set to Point" << G4endl;
}


void RadiogenicSources::GeneratePointsInVolume()
{
  G4ThreeVector RandPos;
  G4double x=0., y=0., z=0.;
  
  if(SourcePosType != "Volume" && verboseLevel >= 1)
    G4cout << "Error SourcePosType not Volume" << G4endl;
  
  if(Shape == "Sphere")
  {
    x = Radius*2.;
    y = Radius*2.;
    z = Radius*2.;
    while(((x*x)+(y*y)+(z*z)) > (Radius*Radius)) 
    {
      x = G4UniformRand();
      y = G4UniformRand();
      z = G4UniformRand();
      
      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
      z = (z*2.*Radius) - Radius;
    }
  }

  else if(Shape == "Cylinder") 
  {
    x = Radius*2.;
    y = Radius*2.;
    while(((x*x)+(y*y)) > (Radius*Radius)) 
    {
      x = G4UniformRand();
      y = G4UniformRand();
      z = G4UniformRand();
      x = (x*2.*Radius) - Radius;
      y = (y*2.*Radius) - Radius;
      z = (z*2.*halfz) - halfz;
    }
  }

  else
    G4cout << "Error: Volume Shape Does Not Exist" << G4endl;

  RandPos.setX(x);
  RandPos.setY(y);
  RandPos.setZ(z);
  particle_position = CentreCoords + RandPos;
}


G4bool RadiogenicSources::IsSourceConfined()
{
  // Method to check point is within the volume specified
  if(Confine == false)
    G4cout << "Error: Confine is false" << G4endl;
  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector *ptr;
  ptr = &null;

  // Check particle_position is within VolName
  G4VPhysicalVolume *theVolume;
  theVolume=gNavigator->LocateGlobalPointAndSetup(particle_position,ptr,true);
  G4String theVolName = theVolume->GetName();
  if(theVolName == VolName) 
  {
    if(verboseLevel >= 1)
      G4cout << "Particle is in volume " << VolName << G4endl;
    return(true);
  }
  else
    return(false);
}

void RadiogenicSources::SetRandomPosInMaterial()
{
/*  G4String foundMaterial,VolName;
  G4bool happy = false;
  G4ThreeVector randomPosition;
  G4double rand1,rand2,rand3,rand4;
  G4double radius,angle,x,y,z;
  G4double Rmin=0,DeltaR=0,Zmin=0,DeltaZ=0;

  // Below: for Xe experiment (create a particule in different volumes of the experiment 
  // TO DO : change this to create same kind of things for cdms...
  
  // for outer cryostat outer shell only

  G4VPhysicalVolume *theVolume;

    rand1 = G4UniformRand();

    if (rand1>=0. && rand1<=0.178)
    {// Parameter settings for bottom flange
	bot++;
	//	G4cout << "BOT: Random Number: " << rand1 << G4endl;
	Rmin   =  0.0*cm;
	DeltaR = 21.0*cm;
	Zmin   = -1.45*cm;
	DeltaZ = 2.9*cm;}
    else if(rand1>0.178 && rand1<=0.754)
    {// Parameter settings for cryostat cylinder
      zyl++;
      //      G4cout << "ZYL: Random Number: " << rand1 << G4endl;
      Rmin   = 16.8*cm;
      DeltaR =  8.7*cm;
      Zmin   =  1.5*cm;
      DeltaZ = 62.2*cm;}
    else if(rand1>0.754 && rand1<=0.921)
    { // Parameter settings for top flange      
      top++;
      //      G4cout << "TOP: Random Number: " << rand1 << G4endl;
      Rmin   =  0.0*cm;
      DeltaR = 21.0*cm;
      Zmin   = 63.7*cm;
      DeltaZ = 2.9*cm;}
    else if(rand1>0.921 && rand1<=1.0)
    {// Parameter settings for top flange assembly	
	ass++;
	//	G4cout << "ASS: Random Number: " << rand1 << G4endl;
	Rmin   =  0.0*cm;
	DeltaR = 16.2*cm;
	Zmin   = 66.7*cm;
	DeltaZ = 20.6*cm;}
    } 
    //    G4cout << "Random Number: " << rand1 << G4endl;
    //    G4cout << "bot,zyl,top,ass: " << bot << " " << zyl << " " << top << " " << ass << G4endl; 

  G4int ihap = 0;
  x=0.;y=0.;z=0.;

  while (!happy)
  {      
      ihap++;
      
      rand1 = G4UniformRand();
      rand2 = G4UniformRand();
      rand3 = G4UniformRand();
      rand4 = G4UniformRand();

      radius = Rmin + rand1*DeltaR;
      angle  = rand2*MaxPhi;
  
      x = radius*cos(angle);
      y = radius*sin(angle);

      z = Zmin + rand4*DeltaZ;
      	
      randomPosition.setX(x);
      randomPosition.setY(y);
      randomPosition.setZ(z);

  // Check particle_position is within Material
 
      theVolume=gNavigator->LocateGlobalPointAndSetup(randomPosition,0,true,true);
      VolName = theVolume->GetName();
      foundMaterial = theVolume->GetLogicalVolume()->GetMaterial()->GetName();							
	//     G4cout  << "Volume: " << VolName << " Material: " << foundMaterial 
	// << " at Position: " << randomPosition << G4endl;				
      happy = (foundMaterial == myMaterial);
      //      G4cout << "Happy? " << happy << G4endl;
  }

  //  G4cout << "********************* Finished Looping ********************" << G4endl;
  //    G4cout << "Number of loops: " << ihap << G4endl;

  particle_position.setX(x);
  particle_position.setY(y);
  particle_position.setZ(z);

  //  G4cout << "Volume: " << VolName << " Material: " << foundMaterial << " Particle Position: " << particle_position << G4endl; 
*/
}


void RadiogenicSources::GenerateIsotropicFlux()
{
  G4double rndm, rndm2;
  G4double px, py, pz;

  G4double sintheta, sinphi, costheta, cosphi;
  rndm = G4UniformRand();
  costheta = cos(MinTheta) - rndm * (cos(MinTheta) - cos(MaxTheta));
  sintheta = sqrt(1. - costheta*costheta);
  
  rndm2 = G4UniformRand();
  Phi = MinPhi + (MaxPhi - MinPhi) * rndm2; 
  sinphi = sin(Phi);
  cosphi = cos(Phi);

  px = -sintheta * cosphi;
  py = -sintheta * sinphi;
  pz = -costheta;

  G4double ResMag = sqrt((px*px) + (py*py) + (pz*pz));
  px = px/ResMag;
  py = py/ResMag;
  pz = pz/ResMag;

  particle_momentum_direction.setX(px);
  particle_momentum_direction.setY(py);
  particle_momentum_direction.setZ(pz);

  // particle_momentum_direction now holds unit momentum vector.
  if(verboseLevel >= 2)
    G4cout << "Generating isotropic vector: " << particle_momentum_direction << G4endl;
}


void RadiogenicSources::SetCEflag(G4int generateCEFlag)
{
  ceFlag = generateCEFlag;
  
  if (ceFlag == 1) {
      G4cout << "In this simulation conversion and beta electrons will be started" << G4endl;
      G4cout << "In the following decay chains we create only the photons, not the conversion electrons. " << G4endl;
      G4cout << "In an analysis one should use the following multiplication factor to estimate the actual" << G4endl;
      G4cout << "number of events:" << G4endl;
      G4cout << " U234 multiply No. of Events by                       3.6" << G4endl;
      G4cout << "Th230 multiply No. of Events by                       4.2" << G4endl;
      G4cout << "Th228 multiply No. of Events by                       3.7" << G4endl;
      G4cout << "Ra226 multiply No. of Events by                      18.0" << G4endl;
      G4cout << "Ra224 multiply No. of Events by                      20.0" << G4endl;
    }
  else {
      G4cout << "This simulation starts ONLY GAMMAs" << G4endl; 
      G4cout << "In the following decay chains we create only the photons, not the conversion electrons. " << G4endl;
      G4cout << "In an analysis one should use the following multiplication factor to estimate the actual" << G4endl;
      G4cout << "number of events:" << G4endl;
      G4cout << " U234 multiply No. of Events by                      1212" << G4endl;
      G4cout << "Th230 multiply No. of Events by                      1055" << G4endl;
      G4cout << "Th228 multiply No. of Events by                      84.6" << G4endl;
      G4cout << "Ra226 multiply No. of Events by                      28.5" << G4endl;
      G4cout << "Ra224 multiply No. of Events by                      25.6" << G4endl;
      G4cout << "Pb210 multiply No. of Events by                      24.7" << G4endl;
      G4cout << " K40  multiply No. of Events by                      84.6" << G4endl;
    }

}


void RadiogenicSources::SetEnergyDisType(const G4String& DisType)
{ 
  G4int nf=0;
  EnergyDisType = DisType;
  for(G4int ii=0;ii<2000;ii++) {
	beta_eng[ii]=0.;beta_prob[ii]=1.;}
  count = 0; nlfile[0]=0;
  
   // If a beta spectrum file needs to be read:
     // call ReadEnergyProbaFile for a 2 columns (Energy vs Proba) file, to calculate integral of proba 
     // call ReadIntProbaFile for 1 column (integrated Proba) files with energy bin deduced from 1st line

  if ((EnergyDisType == "Bi210") | (EnergyDisType == "PbBiPo")) {
	ReadEnergyProbaFile("beta_spectra/BI210_PO210.DAT",nf+1);
	Bi210rank = 2*nf; nf+=1;}
	
  if (EnergyDisType == "C14") {
	ReadEnergyProbaFile("beta_spectra/14C_BETA.DAT",nf+1);
	C14rank = 2*nf; nf+=1;}

  if (EnergyDisType == "Cf252") {
        ReadEnergyProbaFile("neutron_spectra/cf252_neutron.dat",nf+1);
        Cf252rank =2*nf; nf+=1;}

  if (EnergyDisType == "Cs137") {
	ReadEnergyProbaFile("beta_spectra/137CS_BETA.DAT",nf+1);
	Cs137rank = 2*nf; nf+=1;}
	
  if (EnergyDisType == "Ga68") {
	ReadEnergyProbaFile("beta_spectra/GA68.DAT",nf+1);
	Ga68rank = 2*nf; nf+=1;}
	
  if (EnergyDisType == "K40") {
	ReadEnergyProbaFile("beta_spectra/40K.DAT",nf+1);
	ReadEnergyProbaFile("beta_spectra/40K_BETAP.DAT",nf+2);
	K40rank1 = 2*nf; K40rank2 = 2*(nf+1); nf+=2;}
	
  if ((EnergyDisType == "Pb210") | (EnergyDisType == "PbBiPo")) {
	// ReadEnergyProbaFile("beta_spectraPB210_BI210_BETAS.DAT",nf+1);
	// Pb210rank1 = nf; Pb210rank2 = nf; nf+=1;}
  	ReadIntProbaFile("beta_spectra/Pb210_17keV_branch.DAT",nf+1);
	ReadIntProbaFile("beta_spectra/Pb210_64keV_branch.DAT",nf+2);
	Pb210rank1 = 2*nf; Pb210rank2 = 2*(nf+1); nf+=2;}
  		
  if (EnergyDisType == "Sb125") {
	ReadEnergyProbaFile("beta_spectra/125SB_BETA.DAT",nf+1);
	Bi210rank = 2*nf; nf+=1;}
	
  if (EnergyDisType == "Tl204") {
	ReadEnergyProbaFile("beta_spectra/204Tl_BETA.DAT",nf+1);
	Bi210rank = 2*nf; nf+=1;}
	
  // For the other nuclides, no beta spectrum ==> go directly generate their gammas
}


void RadiogenicSources::ReadEnergyProbaFile(const G4String& EPFile, G4int nf)
{
  G4int iBeta;
  G4double e=0., p=0.;
  G4double dx,dy;
  G4double energy[2000], prob[2000];
  G4int cnt = 0;
  count=nlfile[2*nf-2];
 
  G4cout << "Read 2 columns Energy [MeV] vs Probability in file" << EPFile << G4endl;
  ifstream input_file (EPFile);
  if (! input_file.is_open())
     {G4cout << "Error opening " << EPFile << G4endl;}
  while (! input_file.eof()) {
       input_file >> e >> p;
       energy[cnt] = e; prob[cnt] = p;
       //G4cout << "index:" << cnt << " E:" << energy[count] MeV << " proba:" << prob[cnt] << G4endl;
       cnt++;
  }
  input_file.close();
  G4cout << cnt << " lines of data" << G4endl;

  beta_eng[count] = 0.; beta_prob[count] = 0.;
  for(iBeta=0;iBeta<cnt-1;iBeta++) {
     dx = energy[iBeta+1] - energy[iBeta];
     dy = prob[iBeta+1] - prob[iBeta];
     beta_eng[count+iBeta+1] = energy[iBeta] * 1000; // convert energy from MeV to keV (check data files !)
     beta_prob[count+iBeta+1] = beta_prob[count+iBeta] + dx*prob[iBeta] + 0.5*dx*dy;
  }
  for(iBeta=0;iBeta<cnt;iBeta++) {
    beta_prob[count+iBeta] /= beta_prob[count+cnt-1];
    G4cout << "index:" << count+iBeta << " beta_eng: " << beta_eng[count+iBeta] << "keV ; beta_prob: " << beta_prob[count+iBeta] << G4endl;
  }
  count+=cnt; nlfile[2*nf-1] = count-2; nlfile[2*nf] = count-1;
}


void RadiogenicSources::ReadIntProbaFile(const G4String& EPFile, G4int nf)
{
  G4double eb=0., lb=0., p=0.;

  G4cout << "Read Energy bin [keV] plus left bin [keV] plus Probabilities in 1 column file" << EPFile << G4endl;
  count=nlfile[2*nf-2];

  ifstream input_file (EPFile);
  if (! input_file.is_open())
     {G4cout << "Error opening " << EPFile << G4endl;}
  while (! input_file.eof()) {
       count++;
       if ((count-nlfile[2*nf-2])==1) {
	   input_file >> eb;}     // 1st line = bin width in keV
       else if ((count-nlfile[2*nf-2])==2) {
	   input_file >> lb;       // 2nd line = lefthand edge of bin 1 in keV
	   beta_eng[count-2] = lb; beta_prob[count-2] = 0.;
	   G4cout << "index:" << count-2 << " E:" << beta_eng[count-2] << " proba:" << beta_prob[count-2] << G4endl;}
       else {
	   input_file >> p;          // other lines = list of integrated probabilities
	   beta_eng[count-2] = ((count-2-nlfile[2*nf-2]) * eb + lb );
	   beta_prob[count-2] = p;
	   G4cout << "index:" << count-2 << " E:" << beta_eng[count-2] << " proba:" << beta_prob[count-2] << G4endl;
       }
  }
  input_file.close();
  G4cout << count-2-nlfile[2*nf-2] << " lines of data" << G4endl;
  nlfile[2*nf-1] = count-2; nlfile[2*nf] = count-1;
}


double RadiogenicSources::BetaEnergyFromSpectrum(G4int cstart, G4int cend)
{
   // Emit the initial beta from the correct spectrum
   G4double r_eb, r_be;
   G4double part_eng = 0.;
   // random numbers to select energy bin (r_eb) and an interpolated energy in that bin to avoid discretisation (r_be)
   r_eb = G4UniformRand(); r_be = G4UniformRand();
   // G4cout << "cstart=" << cstart << " cend=" << cend << " r_eb=" << r_eb << G4endl;
   for (int k=cstart;k<=cend;k++)  {
		if ( (r_eb >= beta_prob[k]) and (r_eb < beta_prob[k+1]) ) {
		// G4cout << "k=" << k << " beta_eng[k]=" << beta_eng[k] << G4endl;
		part_eng = beta_eng[k] + r_be*(beta_eng[k+1]-beta_eng[k]); k=cend+1;}
   }
   // G4cout << "energy of primary beta: " << part_eng << G4endl;
  return part_eng; 
}


void RadiogenicSources::SetParticleDefinition
  (G4ParticleDefinition* aParticleDefinition)
{
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
}

void RadiogenicSources::CreatePrimaryParticle(G4PrimaryVertex* vt)
{
  G4double mass =  particle_definition->GetPDGMass();
  G4double energy = particle_energy + mass;
  G4double pmom = sqrt(energy*energy-mass*mass);
  G4double px = pmom*particle_momentum_direction.x();
  G4double py = pmom*particle_momentum_direction.y();
  G4double pz = pmom*particle_momentum_direction.z();

  G4PrimaryParticle* particle =
    new G4PrimaryParticle(particle_definition,px,py,pz);
  particle->SetMass( mass );
  particle->SetCharge( particle_charge );
  particle->SetPolarization(particle_polarization.x(),
			    particle_polarization.y(),
			    particle_polarization.z());
  vt->SetPrimary( particle );
}


void RadiogenicSources::GenerateAc228Decay(G4PrimaryVertex* vert)
{
 G4double rndm;

 SetParticleDefinition(particleTable->FindParticle("gamma"));
 
 top:  rndm = G4UniformRand();
  if (rndm < 0.0030)
    goto two010;
  else if (rndm < 0.0054)
    goto one945;
  else if (rndm < 0.0070)
    goto  one760;
  else if (rndm < 0.0104)
    goto  one744;
  else if (rndm < 0.0116)
    goto  one735;
  else if (rndm < 0.0273)
    goto  one724;
  else if (rndm < 0.0533)
    goto  one688;
  else if (rndm < 0.0303)
    goto  one684;
  else if (rndm < 0.0549)
    goto  one683;
  else if (rndm < 0.0667)
    goto  one646;
  else if (rndm < 0.1085)
    goto  one643;
  else if (rndm < 0.1163)
    goto  one638;
  else if (rndm < 0.1278)
    goto  one539;
  else if (rndm < 0.1296)
    goto  one531;
  else if (rndm < 0.2106)
    goto  one450;
  else if (rndm < 0.2150)
    goto  one432;
  else if (rndm < 0.2240)
    goto  one344;
  else if (rndm < 0.2261)
    goto  one227;
  else if (rndm < 0.2241)
    goto  one175;
  else if (rndm < 0.2364)
    goto one168;
  else if (rndm < 0.2718)
    goto  one153;
  else if (rndm < 0.3278)
    goto  one123;
  else if (rndm < 0.3860)
    goto one091;
  else if (rndm < 0.3877)
    goto  one023;
  else if (rndm < 0.4177)
    goto  one016;
  else if (rndm < 0.7307)
    goto nine69;
  else if (rndm < 0.7329)
    goto eight74;
  else if (rndm < 0.8489)
    goto three96;
  else if (rndm < 0.8543)
    goto three28;
  else if (rndm < 0.8733)
    goto one87;
  else if (rndm < 0.9733)
    goto  fifty8;
  else 
    goto top;

 two010: rndm = G4UniformRand();
  if (rndm < 0.2062) {
	particle_energy = 1955.52*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.3582) {
	particle_energy = 1823.22*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.5137) {
	particle_energy = 1040.92*keV; CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.6860) {
	particle_energy =  987.41*keV; CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.6978) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 918.97*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one091;}
  else if (rndm < 0.7700) {
	particle_energy =  918.97*keV; CreatePrimaryParticle(vert);
	goto one091;}
  else if (rndm < 0.8633) {
	particle_energy =  887.33*keV; CreatePrimaryParticle(vert);
	goto one123;}  
  else if (rndm < 0.8863) {
	particle_energy =  372.57*keV; CreatePrimaryParticle(vert);
	goto one638;}
  else if (rndm < 0.9863) {
	particle_energy =  214.85*keV; CreatePrimaryParticle(vert);
	goto one795;}

 one945: rndm = G4UniformRand();
  if (rndm < 0.3928) {
	particle_energy = 1887.10*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.5432) {
	particle_energy = 1758.12*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.5656) {
	particle_energy = 1000.69*keV; CreatePrimaryParticle(vert);
	goto nine44;}
  else if (rndm < 0.7811) {
	particle_energy =  975.95*keV; CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.7856) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 975.95*keV;  CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.8373) {
	particle_energy =  853.71*keV; CreatePrimaryParticle(vert);
	goto one091;}
  else if (rndm < 0.9192) {
	particle_energy =  776.56*keV; CreatePrimaryParticle(vert);
	goto one168;}  
  else {
	particle_energy =  718.48*keV; CreatePrimaryParticle(vert);
	goto one227;}

 one795: rndm = G4UniformRand();
  if (rndm < 0.4474) {
	particle_energy = 1738.22*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.6447) {
	particle_energy = 1609.47*keV; CreatePrimaryParticle(vert);
	goto one87;}  
  else {
	particle_energy = 1276.69*keV; CreatePrimaryParticle(vert);
	goto five19;}

 one760: rndm = G4UniformRand();
  if (rndm < 0.3062) {
	particle_energy = 1702.40*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.5154) {
	particle_energy = 1573.20*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.6619) {
	particle_energy =  791.49*keV; CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.8973) {
	particle_energy =  737.72*keV; CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.9160) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 737.72*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}
  else {
	particle_energy =  416.30*keV; CreatePrimaryParticle(vert);
	goto one344;}

 one744: rndm = G4UniformRand();
  if (rndm < 0.2875) {
	particle_energy = 1686.21*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.8263) {
	particle_energy = 1557.11*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.8648) {
	particle_energy = 1365.70*keV; CreatePrimaryParticle(vert);
	goto three78;}
  else if (rndm < 0.9114) {
	particle_energy = 1347.50*keV; CreatePrimaryParticle(vert);
	goto three96;}
  else {
	particle_energy =  399.62*keV; CreatePrimaryParticle(vert);
	goto one344;}

 one735: rndm = G4UniformRand();
  if (rndm < 0.4054) {
	particle_energy = 1677.67*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.6892) {
	particle_energy = 1548.64*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.8417) {
	particle_energy = 1357.78*keV; CreatePrimaryParticle(vert);
	goto three78;}
  else {
	particle_energy = 1217.03*keV; CreatePrimaryParticle(vert);
	goto five19;}

 one724: rndm = G4UniformRand();
  if (rndm < 0.0189) {
	particle_energy = 1724.21*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.1357) {
	particle_energy = 1666.52*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.1663) {
	particle_energy = 1537.89*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.8058) {
	particle_energy =  755.32*keV; CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.8536) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 755.32*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.9669) {
	particle_energy =  701.75*keV; CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.9772) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 701.75*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}
  else if (rndm < 0.9922) {
	particle_energy =  548.73*keV; CreatePrimaryParticle(vert);
	goto one175;}
  else if (rndm < 0.9937) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 548.73*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one175;}
  else if (rndm < 0.9976) {
	particle_energy =  497.49*keV; CreatePrimaryParticle(vert);
	goto one227;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 497.49*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one227;}

 one688: rndm = G4UniformRand();
  if (rndm < 0.6220) {
	particle_energy = 1630.63*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.8054) {
	particle_energy = 1501.58*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.8082) {
	particle_energy =  813.77*keV; CreatePrimaryParticle(vert);
	goto eight74;}
  else if (rndm < 0.8186) {
	particle_energy =  672.00*keV; CreatePrimaryParticle(vert);
	goto one016;}
  else if (rndm < 0.8222) {
	particle_energy =   42.46*keV; CreatePrimaryParticle(vert);
	goto one643;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 42.46*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one643;}

 one684: rndm = G4UniformRand();
  if (rndm < 0.4924) {
	particle_energy = 1287.68*keV; CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.8895) {
	particle_energy = 1164.50*keV; CreatePrimaryParticle(vert);
	goto five19;}
  else if (rndm < 0.9848) {
	particle_energy =  457.17*keV; CreatePrimaryParticle(vert);
	goto one227;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 457.17*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one227;}

 one683: rndm = G4UniformRand();
  if (rndm < 0.2230) {
	particle_energy = 1624.99*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.9797) {
	particle_energy = 1495.91*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else {
	particle_energy = 1022.53*keV; CreatePrimaryParticle(vert);
	goto one023;}

 one646: rndm = G4UniformRand();
  if (rndm < 0.7037) {
	particle_energy = 1588.21*keV; CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.8754) {
	particle_energy = 1459.14*keV; CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.8891) {
	particle_energy = 1250.00*keV; CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.9028) {
	particle_energy =  677.11*keV; CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.9031) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 677.11*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.7700) {
	particle_energy =  918.97*keV; CreatePrimaryParticle(vert);
	goto one091;}
  else if (rndm < 0.9169) {
	particle_energy =  666.45*keV; CreatePrimaryParticle(vert);
	goto nine80;}  
  else if (rndm < 0.9177) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 666.45*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine80;}
  else if (rndm < 0.9276) {
	particle_energy =  629.45*keV; CreatePrimaryParticle(vert);
	goto one016;}
  else if (rndm < 0.9283) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 629.45*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one016;}
  else if (rndm < 0.9314) {
	particle_energy =  623.27*keV; CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.9416) {
	particle_energy =  555.12*keV; CreatePrimaryParticle(vert);
	goto one091;}
 else if (rndm < 0.9426) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 555.12*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one091;}
  else if (rndm < 0.9655) {
	particle_energy =  523.13*keV; CreatePrimaryParticle(vert);
	goto one123;}
  else if (rndm < 0.9658) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 523.13*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one123;}
  else if (rndm < 0.9710) {
	particle_energy =  492.37*keV; CreatePrimaryParticle(vert);
	goto one153;}
  else if (rndm < 0.9717) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 492.37*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one153;}
  else if (rndm < 0.9763) {
	particle_energy =  419.42*keV; CreatePrimaryParticle(vert);
	goto one227;}
  else if (rndm < 0.9785) {
	particle_energy =  114.56*keV; CreatePrimaryParticle(vert);
	goto one531;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 114.56*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one531;}

 one643: rndm = G4UniformRand();
  if (rndm < 0.0195) {
	particle_energy = 1313.34*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.6593) {
	particle_energy = 1247.08*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.7078) {
	particle_energy =  699.08*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine44;}
  else if (rndm < 0.8413) {
	particle_energy =  674.63*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine68;}
  else if (rndm < 0.9487) {
	particle_energy =  520.15*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one123;}
  else if (rndm < 0.9666) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 520.15*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one123;}
  else if (rndm < 0.9956) {
	particle_energy =  474.75*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one168;}  
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 474.75*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one168;}

 one638: rndm = G4UniformRand();
  if (rndm < 0.4058) {
	particle_energy = 1638.28*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9257) {
	particle_energy = 1580.54*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.9350) {
	particle_energy = 1451.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.9457) {
	particle_energy = 1309.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.9886) {
	particle_energy =  510.65*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one123;}
  else {
	particle_energy =  470.25*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one168;}

 one539: rndm = G4UniformRand();
  if (rndm < 0.0482) {
	particle_energy = 1142.85*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.8067) {
	particle_energy =  570.91*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine68;}
  else if (rndm < 0.9251) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 570.91*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one87;}
  else if (rndm < 0.9895) {
	particle_energy =  416.30*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one123;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 416.30*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one168;}

 one531: rndm = G4UniformRand();
  if (rndm < 0.0011) {
	particle_energy = 1344.59*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.0024) {
	particle_energy = 1135.24*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.1089) {
	particle_energy =  562.50*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.1157) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 562.50*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.1731) {
	particle_energy =  508.96*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.1799) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 508.96*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}
  else if (rndm < 0.1952) {
	particle_energy =  440.45*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one091;}  
  else if (rndm < 0.2000) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 440.45*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one091;}
  else if (rndm < 0.2031) {
	particle_energy =  377.99*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one153;}
  else if (rndm < 0.2040) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 377.99*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one153;}
  else if (rndm < 0.3604) {
	particle_energy =   99.49*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one432;}  
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 99.49*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one432;}

 one450: rndm = G4UniformRand();
  if (rndm < 0.0046) {
	particle_energy = 1054.11*keV; CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.3229) {
	particle_energy =  327.99*keV; CreatePrimaryParticle(vert);
	goto one123;}
  else if (rndm < 0.3594) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 327.99*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one123;}
  else if (rndm < 0.5176) {
	particle_energy =  281.92*keV; CreatePrimaryParticle(vert);
	goto one168;}
  else if (rndm < 0.6346) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 281.92*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one168;}
  else if (rndm < 0.7730) {
	particle_energy =  223.85*keV; CreatePrimaryParticle(vert);
	goto one227;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 223.95*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one227;}  

 one432: rndm = G4UniformRand();
  if (rndm < 0.0019) {
	particle_energy = 1374.19*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.0153) {
	particle_energy = 1245.05*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.1740) {
	particle_energy = 1103.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.6215) {
	particle_energy =  463.01*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.6530) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 463.01*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.6552) {
	particle_energy =  452.47*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine80;}
  else if (rndm < 0.6556) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 452.47*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine80;}  
  else if (rndm < 0.9189) {
	particle_energy =  409.46*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.9413) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 409.46*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}
  else if (rndm < 0.9926) {
	particle_energy =  340.97*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one091; }
  else if (rndm < 0.9941) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 340.97*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one432;}
  else if (rndm < 0.9997) {
	particle_energy =  263.58*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one168;}  
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 263.58*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
      goto one168;}

 one344: rndm = G4UniformRand();
  if (rndm < 0.2017) {
		particle_energy = 1268.27*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.2299) {
	particle_energy = 1157.14*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.3051) {
	particle_energy = 1016.44*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.7336) {
	particle_energy =  957.98*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.9374) {
	particle_energy =  824.94*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto five19;}
  else if (rndm < 0.9404) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 824.94*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto five19;}  
  else if (rndm < 0.9927) {
	particle_energy =  168.65*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one175;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 168.65*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one175;}

 one227: rndm = G4UniformRand();
  if (rndm < 0.0455) {
	particle_energy = 1039.65*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.5921) {
	particle_energy =  830.49*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.6011) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 830.49*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three96;}
  else if (rndm < 0.7588) {
	particle_energy =  707.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto five19;}
  else if (rndm < 0.7656) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 707.41*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto five19;}
  else if (rndm < 0.8838) {
	particle_energy =  204.03*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one023;}  
  else if (rndm < 0.8946) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 204.03*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}
  else if (rndm < 0.9124) {
	particle_energy =  135.54*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one091;}  
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 135.54*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one091;}

 one175: rndm = G4UniformRand();
  if (rndm < 0.0925) {
	particle_energy = 1175.31*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.3016) {
	particle_energy = 1117.64*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.7993) {
	particle_energy =  988.43*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto  one87;}
  else if (rndm < 0.8949) {
	particle_energy =  231.42*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine44;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 231.42*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine44;}

 one168: rndm = G4UniformRand();
  if (rndm < 0.0842) {
	particle_energy = 1110.61*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.3362) {
	particle_energy =  840.38*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.3397) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 840.38*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
      goto three28;}
  else if (rndm < 0.7430) {
	particle_energy =  772.29*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.7539) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 772.29*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three96;}
  else if (rndm < 0.7650) {
	particle_energy =  648.84*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto five19;}
  else if (rndm < 0.7652) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 648.82*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto five19;}  
  else if (rndm < 0.8523) {
	particle_energy =  199.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.8607) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 199.41*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.9042) {
	particle_energy =  145.85*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one023;}
  else if (rndm < 0.9913) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 145.85*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}
  else if (rndm < 0.9983) {
	particle_energy =   77.34*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one123;}  
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 77.34*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one123;}

 one153: rndm = G4UniformRand();
  if (rndm < 0.1348) {
	particle_energy = 1153.52*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.2596) {
	particle_energy = 1095.68*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.4792) {
	particle_energy =  321.65*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto eight32;}
  else if (rndm < 0.5099) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 321.65*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto eight32;}
  else if (rndm < 0.6946) {
	particle_energy =  278.95*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto eight74;}
  else if (rndm < 0.8147) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 278.95*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto eight74;}
  else if (rndm < 0.8821) {
	particle_energy =  278.95*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.9158) {
	particle_energy =  173.96*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine80;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 173.96*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
      goto nine80;}

 one123: rndm = G4UniformRand();
  if (rndm < 0.0218) {
	particle_energy = 1065.17*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.7182) {
	particle_energy =  794.95*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.7321) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 794.95*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three28;}
  else if (rndm < 0.8346) {
	particle_energy =  726.86*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.8366) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 726.86*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three96;}
  else if (rndm < 0.9549) {
	particle_energy =  153.98*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto nine69;}
  else if (rndm < 0.9751) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 153.98*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto nine69;}
  else if (rndm < 0.9904) {
	particle_energy =  100.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one023;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 100.41*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one023;}

one091: rndm = G4UniformRand();
  if (rndm < 0.2077) {
	particle_energy = 1033.25*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.9904) {
	particle_energy =  904.20*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 904.20*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one87;}

one023: rndm = G4UniformRand();
  if (rndm < 0.7443) {
	particle_energy =  964.77*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.7523) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 964.77*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty8;}
  else if (rndm < 0.9965) {
	particle_energy =  835.71*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 835.71*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
    goto one87;}

 one016: rndm = G4UniformRand();
  if (rndm < 0.0469) {
	particle_energy = 1016.44*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.7640) {
	particle_energy =  958.61*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.7962) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 958.61*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty8;}
  else if (rndm < 0.9983) {
	particle_energy =  620.37*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 620.37*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three96;}

 nine80: rndm = G4UniformRand();
  if (rndm < 0.1086) {
	particle_energy =  979.48*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.1693) {
	particle_energy =  921.98*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.5421) {
	particle_energy =  651.94*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else {
	particle_energy =  583.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}

 nine69: rndm = G4UniformRand();
  if (rndm < 0.3694) {
	particle_energy =  968.97*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.3733) {
	return;}
  else if (rndm < 0.9809) {
	particle_energy =  911.21*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.9882) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 911.21*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty8;}
  else if (rndm < 0.9996) {
	particle_energy =  782.14*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else  {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 782.14*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one87;}

 nine68: rndm = G4UniformRand();
  if (rndm < 0.1976) {
	particle_energy =  640.34*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else if (rndm < 0.2213) {
	goto three28;}
  else if (rndm < 0.7671) {
	particle_energy =  572.29*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}
  else if (rndm < 0.8162) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 572.29*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three96;}
  else if (rndm < 0.9903) {
	particle_energy =  449.15*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto five19;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 449.15*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto five19;}

 nine44: rndm = G4UniformRand();
  if (rndm < 0.5441) {
	particle_energy =  944.20*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	particle_energy = 616.22*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}

 eight74: rndm = G4UniformRand();
  if (rndm < 0.0856) {
	particle_energy =  874.44*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.1395) {
	particle_energy =  816.71*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8; }
  else if (rndm < 0.2605) {
	particle_energy =  688.10*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.6233) {
	particle_energy =  546.45*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28;}
  else {
	particle_energy =  478.30*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three96;}

 eight32: rndm = G4UniformRand();
  if (rndm < 0.2974) {
	particle_energy =  774.10*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.9913) {
	particle_energy =  503.82*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three28; }
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 503.82*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three28;}

 five19: rndm = G4UniformRand();
  if (rndm < 0.8335) {
	particle_energy =  332.37*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else if (rndm < 0.8584) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 332.37*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one87; }
  else if (rndm < 0.9745) {
	particle_energy =  141.02*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three78;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 141.02*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three78;}

 three96: rndm = G4UniformRand();
  if (rndm < 0.7125) {
	particle_energy =  338.32*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else if (rndm < 0.7330) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 338.32*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty8; }
  else if (rndm < 0.9789) {
	particle_energy =  209.25*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 186.82*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one87;}

 three78: rndm = G4UniformRand();
  if (rndm < 0.5883) {
	particle_energy =  191.35*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one87;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-")); 
		particle_energy = 209.25*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one87;}

 three28: rndm = G4UniformRand();
  if (rndm < 0.4449) {
	particle_energy =  327.99*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.4586) {
	return;}
  else if (rndm < 0.9757) {
	particle_energy =  270.24*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 270.24*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty8;}

 one87: rndm = G4UniformRand();
  if (rndm < 0.2079) {
	particle_energy =  129.07*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty8;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 129.07*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty8;}

 fifty8:rndm = G4UniformRand();
  if (rndm < 0.0064) {
	particle_energy =  129.07*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	return;}
}



int RadiogenicSources::GenerateBa133Decay(G4PrimaryVertex* vert)
{
  G4double rndm,rndm2,rndm3,rndm4,rndm5;
  G4int numberOfGamma = 0;

  SetParticleDefinition(particleTable->FindParticle("gamma"));

  rndm = G4UniformRand();
  if (rndm<=0.86)
    goto decay1;
  else
    goto decay2;

 decay1: rndm2 = G4UniformRand();
  if(rndm2>=0 && rndm2<=0.7165) {
	particle_energy = 356.02*keV; numberOfGamma +=1; CreatePrimaryParticle(vert);
	goto decay4;}
  else if(rndm2>0.7165 && rndm2<=0.7348) {
	goto decay4;}
  else if (rndm2>0.7348 && rndm2<=0.8175) {
	particle_energy = 276.40*keV; numberOfGamma +=1; CreatePrimaryParticle(vert);
	goto decay3;}
  else if (rndm2>0.8175 && rndm2<=0.8222) {
	goto decay3;}
  else if (rndm2>0.8222 && rndm2<=0.8476) {
	particle_energy = 53.16*keV; numberOfGamma +=1; CreatePrimaryParticle(vert);
	goto decay2;}
  else if (rndm2>0.8476 && rndm2<=1.0) {
	goto decay2;}

 decay2: rndm3 = G4UniformRand();
  if (rndm3<=0.3110) {
	particle_energy = 383.85*keV; numberOfGamma +=1; CreatePrimaryParticle(vert);
	return(numberOfGamma);}
  else if (rndm3<=0.3173) {
	return(numberOfGamma);}
  else if (rndm3<=0.9549) {
	particle_energy = 302.85*keV; numberOfGamma +=1; CreatePrimaryParticle(vert);
	goto decay4;}
  else if (rndm3<=0.9828) {
	goto decay4;}
  else if (rndm3<=0.9985) {
	particle_energy = 223.23*keV; numberOfGamma +=1; CreatePrimaryParticle(vert);
	goto decay3;}
  else {
	goto decay3;}

 decay3: rndm4 = G4UniformRand();
  if (rndm4<=0.0815) {
	particle_energy = 160.61*keV; GenerateIsotropicFlux(); numberOfGamma +=1; CreatePrimaryParticle(vert);
	return(numberOfGamma);}
  else if (rndm4<=0.1057) {
	return(numberOfGamma);}
  else if (rndm4<=0.4369) {
	particle_energy = 79.62*keV; GenerateIsotropicFlux(); numberOfGamma +=1; CreatePrimaryParticle(vert);
	goto decay4;}
  else {
	goto decay4;}

 decay4: rndm5 = G4UniformRand();
  if (rndm5<=0.3802) {
	particle_energy = 81.00*keV; GenerateIsotropicFlux(); numberOfGamma +=1; CreatePrimaryParticle(vert);
	return(numberOfGamma);}
  else {
	return(numberOfGamma);}
}



void RadiogenicSources::GenerateBi210Decay(G4PrimaryVertex* vert)
{
  if (ceFlag == 1) {
	SetParticleDefinition(particleTable->FindParticle("e-"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Bi210rank],nlfile[Bi210rank+1])*keV;
	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);}
  return;
}



void RadiogenicSources::GenerateBi212Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();
  
  SetParticleDefinition(particleTable->FindParticle("gamma"));

  if (rndm < 0.3594)
    goto alpha;    // Bi212 -> Tl208
  else
    goto beta;     // Bi212 -> Po212


 beta: rndm = G4UniformRand();
  if (rndm < 0.0101)
    goto one806;
  else if (rndm < 0.0106)
    goto one801;
  else if (rndm < 0.0144)
    goto one679;
  else if (rndm < 0.0439)
    goto one621;
  else if (rndm < 0.0666)
    goto one513;
  else if (rndm < 0.1367)
    goto seven27;
  else 
    return;

 one806: rndm = G4UniformRand();
  if (rndm < 0.17222) {
	particle_energy = 1806.00*keV; CreatePrimaryParticle(vert);
	return;}
  else {
	particle_energy = 1078.62*keV; CreatePrimaryParticle(vert);
      goto seven27; }
     
 one801: rndm = G4UniformRand();
  if (rndm < 0.3846) {
	particle_energy = 1800.20*keV; CreatePrimaryParticle(vert);
	return;}
  else {
	particle_energy = 1074.00*keV; CreatePrimaryParticle(vert);
	goto seven27; }

 one679: rndm = G4UniformRand();
  if (rndm < 0.2794) {
	particle_energy = 1679.50*keV; CreatePrimaryParticle(vert);
	return;}
  else {
	particle_energy =  952.10*keV; CreatePrimaryParticle(vert);
	goto seven27; }

 one621: rndm = G4UniformRand();
  if (rndm < 0.8001) {
	particle_energy = 1620.58*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9944) {
	particle_energy =  893.39*keV; CreatePrimaryParticle(vert);
	goto seven27; }
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 893.39*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto seven27;}

 one513: rndm = G4UniformRand();
  if (rndm < 0.2134) {
	particle_energy = 1512.75*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9690) {
	particle_energy =  785.42*keV; CreatePrimaryParticle(vert);
	goto seven27; }
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 785.42*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
      goto seven27;}
  
 seven27: rndm = G4UniformRand();
  if (rndm < 0.9861) {
	particle_energy =  727.18*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 727.18*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));} 
	return;}


 alpha: rndm = G4UniformRand();
  if (rndm < 0.2724)
    return;
  else if (rndm < 0.9723)
    goto fourty;
  else if (rndm < 0.9901)
    goto three28;
  else 
    goto four93;

 fourty: rndm = G4UniformRand();
  if (rndm < 0.9606) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy =  39.85*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else {
	particle_energy = 39.84*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}

 three28: rndm = G4UniformRand();
  if (rndm < 0.5052) {
	particle_energy = 288.07*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fourty;}
  else if (rndm < 0.7335) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 288.07*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fourty;}
  else if (rndm < 0.9355) {
	particle_energy = 327.96*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 327.96*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}

 four93: rndm = G4UniformRand();
  if (rndm < 0.1271) {
	particle_energy = 164.00*keV; CreatePrimaryParticle(vert);
	goto three28;}
  else {
	particle_energy = 452.80*keV; CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GenerateBi214Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();

  SetParticleDefinition(particleTable->FindParticle("gamma"));

  if (rndm < 0.0059)
    goto two728;
  else if (rndm < 0.0083)
    goto two719;
  else if (rndm < 0.0095)
    goto  two700;
  else if (rndm < 0.0113)
    goto  two695;
  else if (rndm < 0.0125)
    goto  two662;
  else if (rndm < 0.0137)
    goto  two508;
  else if (rndm < 0.0158)
    goto  two506;
  else if (rndm < 0.0303)
    goto  two482;
  else if (rndm < 0.0576)
    goto  two448;
  else if (rndm < 0.0633)
    goto  two293;
  else if (rndm < 0.0654)
    goto  two267;
  else if (rndm < 0.0675)
    goto  two209;
  else if (rndm < 0.1228)
    goto  two204;
  else if (rndm < 0.1311)
    goto  two193;
  else if (rndm < 0.1366)
    goto  two148;
  else if (rndm < 0.1780)
    goto  two119;
  else if (rndm < 0.1792)
    goto two089;
  else if (rndm < 0.2081)
    goto  two017;
  else if (rndm < 0.2247)
    goto  two010;
  else if (rndm < 0.2384)
    goto one994;
  else if (rndm < 0.2543)
    goto  one890;
  else if (rndm < 0.3367)
    goto  one847;
  else if (rndm < 0.5054)
    goto one764;
  else if (rndm < 0.5076)
    goto  one743;
  else if (rndm < 0.6822)
    goto  one730;
  else if (rndm < 0.6886)
    goto one661;
  else if (rndm < 0.7191)
    goto  one543;
  else if (rndm < 0.7278)
    goto  one415;
  else if (rndm < 0.7994)
    goto one378;
  else if (rndm < 0.8014)
    goto  six09;
  else 
    return;

 two728: rndm = G4UniformRand();
  if (rndm < 0.0131) {
	particle_energy = 2119.93*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.0213) {
	particle_energy = 1351.00*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.0735) {
	particle_energy = 1066.90*keV; CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.7825) {
	particle_energy = 964.08*keV; CreatePrimaryParticle(vert);
	goto one764;}
  else if (rndm < 0.8552) {
	particle_energy = 733.65*keV; CreatePrimaryParticle(vert);
	goto one994;}
  else if (rndm < 0.8832) {
	particle_energy = 525.00*keV; CreatePrimaryParticle(vert);
	goto two204;}
  else if (rndm < 0.8862) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 525.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two204;}
  else  {
	particle_energy = 280.94*keV; CreatePrimaryParticle(vert);
	goto two448;}
    
 two719: rndm = G4UniformRand();
  if (rndm < 0.0067) {
	particle_energy = 2719.40*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.3118) {
	particle_energy = 2109.90*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.3936) {
	particle_energy = 1341.50*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.8028) {
	particle_energy =1303.76*keV; CreatePrimaryParticle(vert);
	goto one415;}
  else if (rndm < 0.8400) {
	particle_energy =1058.10*keV; CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.9219) {
	particle_energy = 976.20*keV; CreatePrimaryParticle(vert);
	goto one743;}
  else if (rndm < 0.9702) {
	particle_energy = 708.70*keV; CreatePrimaryParticle(vert);
	goto two010;}
  else {
	particle_energy = 600.00*keV; CreatePrimaryParticle(vert);
	goto two119;}

 two700: rndm = G4UniformRand();
  if (rndm < 0.0188) {
	particle_energy = 2699.94*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.4881) {
	particle_energy = 2109.90*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.5734) {
	particle_energy = 1321.50*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.7014) {
	particle_energy =1155.60*keV; CreatePrimaryParticle(vert);
	goto one543;}
  else if (rndm < 0.7952) {
	particle_energy =1038.00*keV; CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.8805) {
	particle_energy = 934.50*keV; CreatePrimaryParticle(vert);
	goto one764;}
  else if (rndm < 0.9317) {
	particle_energy = 687.70*keV; CreatePrimaryParticle(vert);
	goto two010;}
  else {
	particle_energy = 600.00*keV; CreatePrimaryParticle(vert);
	goto two204;}

 two695: rndm = G4UniformRand();
  if (rndm < 0.1990) {
	particle_energy = 2694.80*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.2482) {
	particle_energy = 2085.00*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.2758) {
	particle_energy = 1419.70*keV; CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.7678) {
	particle_energy =1316.96*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.8947) {
	particle_energy =1033.20*keV; CreatePrimaryParticle(vert);
	goto one661;}
  else {
	particle_energy = 930.50*keV; CreatePrimaryParticle(vert);
	goto one764;}

 two662: rndm = G4UniformRand();
  if (rndm < 0.0023) {
	particle_energy = 2662.00*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.6056) {
	particle_energy = 2052.94*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.6906) {
	particle_energy = 1284.00*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else {
	particle_energy =1118.90*keV; CreatePrimaryParticle(vert);
	goto one543;}

 two508: rndm = G4UniformRand();
  if (rndm < 0.5304) {
	particle_energy = 1898.70*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9130) {
	particle_energy = 1130.80*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else {
	particle_energy = 965.00*keV; CreatePrimaryParticle(vert);
	goto one543;}

 two506: rndm = G4UniformRand();
  if (rndm < 0.0323) {
	particle_energy = 2505.56*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.8507) {
	particle_energy = 1896.30*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9519) {
	particle_energy = 1230.50*keV; CreatePrimaryParticle(vert);
      goto one378;}
  else {
	particle_energy = 962.00*keV; CreatePrimaryParticle(vert);
	goto one543;}

 two482: rndm = G4UniformRand();
  if (rndm < 0.0018) {
	particle_energy = 2482.80*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.2187) {
	particle_energy = 1873.20*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.6437) {
	particle_energy = 1207.68*keV; CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.7131) {
	particle_energy = 1104.80*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.8693) {
	particle_energy =  821.18*keV; CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.8750) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 821.18*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one661;}
  else if (rndm < 0.9965) {
	particle_energy =  752.84*keV; CreatePrimaryParticle(vert);
	goto one730;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 752.84;CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one730;}

 two448: rndm = G4UniformRand();
  if (rndm < 0.5569) {
	particle_energy = 2447.76*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.7054) {
	particle_energy = 1838.36*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.7239) {
	particle_energy = 1173.10*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.7240) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1173.10*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one275;}
  else if (rndm < 0.8280) {
	particle_energy = 1070.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.8282) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1070.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.8549) {
	particle_energy = 1032.40*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one415;}
  else if (rndm < 0.8550) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1032.40*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one415;}
  else if (rndm < 0.8880) {
	particle_energy =  904.20*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one543;}
  else if (rndm < 0.8882) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 904.20*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one543;}
  else if (rndm < 0.9995) {
	particle_energy =  786.10*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one661;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 786.10*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one661;}

 two293: rndm = G4UniformRand();
  if (rndm < 0.5253) {
	particle_energy = 2293.36*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9616) {
	particle_energy = 1683.99*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else {
	particle_energy = 915.80*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}

two267: rndm = G4UniformRand();
  if (rndm < 0.0948) {
	particle_energy = 2266.67*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.4265) {
	particle_energy = 1657.40*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.6351) {
	particle_energy = 1657.40*keV; CreatePrimaryParticle(vert);
	goto one543;}
  else if (rndm < 0.8957) {
	particle_energy = 536.90*keV; CreatePrimaryParticle(vert);
	goto one730;}
  else if (rndm < 0.9763) {
	particle_energy = 502.20*keV; CreatePrimaryParticle(vert);
	goto one764;}
  else {
	particle_energy = 376.60*keV; CreatePrimaryParticle(vert);
	goto one890;}

 two209: rndm = G4UniformRand();
  if (rndm < 0.8225) {
	particle_energy = 1599.31*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9307) {
	particle_energy =  934.00*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else {
	particle_energy = 547.10*keV; CreatePrimaryParticle(vert);
	goto one661;}

 two204: rndm = G4UniformRand();
  if (rndm < 0.8732) {
	particle_energy = 2204.21*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9289) {
	particle_energy = 1594.73*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9451) {
	particle_energy =  826.12*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.9456) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 816.12*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one275;}
  else if (rndm < 0.9478) {
	particle_energy =  789.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9479) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 789.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.9554) {
	particle_energy =  661.40*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one415;}
  else if (rndm < 0.9557) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 661.40*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one415;}
  else if (rndm < 0.9719) {
	particle_energy =  543.40*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one543;}
  else if (rndm < 0.9731) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 543.40*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one543;}
  else if (rndm < 0.9910) {
	particle_energy =  474.38*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.9928) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 474.38*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one661;}
  else {
	particle_energy =  460.90*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one661;}

 two193: rndm = G4UniformRand();
  if (rndm < 0.0365) {
	particle_energy = 2192.60*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.8878) {
	particle_energy = 1583.22*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9146) {
	particle_energy =  815.08*keV; CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.9152) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 815.08*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one275;}
  else if (rndm < 0.9821) {
	particle_energy =  649.18*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9850) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 649.18*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
   else if (rndm < 0.9984) {
	particle_energy =  428.00*keV; CreatePrimaryParticle(vert);
	goto one415;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 428.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one415;}

two148: rndm = G4UniformRand();
  if (rndm < 0.0266) {
	particle_energy = 2147.80*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.8738) {
	particle_energy = 1538.50*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9003) {
	particle_energy =  873.00*keV; CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.9502) {
	particle_energy = 769.70*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else {
	particle_energy = 486.30*keV; CreatePrimaryParticle(vert);
	goto one415;}

 two119: rndm = G4UniformRand();
  if (rndm < 0.2756) {
	particle_energy = 2118.55*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.7849) {
	particle_energy = 1509.23*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.7945) {
	particle_energy =  741.50*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.8003) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 741.50*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.8828) {
	particle_energy =  703.11*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one415;}
  else if (rndm < 0.8877) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 703.11*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one415;}
  else if (rndm < 0.9768) {
	particle_energy =  389.10*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one730;}
   else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 389.10*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one730;}

 two089: rndm = G4UniformRand();
  if (rndm < 0.4331) {
	particle_energy = 1479.22*keV; CreatePrimaryParticle(vert);
	return;}
  else {
	particle_energy = 710.80*keV; CreatePrimaryParticle(vert);
	goto six09;}

 two017: rndm = G4UniformRand();
  if (rndm < 0.9799) {
	particle_energy = 1407.98*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9836152) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1407.98*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9952) {
	particle_energy = 639.37*keV; CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9954) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 639.37*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.9979) {
	particle_energy =  356.50*keV; CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.9981) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 356.50*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one661;}
  else if (rndm < 0.9991) {
	particle_energy =  253.00*keV; CreatePrimaryParticle(vert);
	goto one730;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 253.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one730;}

two010: rndm = G4UniformRand();
  if (rndm < 0.0298) {
	particle_energy = 2010.71*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9542) {
	particle_energy = 1401.50*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9592) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1401.50*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9920) {
	particle_energy = 633.14*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9592) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 633.14*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else {
	particle_energy = 596.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one415;}

 one994: rndm = G4UniformRand();
  if (rndm < 0.6248) {
	particle_energy = 1385.31*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9198) {
	particle_energy =  719.86*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.9241) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 719.86*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
      goto one275;}
  else if (rndm < 0.9472) {
	particle_energy =  617.10*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9474) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 617.10*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.9966) {
	particle_energy =  333.61*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one661;}
  else if (rndm < 0.9978) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 333.61*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one661;}
  else if (rndm < 0.9999) {
	particle_energy =  230.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one730;}
   else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 230.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one730;}

 one890: rndm = G4UniformRand();
  if (rndm < 0.0566) {
	particle_energy = 1890.35*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9616) {
	particle_energy = 1280.96*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9722) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1280.96*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9998) {
	particle_energy =  615.78*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one275;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 615.78*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one275;}

one847: rndm = G4UniformRand();
  if (rndm < 0.2466) {
	particle_energy = 1847.42*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9552) {
	particle_energy = 1238.11*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9642) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 918.97*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9675) {
	particle_energy =  832.35*keV; CreatePrimaryParticle(vert);
	goto one015;}
  else if (rndm < 0.9762) {
	particle_energy =  572.83*keV; CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.9931) {
	particle_energy =  469.69*keV;CreatePrimaryParticle(vert);
	goto one378;}  
  else if (rndm < 0.9947) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 469.69*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.9987) {
	particle_energy =  394.43*keV; CreatePrimaryParticle(vert);
	goto one543;}
   else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 394.43*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one543;}

 one764: rndm = G4UniformRand();
  if (rndm < 0.8791) {
	particle_energy = 1764.50*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return ;}
  else if (rndm < 0.9729) {
	particle_energy = 1155.19*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9743) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1155.19*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9909) {
	particle_energy =  387.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else if (rndm < 0.9935) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 387.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}
  else if (rndm < 0.9981) {
	particle_energy =  348.80*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one415;}
  else if (rndm < 0.9997) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 348.80*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one415;}
  else if (rndm < 0.9999) {
	particle_energy =  221.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one543;}
   else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 221.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one543;}

 one743: rndm = G4UniformRand();
  if (rndm < 0.8485) {
	particle_energy = 1133.66*keV; CreatePrimaryParticle(vert);
	goto six09;}
  else {
	particle_energy = 727.80*keV; CreatePrimaryParticle(vert);
	goto one015;}

one730: rndm = G4UniformRand();
  if (rndm < 0.1575) {
	particle_energy = 1729.60*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9667) {
	particle_energy = 1120.29*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9798) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1120.29*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9952) {
	particle_energy =  454.77*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one275;}
  else if (rndm < 0.9954) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 454.77*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one275;}
  else if (rndm < 0.9992) {
	particle_energy =  351.90*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one378;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 351.90*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one378;}

 one661: rndm = G4UniformRand();
  if (rndm < 0.7703) {
	particle_energy = 1661.28*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	particle_energy =1051.96*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}

 one543: rndm = G4UniformRand();
  if (rndm < 0.0951) {
	particle_energy = 1543.32*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9682) {
	particle_energy =  934.06*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else if (rndm < 0.9899) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 934.06*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
  else if (rndm < 0.9910) {
	particle_energy =  528.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one015;}
  else if (rndm < 0.9996) {
	particle_energy =  269.00*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one275;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 269.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one275;}

 one415: rndm = G4UniformRand();
  if (rndm < 0.2919) {
	particle_energy = 1415.80*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9920) {
	particle_energy =  806.17*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 806.17*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}

 one378: rndm = G4UniformRand();
  if (rndm < 0.4448) {
	particle_energy = 1377.70*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.4466) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 1377.70*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else if (rndm < 0.9912) {
	particle_energy =  768.36*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 768.36*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}

 one275: rndm = G4UniformRand();
  if (rndm < 0.9943) {
	particle_energy =  665.45*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto six09;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 665.45*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto six09;}
	
 one015: 
  particle_energy =  665.45*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
  goto six09;

 six09: rndm = G4UniformRand();
  if (rndm < 0.9798) {
	particle_energy =  609.31*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 609.31*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
}



void RadiogenicSources::GenerateC14Decay(G4PrimaryVertex* vert)
{
  if (ceFlag == 1) {
	SetParticleDefinition(particleTable->FindParticle("e-"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[C14rank],nlfile[C14rank+1])*keV;
	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);}
  return;
}

void RadiogenicSources::GenerateCf252Decay(G4PrimaryVertex* vert)
{  G4double rndm;
  if (ceFlag == 1) {

       rndm=G4UniformRand();
    if(rndm<0.89){
//	SetParticleDefinition(particleTable->FindParticle("neutron"));
//	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
//	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//        SetParticleDefinition(particleTable->FindParticle("gamma"));
//        particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        SetParticleDefinition(particleTable->FindParticle("gamma"));
//        particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//       return;}
//      else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
//     else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
//     else if() {
         SetParticleDefinition(particleTable->FindParticle("neutron"));
 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
        SetParticleDefinition(particleTable->FindParticle("neutron"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
        SetParticleDefinition(particleTable->FindParticle("neutron"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
          SetParticleDefinition(particleTable->FindParticle("neutron"));
 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
        SetParticleDefinition(particleTable->FindParticle("neutron"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("neutron"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("gamma"));
         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("gamma"));
         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
        return;}
//     else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
//     else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
// 
//     else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
//    else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
//     else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
//      else if() {
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("neutron"));
// 	particle_energy = BetaEnergyFromSpectrum(nlfile[Cf252rank],nlfile[Cf252rank+1])*keV;
// 	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//         SetParticleDefinition(particleTable->FindParticle("gamma"));
//         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
//        return;}
     else{
         SetParticleDefinition(particleTable->FindParticle("gamma"));
         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
         SetParticleDefinition(particleTable->FindParticle("gamma"));
         particle_energy = 100.00*keV; GenerateIsotropicFlux();CreatePrimaryParticle(vert);
        return;}
 }
}




void RadiogenicSources::GenerateCo57Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();

  SetParticleDefinition(particleTable->FindParticle("gamma"));

  if (rndm < 0.9982)
    goto one36;
  else 
    goto seven07;

 seven07: rndm = G4UniformRand();
  if (rndm < 0.0259) {
	particle_energy =  706.54*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;} 
  else if (rndm < 0.9115) {
	particle_energy =  692.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto ten4;} 
  else if (rndm < 0.9865) {
	particle_energy =  570.09*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one36;} 
  else {
	particle_energy =  339.69*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three67;} 

 three67: rndm = G4UniformRand();
  if (rndm < 0.1381) {
	particle_energy = 366.80*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;} 
  else if (rndm < 0.9205) {
	particle_energy =  352.33*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto ten4;} 
  else {
	particle_energy = 230.40*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto one36;} 

 one36: rndm = G4UniformRand();
  if (rndm < 0.1036) {
	particle_energy = 136.47*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;} 
  else if (rndm < 0.1177) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 136.47*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;} 
  else if (rndm < 0.9793) {
	particle_energy = 122.06*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto ten4;} 
 else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 122.06*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto ten4;}

 ten4: rndm = G4UniformRand();
  if (rndm < 0.1098) {
	particle_energy = 14.41*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;} 
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 14.41*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
}



void RadiogenicSources::GenerateCo60Decay(G4PrimaryVertex* vert)
{
  SetParticleDefinition(particleTable->FindParticle("gamma"));

  particle_energy = 1173.24*keV; CreatePrimaryParticle(vert);
  particle_energy = 1332.50*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
  return;
}



void RadiogenicSources::GenerateCs137Decay(G4PrimaryVertex* vert)
{ 
  G4double rndm;
//
//  rndm = G4UniformRand();
//  if (rndm < 0.557)
//
  if (ceFlag == 1) {
	SetParticleDefinition(particleTable->FindParticle("e-"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Cs137rank],nlfile[Cs137rank+1])*keV;
	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);

	rndm = G4UniformRand();
	if (rndm < 0.1097) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 661.66*keV; CreatePrimaryParticle(vert);}
	else if (rndm < 0.851) {
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 661.66*keV; CreatePrimaryParticle(vert);}
	}
  else {
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 661.66*keV; CreatePrimaryParticle(vert);}
  return;
 }



void RadiogenicSources::GenerateGa68Decay(G4PrimaryVertex* vert)
{
  if (ceFlag == 1) {
	SetParticleDefinition(particleTable->FindParticle("e+"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Ga68rank],nlfile[Ga68rank+1])*keV;
	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);}
  return;
}



void RadiogenicSources::GenerateK40Decay(G4PrimaryVertex* vert)
{
  G4double r_bd;   // random number to select the beta decay path (K40 -> Ca40 or K40 -> Ar40)
   
  if (ceFlag == 1) {
	GenerateIsotropicFlux(); 
	r_bd = G4UniformRand();
	if (r_bd < 0.893) {  // K40 -> Ca40
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = BetaEnergyFromSpectrum(nlfile[K40rank1],nlfile[K40rank1+1])*keV;
		CreatePrimaryParticle(vert);}
	else {  // K40 -> Ca40
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 1460.81*keV;
		CreatePrimaryParticle(vert);
		if (r_bd > (1-0.001)) {
			SetParticleDefinition(particleTable->FindParticle("e+"));
			particle_energy = BetaEnergyFromSpectrum(nlfile[K40rank2],nlfile[K40rank2+1])*keV;
			CreatePrimaryParticle(vert);}
		}
	}
  else {	
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 1460.81*keV;
	CreatePrimaryParticle(vert);}
  return;
}



void RadiogenicSources::GeneratePb210Decay(G4PrimaryVertex* vert)
{
  G4double r_bd; // r_a;   // random number to select the beta decay path (Pb210 -> Bi210 or Pb210 -> Bi210*+CE ...) ;  alpha
   
  if (ceFlag == 1) {
	// emission of the beta (100% of the time)
	SetParticleDefinition(particleTable->FindParticle("e-"));   
	r_bd = G4UniformRand();
	if (r_bd < 0.84) { // Pb210 -> Bi210*
		GenerateIsotropicFlux();
		particle_energy = BetaEnergyFromSpectrum(nlfile[Pb210rank1],nlfile[Pb210rank1+1])*keV;
		CreatePrimaryParticle(vert);
		//G4cout << "energy of the beta" << particle_energy << G4endl;
		if (r_bd < 0.05) {
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("gamma"));
			particle_energy = 46.5*keV; CreatePrimaryParticle(vert);} // gamma from Bi210* nucleus
		else if (r_bd < 0.176) {
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("e-"));
			particle_energy = 31.00*keV; CreatePrimaryParticle(vert); // conversion electron
			GenerateIsotropicFlux();
			particle_energy = 15.53*keV; CreatePrimaryParticle(vert);} // auger electron
		else if (r_bd < 0.6462) {
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("e-"));
			particle_energy = 31.00*keV; CreatePrimaryParticle(vert); // conversion electron
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("gamma"));
			particle_energy = 15.53*keV; CreatePrimaryParticle(vert);} // fluorescence
		else if (r_bd < 0.7842) {
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("e-"));
			particle_energy = 43.00*keV; CreatePrimaryParticle(vert); // conversion electron
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("gamma"));
			particle_energy = 3.53*keV; CreatePrimaryParticle(vert);} // fluorescence
		else if (r_bd < 0.84) {
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("e-"));
			particle_energy = 45.00*keV; CreatePrimaryParticle(vert); // conversion electron
			GenerateIsotropicFlux();
			SetParticleDefinition(particleTable->FindParticle("gamma"));
			particle_energy = 1.53*keV; CreatePrimaryParticle(vert);} // fluorescence
			}
	else { // Pb210 -> Bi210
		GenerateIsotropicFlux();
		particle_energy = BetaEnergyFromSpectrum(nlfile[Pb210rank2],nlfile[Pb210rank2+1])*keV;
		CreatePrimaryParticle(vert);}
	}
  else { //ceFlag ~= 1
	GenerateIsotropicFlux();
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 46.5*keV;
	CreatePrimaryParticle(vert);}
  return;
}



void RadiogenicSources::GeneratePb212Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();

  SetParticleDefinition(particleTable->FindParticle("gamma"));
   
  if (rndm < 0.1200)
    return;
  else if (rndm < 0.17)
    goto four15;
  else 
    goto two39;

 two39: rndm = G4UniformRand();
  if (rndm < 0.4762) {
	return;}
  else {
	particle_energy = 238.63*keV; CreatePrimaryParticle(vert);
	return;}

 four15: rndm = G4UniformRand();
  if (rndm < 0.3126) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 300.09*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto one15;}
  else if (rndm < 0.9849) {
	particle_energy = 300.09*keV; CreatePrimaryParticle(vert);
	goto one15;}
  else if (rndm < 0.9952 ) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 176.66*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two39;}
  else {
	particle_energy = 176.66*keV; CreatePrimaryParticle(vert);
	goto two39;}

 one15: rndm = G4UniformRand();
  if (rndm < 0.8768) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 115.18*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else {
	particle_energy = 115.18*keV; CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GeneratePb214Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();

  SetParticleDefinition(particleTable->FindParticle("gamma"));

  if (rndm < 0.0237)
    goto eight39;
  else if (rndm < 0.0342)
    goto five33;
  else if (rndm < 0.4979)
    goto  three52;
  else if (rndm < 0.9062)
    goto  two95;
  else 
    return;

 eight39: rndm = G4UniformRand();
  if (rndm < 0.2676) {
	particle_energy = 839.03*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.2686) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 839.03*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else if (rndm < 0.6296) {
	particle_energy = 785.91*keV; CreatePrimaryParticle(vert);
	goto fifty3;}
  else if (rndm < 0.6311) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 785.91*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty3;}
  else if (rndm < 0.7798) {
	particle_energy = 580.15*keV; CreatePrimaryParticle(vert);
	goto two58;}
  else if (rndm < 0.7809) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 580.15*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two58;}
  else if (rndm < 0.7979) {
	particle_energy = 544.00*keV; CreatePrimaryParticle(vert);
	goto two95;}
  else if (rndm < 0.7980) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 544.00*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two95;}
  else if (rndm < 0.9849) {
	particle_energy = 487.08*keV; CreatePrimaryParticle(vert);
	goto three52;}
  else if (rndm < 0.9869) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 487.08*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three52;}
  else if (rndm < 0.9996) {
	particle_energy = 305.40*keV; CreatePrimaryParticle(vert);
	goto five33;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 305.40*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto five33;}

 five33: rndm = G4UniformRand();
  if (rndm < 0.1775) {
	particle_energy = 533.69*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.1890) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 533.69*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else if (rndm < 0.4972) {
	particle_energy = 480.42*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty3;}
  else if (rndm < 0.5234) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 480.42*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty3;}
  else if (rndm < 0.8504) {
	particle_energy = 274.53*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto two58;}
  else if (rndm < 0.9779) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 274.53*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two58;}
  else if (rndm < 0.9919) {
	particle_energy = 238.40*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto two95;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 238.40*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two95;}

 three52: rndm = G4UniformRand();
  if (rndm < 0.7706) {
	particle_energy = 351.92*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9994) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 351.92*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else if (rndm < 0.9999) {
	particle_energy = 298.76*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty3;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 298.76*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty3;}

 two95: rndm = G4UniformRand();
  if (rndm < 0.4554) {
	particle_energy = 295.21*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.6695) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 295.21*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else if (rndm < 0.8541) {
	particle_energy = 241.98*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty3;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 241.98*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto fifty3;}

 two58: rndm = G4UniformRand();
  if (rndm < 0.5622) {
	particle_energy = 258.79*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9709) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 258.79*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
  else if (rndm < 0.9862) {
	particle_energy = 205.59*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto fifty3;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 205.59*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}

 fifty3: rndm = G4UniformRand();
  if (rndm < 0.0728) {
	particle_energy = 53.21*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy =  53.21*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
}



void RadiogenicSources::GeneratePo210Decay(G4PrimaryVertex* vert)
{
  SetParticleDefinition(particleTable->FindParticle("alpha"));
  //particle_energy = 5304.4*keV;
  particle_energy = 5389.2*keV;
  CreatePrimaryParticle(vert);
  return;
}

void RadiogenicSources::GeneratePo210RNDecay(G4PrimaryVertex* vert)
{

  GenerateIsotropicFlux();

  SetParticleDefinition(particleTable->FindParticle("alpha"));
  particle_energy = 5389.2*keV;
  CreatePrimaryParticle(vert);

// Pb206

	G4int AtomicNum = 82;
	G4int AtomicMass = 206;
	G4double excit = 0.0;
	
	G4double px = particle_momentum_direction.x();
	G4double py = particle_momentum_direction.y();
	G4double pz = particle_momentum_direction.z();
	
	particle_momentum_direction.setX(-px);
	particle_momentum_direction.setY(-py);
	particle_momentum_direction.setZ(-pz);

    	
  SetParticleDefinition(particleTable->GetIon(AtomicNum , AtomicMass, excit));
  particle_energy = 104.7*keV;
  CreatePrimaryParticle(vert);

  return;
}



void RadiogenicSources::GeneratePo212Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();

  SetParticleDefinition(particleTable->FindParticle("gamma"));

  if (rndm < 0.3594)
    goto Tl208;
  else
    goto Po208;

 Tl208: rndm = G4UniformRand();
  if (rndm < 0.0023)
    goto four180;
  else if (rndm < 0.0040)
    goto four125; 
  else if (rndm < 0.0352) 
    goto three961;
  else if (rndm < 0.0415)
    goto three920;
  else if (rndm < 0.2887)
    goto three708;
  else if (rndm < 0.5087)
    goto three475;
  else
    goto three198;

Po208: rndm = G4UniformRand();
  if (rndm < 0.97)
    return;
  else if (rndm < 0.99)
    goto three198;
  else
    goto two615;

 four180: rndm = G4UniformRand();
  if (rndm < 0.9031) {
	particle_energy = 982.7*keV; CreatePrimaryParticle(vert);
	goto three198;}
  else {
	particle_energy = 705.2*keV; CreatePrimaryParticle(vert);
	goto three475;}

 four125: rndm = G4UniformRand();
  if (rndm < 0.7857) {
	particle_energy = 927.6*keV; CreatePrimaryParticle(vert);
	goto three198;}
  else {
	particle_energy = 650.1*keV; CreatePrimaryParticle(vert);
	goto three475;}

 three961: rndm = G4UniformRand();
  if (rndm < 0.5815) {
	particle_energy = 763.13*keV; CreatePrimaryParticle(vert);
	goto three198;}
  else if (rndm < 0.6033) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 763.13*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three198;}
  else if (rndm < 0.6192) {
	particle_energy = 485.95*keV; CreatePrimaryParticle(vert);
	goto three475;}
  else if (rndm < 0.8416) {
	particle_energy = 252.61*keV; CreatePrimaryParticle(vert);
	goto three708;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 252.61*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three708;}

 three920: rndm = G4UniformRand();
  if (rndm < 0.3373) {
	particle_energy = 722.04*keV; CreatePrimaryParticle(vert);
	goto three198;}
  else if (rndm < 0.3511) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 722.04*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
      goto three198;}
  else if (rndm < 0.6501) {
	particle_energy = 211.40*keV; CreatePrimaryParticle(vert);
	goto three708;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 211.40*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three708;}

 three708: rndm = G4UniformRand();
  if (rndm < 0.1525) {
	particle_energy = 1093.90*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto two615;}
  else if (rndm < 0.8857) {
	particle_energy =  510.77*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three198;}
  else if (rndm < 0.9777) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 510.77*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three198;}
  else if (rndm < 0.9895) {
	particle_energy = 233.36*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three475;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 233.36*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three475;}

 three475: rndm = G4UniformRand();
  if (rndm < 0.5510) {
	particle_energy = 860.56*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto two615;}
  else if (rndm < 0.5662) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 860.56*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two615;}
  else if (rndm < 0.8459) {
	particle_energy = 277.36*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto three198;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 277.36*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto three198;}

 three198: rndm = G4UniformRand();
  if (rndm < 0.9797) {
	particle_energy = 583.19*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	goto two615;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 583.19*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto two615;}

 two615: rndm = G4UniformRand();
  if (rndm < 0.9157) {
	particle_energy = 2614.55*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
  else {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 2614.55*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	return;}
}



void RadiogenicSources::GenerateRa224Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
//
//  if (rndm < 0.9500)
//
  if (ceFlag == 1) {
	rndm = G4UniformRand();
	if (rndm < 0.2188) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 240.99*keV; CreatePrimaryParticle(vert);
		return;}
	else {
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 240.99*keV; CreatePrimaryParticle(vert);
		return;}
	}
  else {
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 240.99*keV; CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GenerateRa226Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
//
//  if (rndm < 0.9445)
//
  if (ceFlag == 1) {
	rndm = G4UniformRand();
	if (rndm < 0.3676) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 186.55*keV; CreatePrimaryParticle(vert);
		return;}
	else {
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 186.16*keV; CreatePrimaryParticle(vert);
		return;}
	}
  else {
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 186.16*keV; CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GenerateSb125Decay(G4PrimaryVertex* vert)
{
  if (ceFlag == 1) {
	SetParticleDefinition(particleTable->FindParticle("e-"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Sb125rank],nlfile[Sb125rank+1])*keV;
	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);}
  return;
}



void RadiogenicSources::GenerateTh228Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
//
//  if (rndm < 0.733)
//
  if (ceFlag == 1) {
	rndm = G4UniformRand();
	if (rndm < 0.9558) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 84.37*keV; CreatePrimaryParticle(vert);
		return;}
	else {
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 84.37*keV; CreatePrimaryParticle(vert);
		return;}
	}
  else {
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 84.37*keV; CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GenerateTh230Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
//
//  if (rndm < 0.763)
//
  if (ceFlag == 1) {
	rndm = G4UniformRand();
	if (rndm < 0.996) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 67.73*keV; CreatePrimaryParticle(vert);
		return;}
	else {
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 67.73*keV; CreatePrimaryParticle(vert);
		return;}
	}
  else {
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 67.73*keV; CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GenerateTh234Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
  rndm = G4UniformRand();

  SetParticleDefinition(particleTable->FindParticle("gamma"));

  if (rndm < 0.0274)
    goto ninety3;
  else if (rndm < 0.0983)
    goto seven3;
  else if (rndm < 0.2705)
    goto seven2;
  else 
    return;

 ninety3: rndm = G4UniformRand();
  if (rndm < 0.1517) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 112.81*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
		return;}
  else if (rndm < 0.7836) {
	particle_energy = 112.81*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.8195) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 83.30*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto thirty;}
  else {
	particle_energy = 83.30*keV; CreatePrimaryParticle(vert);
	goto thirty;}
  
 seven3: rndm = G4UniformRand();
  if (rndm < 0.0414) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 92.80*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
		return;}
  else if (rndm < 0.3182) {
	particle_energy = 92.8*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.5168) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 63.29*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto thirty;}
  else {
	particle_energy = 63.29*keV; CreatePrimaryParticle(vert);
	goto thirty;}

 seven2: rndm = G4UniformRand();
  if (rndm < 0.8228) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 92.38*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
		return;}
  else if (rndm < 0.9692) {
	particle_energy = 92.38*keV; CreatePrimaryParticle(vert);
	return;}
  else if (rndm < 0.9989) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 62.83*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
	goto thirty;}
  else {
	particle_energy = 62.83*keV; CreatePrimaryParticle(vert);
	goto thirty;}

 thirty: rndm = G4UniformRand();
  if (rndm < 0.1517) {
	if (ceFlag == 1) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 29.49*keV; CreatePrimaryParticle(vert);
		SetParticleDefinition(particleTable->FindParticle("gamma"));}
		return;}
  else {
	particle_energy = 29.49*keV; GenerateIsotropicFlux(); CreatePrimaryParticle(vert);
	return;}
}



void RadiogenicSources::GenerateTl204Decay(G4PrimaryVertex* vert)
{
  if (ceFlag == 1) {
	SetParticleDefinition(particleTable->FindParticle("e-"));
	particle_energy = BetaEnergyFromSpectrum(nlfile[Tl204rank],nlfile[Tl204rank+1])*keV;
	GenerateIsotropicFlux(); CreatePrimaryParticle(vert);}
  return;
}



void RadiogenicSources::GenerateU234Decay(G4PrimaryVertex* vert)
{
  G4double rndm;
//
//  if (rndm < 0.725)
//
  if (ceFlag == 1) {
	rndm = G4UniformRand();
	if (rndm < 0.997) {
		SetParticleDefinition(particleTable->FindParticle("e-"));
		particle_energy = 52.30*keV; CreatePrimaryParticle(vert);
		return;}
	else {
		SetParticleDefinition(particleTable->FindParticle("gamma"));
		particle_energy = 52.30*keV; CreatePrimaryParticle(vert);
		return;}
	}
  else {
	SetParticleDefinition(particleTable->FindParticle("gamma"));
	particle_energy = 52.30*keV; CreatePrimaryParticle(vert);
	return;}
}


void RadiogenicSources::GeneratePrimaries(G4Event *evt)
{
  if(particle_definition==NULL) {
    G4cout << "No particle has been defined!" << G4endl;
    return;}
  
  // Position
  G4bool srcconf = false;
  G4int LoopCount = 0;

  while(srcconf == false)  {
    if(SourcePosType == "Point")
      GeneratePointSource();
    else if(SourcePosType == "Volume")
      GeneratePointsInVolume();
/*    else if(SourcePosType == "PointRandomZIPs")
      // generate random position in a given detector material (defined by variable "myMaterial")
      SetRandomPosInMaterial_ZIPs();
    else if(SourcePosType == "PointRandomCopperTower")
      // generate random position in a given detector material (defined by variable "myMaterial")
      SetRandomPosInMaterial_CT();
    else if(SourcePosType == "PointRandomInnerShieldAir")
      // generate random position in a given detector material (defined by variable "myMaterial")
      SetRandomPosInMaterial_SA();
    else if(SourcePosType == "PointRandomCopperCans")
      // generate random position in a given detector material (defined by variable "myMaterial")
      SetRandomPosInMaterial_CC(); */
    else {
      G4cout << "Error: SourcePosType undefined" << G4endl;
      G4cout << "Generating point source" << G4endl;
      GeneratePointSource();
    }
    if(Confine == true) {
      srcconf = IsSourceConfined();
      // if source in confined srcconf = true terminating the loop
      // if source isnt confined srcconf = false and loop continues
    }
    else if(Confine == false)
      srcconf = true; // terminate loop
    
    LoopCount++;
    if(LoopCount == 100000) {
      G4cout << "*************************************" << G4endl;
        G4cout << "LoopCount = 100000" << G4endl;
        G4cout << "Either the source distribution >> confinement" << G4endl;
        G4cout << "or any confining volume may not overlap with" << G4endl;
        G4cout << "the source distribution or any confining volumes" << G4endl;
        G4cout << "may not exist"<< G4endl;
        G4cout << "If you have set confine then this will be ignored" <<G4endl;
        G4cout << "for this event." << G4endl;
        G4cout << "*************************************" << G4endl;
        srcconf = true; //Avoids an infinite loop
      }
  }

  G4PrimaryVertex* vertex = 
    new G4PrimaryVertex(particle_position,particle_time);

  // Angular stuff
  if(AngDistType == "iso")
    GenerateIsotropicFlux();
  else if(AngDistType == "direction")
    SetParticleMomentumDirection(particle_momentum_direction);
  else
    G4cout << "Error: AngDistType has unusual value" << G4endl;
  // Energy stuff

  G4int NOG = 0;

 // list of the decays available ; when you add or modify one please update also the Nuclide_status information file
  if(EnergyDisType == "Mono")
    GenerateMonoEnergetic();
  else if (EnergyDisType == "Ac228")
    GenerateAc228Decay(vertex);  
  else if (EnergyDisType == "Ba133")
    {NOG = GenerateBa133Decay(vertex);} //  G4cout << "Nb of Ba133 Gammas created: " << NOG << G4endl;}
  else if (EnergyDisType == "Bi210")
    GenerateBi210Decay(vertex);
  else if (EnergyDisType == "Bi212")
    GenerateBi212Decay(vertex);
  else if (EnergyDisType == "Bi214")
    GenerateBi214Decay(vertex);
  else if (EnergyDisType == "C14")
    GenerateC14Decay(vertex);
  else if (EnergyDisType == "Cf252")
    GenerateCf252Decay(vertex);
  else if (EnergyDisType == "Co57")
    GenerateCo57Decay(vertex);
  else if (EnergyDisType == "Co60")
    GenerateCo60Decay(vertex);
  else if (EnergyDisType == "Cs137")
    GenerateCs137Decay(vertex);
  else if (EnergyDisType == "Ga68")
    GenerateGa68Decay(vertex);
  else if (EnergyDisType == "K40")
    GenerateK40Decay(vertex);
  else if (EnergyDisType == "Pb210")
    GeneratePb210Decay(vertex);
  else if (EnergyDisType == "Pb212")
    GeneratePb212Decay(vertex);
  else if (EnergyDisType == "Pb214")
    GeneratePb214Decay(vertex);
  else if (EnergyDisType == "Po210")
    GeneratePo210Decay(vertex);
  else if (EnergyDisType == "Po210RN")
    GeneratePo210RNDecay(vertex);
  else if (EnergyDisType == "Po212")
    GeneratePo212Decay(vertex);
  else if (EnergyDisType == "Ra224")
    GenerateRa224Decay(vertex);
  else if (EnergyDisType == "Ra226")
    GenerateRa226Decay(vertex);  
  else if (EnergyDisType == "Sb125")
    GenerateSb125Decay(vertex);  
  else if (EnergyDisType == "Th228")
    GenerateTh228Decay(vertex);
  else if (EnergyDisType == "Th230")
    GenerateTh230Decay(vertex); 
  else if (EnergyDisType == "Th234")
    GenerateTh234Decay(vertex);
  else if (EnergyDisType == "Tl204")
    GenerateTl204Decay(vertex);
  else if (EnergyDisType == "U234")
    GenerateU234Decay(vertex); 
  else if (EnergyDisType == "PbBiPo")
    {G4double rand;   
     rand = G4UniformRand();
     if (rand < 0.3333333) 
	 {GeneratePb210Decay(vertex);}
     else if (rand < 0.6666667) 
	 {GenerateBi210Decay(vertex);}
     else 
	 {GeneratePo210Decay(vertex);}
     }
else
    G4cout << "Error: EnergyDisType has unusual value" << G4endl;
  
  // create a new vertex

  if(verboseLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;
  // create new primaries and set them to the vertex
  
  if(verboseLevel >= 1){
    G4cout << "Particle name: " 
	   << particle_definition->GetParticleName() << G4endl; 
    G4cout << "       Energy: "<<particle_energy << G4endl;
    G4cout << "     Position: "<<particle_position<< G4endl; 
    G4cout << "    Direction: "<<particle_momentum_direction << G4endl;
    G4cout << " NumberOfParticlesToBeGenerated: "
	   << NumberOfParticlesToBeGenerated << G4endl;
  }

  if (EnergyDisType == "Mono"){
    for( G4int i=0; i<NumberOfParticlesToBeGenerated; i++ ) {
      CreatePrimaryParticle(vertex);
    }
  }

  evt->AddPrimaryVertex( vertex );
  if(verboseLevel > 1)
    G4cout << " Primary Vetex generated "<< G4endl;   
}
