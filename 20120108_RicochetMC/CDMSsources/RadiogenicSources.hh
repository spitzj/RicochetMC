#ifndef RadiogenicSources_hh
#define RadiogenicSources_hh
// $Id: RadiogenicSources.hh,v 1.3 2011/01/05 19:16:21 kelsey Exp $
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

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "G4ThreeVector.hh"
#include <globals.hh>

class G4Event;
class G4Navigator;
class G4ParticleDefinition;
class G4ParticleTable;
class G4PrimaryVertex;
class RadiogenicSourcesMessenger;


class RadiogenicSources : public CDMSVSourceConstruction {
public:
  RadiogenicSources();
  virtual ~RadiogenicSources();

  virtual void GeneratePrimaries(G4Event *evt);

  // spatial distribution  
  void SetPosDisType(const G4String& PosType);
  void SetPosDisShape(const G4String& shapeType) { Shape = shapeType; }
  void SetCentreCoords(const G4ThreeVector& coordsOfCentre) {
    CentreCoords = coordsOfCentre;
  }
  void SetCentreR(G4double rofcentre) { centre_r = rofcentre; }
  void SetCentrePhi(G4double phiofcentre) { centre_phi = phiofcentre; }
  void SetCentreZ(G4double zofcentre) { centre_z = zofcentre; }
  void SetHalfZ(G4double zhalf) { halfz = zhalf; }
  void SetRadius(G4double rad) { Radius = rad; }
  void SetCEflag(G4int);
  void GeneratePointSource();
  void GeneratePointsInVolume();
  void SetRandomPosInMaterial();
  void CreatePrimaryParticle(G4PrimaryVertex *vt);
  G4bool IsSourceConfined();
  void ConfineSourceToVolume(const G4String&);
  
  // angular distribution
  void SetAngDistType(const G4String&);
  void SetParticleMomentumDirection(const G4ThreeVector& aDirection) {
    particle_momentum_direction =  aDirection.unit();
  }
  void GenerateIsotropicFlux();

  // energy distribution 
  void SetEnergyDisType(const G4String&);
  void ReadEnergyProbaFile(const G4String&, G4int);
  void ReadIntProbaFile(const G4String&, G4int);
  double BetaEnergyFromSpectrum(G4int, G4int);
  void GenerateAc228Decay(G4PrimaryVertex *vert);
  int  GenerateBa133Decay(G4PrimaryVertex *vert);
  void GenerateBi210Decay(G4PrimaryVertex *vert);
  void GenerateBi212Decay(G4PrimaryVertex *vert);
  void GenerateBi214Decay(G4PrimaryVertex *vert);
  void GenerateC14Decay(G4PrimaryVertex *vert);
  void GenerateCf252Decay(G4PrimaryVertex* vert);
  void GenerateCo57Decay(G4PrimaryVertex *vert);
  void GenerateCo60Decay(G4PrimaryVertex *vert);
  void GenerateCs137Decay(G4PrimaryVertex *vert);
  void GenerateGa68Decay(G4PrimaryVertex *vert);
  void GenerateK40Decay(G4PrimaryVertex *vert);
  void GeneratePb210Decay(G4PrimaryVertex *vert);
  void GeneratePb212Decay(G4PrimaryVertex *vert);
  void GeneratePb214Decay(G4PrimaryVertex *vert);
  void GeneratePo210Decay(G4PrimaryVertex *vert);
  void GeneratePo210RNDecay(G4PrimaryVertex *vert);
  void GeneratePo212Decay(G4PrimaryVertex *vert);
  void GenerateRa224Decay(G4PrimaryVertex *vert);
  void GenerateRa226Decay(G4PrimaryVertex *vert);
  void GenerateSb125Decay(G4PrimaryVertex *vert);
  void GenerateTh228Decay(G4PrimaryVertex *vert);
  void GenerateTh230Decay(G4PrimaryVertex *vert);
  void GenerateTh234Decay(G4PrimaryVertex *vert);
  void GenerateTl204Decay(G4PrimaryVertex *vert);
  void GenerateU234Decay(G4PrimaryVertex *vert);

  void SetMonoEnergy(G4double menergy) { MonoEnergy = menergy; }
  void GenerateMonoEnergetic() { particle_energy = MonoEnergy; }
  inline G4double GetParticleEnergy() { return particle_energy; }

  // particle properties
  void SetParticleDefinition(G4ParticleDefinition * aParticleDefinition);
  inline void SetParticleCharge(G4double aCharge)
  { particle_charge = aCharge; }
  
private:
  // position distribution
  G4String SourcePosType;
  G4String Shape;
  G4double centre_z;
  G4double centre_r;
  G4double centre_phi;
  G4double halfz;
  G4double Radius;
  G4ThreeVector CentreCoords;
  G4bool Confine;
  G4String VolName;
  G4String AngDistType;
  G4double MinTheta, MaxTheta, MinPhi, MaxPhi;
  G4double Theta, Phi;
  G4String EnergyDisType;
  G4int ceFlag;
  G4double MonoEnergy;
  
  // particle properties 
  G4int                  NumberOfParticlesToBeGenerated;
  G4ParticleDefinition*  particle_definition;
  G4ThreeVector          particle_momentum_direction;
  G4double               particle_energy;
  G4double               particle_charge;
  G4ThreeVector          particle_position;
  G4double               particle_time;
  G4ThreeVector          particle_polarization;
  
  G4Navigator *gNavigator;
  G4String myMaterial;
  G4int bot,zyl,top,ass;
  
  G4int count; G4int nlfile[10];
  G4double beta_eng[2000];
  G4double beta_prob[2000];
  
  G4int Bi210rank;
  G4int C14rank;
  G4int Cf252rank;
  G4int Cs137rank;
  G4int Ga68rank;
  G4int K40rank1, K40rank2;
  G4int Pb210rank1, Pb210rank2;
  G4int Sb125rank;
  G4int Tl204rank;

  G4ParticleTable* particleTable;
  RadiogenicSourcesMessenger *theMessenger;
};

#endif	/* RadiogenicSources_hh */
