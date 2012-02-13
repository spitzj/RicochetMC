////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCNeutronWall.hh                                    //
//  Description: Example G4PrimaryGeneratorAction for RMC Simulation //
//                                                                    //
//  Author:      Dennis Wright                                        //
//  Date:        16 July 2010                                         //
//                                                                    //
//  20101026  M. Kelsey -- Use RMCVSourceConstruction as base, add   //
//            data members for all gun parameters (settable?)         //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
////////////////////////////////////////////////////////////////////////

#ifndef RMCNeutronWall_hh
#define RMCNeutronWall_hh 1

#include "RMCg4base/RMCVSourceConstruction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;


class RMCNeutronWall: public RMCVSourceConstruction {
public:
  RMCNeutronWall();    
  virtual ~RMCNeutronWall();
  
  void GeneratePrimaries(G4Event*);

private:
  G4int particlesPerEvent;
  G4String particleName;
};

#endif	/* RMCNeutronWall_hh */
