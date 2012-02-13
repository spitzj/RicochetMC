// $Id: CDMSNeutronWall.hh,v 1.5 2010/12/23 16:49:48 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSNeutronWall.hh                                    //
//  Description: Example G4PrimaryGeneratorAction for CDMS Simulation //
//                                                                    //
//  Author:      Dennis Wright                                        //
//  Date:        16 July 2010                                         //
//                                                                    //
//  20101026  M. Kelsey -- Use CDMSVSourceConstruction as base, add   //
//            data members for all gun parameters (settable?)         //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSNeutronWall_hh
#define CDMSNeutronWall_hh 1

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;


class CDMSNeutronWall: public CDMSVSourceConstruction {
public:
  CDMSNeutronWall();    
  virtual ~CDMSNeutronWall();
  
  void GeneratePrimaries(G4Event*);

private:
  G4int particlesPerEvent;
  G4String particleName;
};

#endif	/* CDMSNeutronWall_hh */
