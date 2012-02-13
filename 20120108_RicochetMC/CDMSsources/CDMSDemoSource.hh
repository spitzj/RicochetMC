// $Id: CDMSDemoSource.hh,v 1.5 2010/12/23 16:49:48 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSDemoSource.hh                                    //
//  Description: Example G4PrimaryGeneratorAction for CDMS Simulation //
//                                                                    //
//  Author:      Dennis Wright                                        //
//  Date:        16 July 2010                                         //
//                                                                    //
//  20101026  M. Kelsey -- Use CDMSVSourceConstruction as base, add   //
//            data members for all gun parameters (settable?)         //
//  20101210  M. Kelsey -- Follow geometry class change to dimensions //
////////////////////////////////////////////////////////////////////////

#ifndef CDMSDemoSource_hh
#define CDMSDemoSource_hh 1

#include "CDMSg4base/CDMSVSourceConstruction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;


class CDMSDemoSource: public CDMSVSourceConstruction {
public:
  CDMSDemoSource();    
  virtual ~CDMSDemoSource();
  
  void GeneratePrimaries(G4Event*);

private:
  G4int particlesPerEvent;
  G4String particleName;
};

#endif	/* CDMSDemoSource_hh */
