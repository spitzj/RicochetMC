#ifndef CDMSMaterialTable_hh
#define CDMSMaterialTable_hh 1
// $Id: CDMSMaterialTable.hh,v 1.2 2011/07/28 05:58:49 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMSMaterialTable.hh                                 //     
//  Description: Singleton manager for all CDMS material construction //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        22 July 2011                                         //
//                                                                    //
//  All materials used in CDMS (regardless of whether they are used   //
//  in any particlar job) must be instantiated at the beginning of    //
//  the job.  NeutronHP depends on the active material database to    //
//  initialize its cross-section tables, and will crash if materials  //
//  are added after that initialization is done (e.g., at geometry    //
//  construction time).                                               //
//////////////////////////////////////////////////////////////////////// 

class G4String;
class G4Material;

class CDMSMaterialTable {
public:
  static G4Material* GetMaterial(const G4String& mat);

private:
  // Unchangeable object will be created in .cc file, forcing material builds
  static const CDMSMaterialTable Instance;

  CDMSMaterialTable();		// All materials will be created here
  ~CDMSMaterialTable() {}

  void CreateNorite();		// Add other custom materials here
  void CreateGreenstone();
  void CreateMuMetal();
  void CreateScintillator();
};

#endif	/* CDMSMaterialTable_hh */
