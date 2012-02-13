#ifndef RMCMaterialTable_hh
#define RMCMaterialTable_hh 1
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        RMCMaterialTable.hh                                  //     
//  Description: Singleton manager for all material construction      //
//                                                                    //
//  Author:      Adam Anderson (MIT)                                  //
//               (adapted from Mike Kelsey (SLAC))                    //
//  Date:        8 January 2012                                       //
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

class RMCMaterialTable {
public:
  static G4Material* GetMaterial(const G4String& mat);

private:
  // Unchangeable object will be created in .cc file, forcing material builds
  static const RMCMaterialTable Instance;

  RMCMaterialTable();		// All materials will be created here
  ~RMCMaterialTable() {}

  void CreateNorite();		// Add other custom materials here
  void CreateGreenstone();
  void CreateMuMetal();
  void CreateScintillator();
};

#endif	/* RMCMaterialTable_hh */
