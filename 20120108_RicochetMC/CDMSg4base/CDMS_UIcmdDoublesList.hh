#ifndef CDMS_UIcmdDoublesList_hh
#define CDMS_UIcmdDoublesList_hh 1
// $Id: CDMS_UIcmdDoublesList.hh,v 1.1 2010/12/23 16:48:42 kelsey Exp $
////////////////////////////////////////////////////////////////////////
//                                                                    //
//  File:        CDMS_UIcmdDoublesList.hh                             //     
//  Description: Instance of templated class for std::vector<double>  //
//                                                                    //
//  Author:      Michael Kelsey (SLAC)                                //
//  Date:        22 December 2010                                     //
//                                                                    //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSg4base/CDMS_TUIcmdWithList.hh"
#include "globals.hh"

typedef CDMS_TUIcmdWithList<G4double> CDMS_UIcmdDoublesList;

#endif	/* CDMS_UIcmdDoublesList_hh */
