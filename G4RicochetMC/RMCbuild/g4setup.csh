#!/bin/csh
#  $Id: g4setup.csh,v 1.9 2011/11/13 14:33:28 kelsey Exp $
###################################################################
#                                                                 #
#  CDMSg4setup.csh - environment setup file for CDMS GEANT4 use   #
#                                                                 #
#  Author: Dennis Wright (SLAC)                                   #
#  Date:   21 July 2010                                           #
#                                                                 #
# 20100729  M. Kelsey -- Rewrite to handle site configurations    #
# 20100830  M. Kelsey -- Check for user override before setup     #
# 20110121  M. Kelsey -- Protect against users' "cd" aliases      #
# 20111112  M. Kelsey -- Remove old ROOTSYS from user's PATH      #
###################################################################

# Sanity check -- this file must be |source|d, not executed

if ("$0" =~ *g4setup.csh) then
   echo "The $0 script must be run using the |source| command."
   exit 2
endif

# Useful magic for changing PATH values

alias delpath 'setenv \!:1 `echo ${\!:1} | sed -e s%^\!:2\$%% -e s%^\!:2\:%% -e s%:\!:2\:%:%g -e s%:\!:2\$%%`'
alias addpath 'if ( $\!:1 != \!:2 && $\!:1 !~ \!:2\:* && $\!:1 !~ *\:\!:2\:* && $\!:1 !~ *\:\!:2 ) setenv \!:1 ${\!:1}\:\!:2'

# Name of library path is different on MacOSX

set LIBPATH = LD_LIBRARY_PATH
if (`uname` == Darwin) set LIBPATH = DYLD_LIBRARY_PATH

# If user has preset ROOTSYS, remove it from paths

if ($?ROOTSYS) then
   if (`env|grep "^${LIBPATH}="` != "") then
      eval delpath $LIBPATH $ROOTSYS/lib
   endif
   delpath PATH $ROOTSYS/bin
endif

# Get site configuration settings if user hasn't preset things

if (! $?G4INSTALL) eval `RMCbuild/CDMSsite.csh`

# If no base settings, warn user and quit

if (! $?G4INSTALL || ! $?ROOTSYS) then
   echo "No local GEANT4 or ROOT installation.  CDMS build will not work."
   exit 1
endif

# GEANT4 configuration, including local G4WORKDIR

if ($?G4WORKDIR) then
   unalias cd				### Don't use user's aliases
   set wkdir = `cd $G4WORKDIR; pwd`
   if ("$wkdir" != `pwd`) echo "$0 : G4WORKDIR must point to CDMS area."
else
   setenv G4WORKDIR `pwd`
endif

if (! $?G4ENV) setenv G4ENV env
source $G4INSTALL/${G4ENV}.csh

# Add ROOT libraries to search paths

eval addpath $LIBPATH $ROOTSYS/lib
addpath PATH $ROOTSYS/bin

# FIXME:  Need to get rid of these user-changing lines!

# Uncomment this to do optimized builds
#
# setenv G4OPTIMIZE 1

# Decide which vizualization driver to use and un-comment
# the appropriate lines
#
setenv G4VIS_USE_OPENGLX 1
#  setenv G4VIS_USE_VRML 1
#  setenv G4VIS_USE_DAWN 1
#  setenv G4VIS_USE_DAWNFILE 1
#  setenv G4VIS_USE_VRMLFILE 1
#  setenv G4VIS_USE_RAYTRACER 1
#  setenv G4VIS_USE_ASCIITREE 1
#  setenv G4VIS_USE_GAGTREE 1
#  setenv G4DAWNFILE_VIEWER david


