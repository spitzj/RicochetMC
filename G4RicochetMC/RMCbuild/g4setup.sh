#!/bin/sh
#  $Id: g4setup.sh,v 1.8 2011/11/13 14:33:28 kelsey Exp $
###################################################################
#                                                                 #
#  g4setup.sh - environment setup file for RMC GEANT4 use        #
#                                                                 #
#  Author: Michael Kelsey (SLAC)                                  #
#  Date:   29 July 2010                                           #
#                                                                 #
# 20100729  M. Kelsey -- Rewrite to handle site configurations    #
# 20100830  M. Kelsey -- Check for user override before setup,    #
#		fix bug in G4ENV pre-check.                       #
# 20110121  M. Kelsey -- Protect against users' "cd" aliases      #
# 20110606  M. Kelsey -- Remove debugging, protect |unalias| use  #
# 20111112  M. Kelsey -- Remove old ROOTSYS from user's PATH      #
###################################################################

# Sanity check -- this file must be |source|d, not executed

if [ -n "`echo $0|grep g4setup.sh`" ]; then
   echo "The $0 script must be run using the . command."
   exit 2
fi

# Useful magic for changing PATH values (works with BASH only)

delpath()
{
  eval "$1=\$(echo \$$1 | sed -e s%^$2\$%% -e s%^$2\:%% -e s%:$2\:%:%g -e s%:$2\\\$%%)"
}

addpath()
{
  if eval test -z \"\${$1}\" ; then
    eval "$1=$2"
  elif ! eval test -z \"\${$1##$2}\" -o -z \"\${$1##*:$2:*}\" -o -z \"\${$1%%*:$2}\" -o -z \"\${$1##$2:*}\" ; then
    eval "$1=\$$1:$2"
  fi
}

# Name of library path is different on MacOSX

LIBPATH=LD_LIBRARY_PATH
[ `uname` = Darwin ] && LIBPATH=DYLD_LIBRARY_PATH

# If user has preset ROOTSYS, remove it from paths

if [ -n "$ROOTSYS" ]; then
   eval test -n \"\${$LIBPATH}\" && eval delpath $LIBPATH $ROOTSYS/lib
   delpath PATH $ROOTSYS/bin
fi

# Get site configuration settings if user hasn't preset things

[ -z "$G4INSTALL" ] && eval `RMCbuild/RMCsite.sh`

# If no base settings, warn user and quit

if [ -z "$G4INSTALL" -o -z "$ROOTSYS" ]; then
   echo "No local GEANT4 or ROOT installation.  RMC build will not work."
   exit 1
fi

# GEANT4 configuration, including local G4WORKDIR

if [ -n "$G4WORKDIR" ]; then
   [ `type -t cd` = "alias" ] && unalias cd	### Don't use user's aliases
   wkdir=`cd $G4WORKDIR; pwd`
   [ $wkdir != `pwd` ] && echo "$0 : G4WORKDIR must point to RMC area."
else
   export G4WORKDIR=`pwd`
fi

[ -z "$G4ENV" ] && export G4ENV=env
. $G4INSTALL/${G4ENV}.sh

# Add ROOT libraries to search paths

eval addpath $LIBPATH $ROOTSYS/lib
addpath PATH $ROOTSYS/bin

# FIXME:  Need to get rid of these user-changing lines!

# Uncomment this to do optimized builds
#
# export G4OPTIMIZE=1

# Decide which vizualization driver to use and un-comment
# the appropriate lines
#
export G4VIS_USE_OPENGLX=1
#  export G4VIS_USE_VRML=1
#  export G4VIS_USE_DAWN=1
#  export G4VIS_USE_DAWNFILE=1
#  export G4VIS_USE_VRMLFILE=1
#  export G4VIS_USE_RAYTRACER=1
#  export G4VIS_USE_ASCIITREE=1
#  export G4VIS_USE_GAGTREE=1
#  export G4DAWNFILE_VIEWER=david


