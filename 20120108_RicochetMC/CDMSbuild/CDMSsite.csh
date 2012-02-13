#!/bin/csh
#
# Usage:  eval `CDMSbuild/CDMSsite.csh`
#
# Selects and executes site-specific CDMSbuild/config.<site> script to
# define local GEANT4 and ROOT installation variables, in CSH format (setenv).
# See CDMSbuild/CDMSsite for actual process.
#
# $Id: CDMSsite.csh,v 1.2 2010/08/31 16:48:57 kelsey Exp $
# 20100729  Michael Kelsey
# 20100831  Allow for second "=" in value string

# CDMSsite puts out "export VAR=value"; use |sed| to reformat for CSH

set cdmspath = `dirname $0`
$cdmspath/CDMSsite.sh | sed 's/export \([^=]*\)=\(.*\) ;/setenv \1 \2 ;/'
