#!/bin/sh
#
# Usage:  eval `CDMSbuild/CDMSsite.sh`
#
# Selects and executes site-specific CDMSbuild/config.<site> script to
# define local GEANT4 and ROOT installation variables.  Selection is done
# as follows:
#
# 1) If CDMS_SITE is defined, use it.  Otherwise,
# 2) Get DNS domain.  If config.* exists, use it.  Otherwise,
# 3) Get local host.  If config.* exists, use it.
#
# For cases (2) and (3), CDMS_SITE will be set accordingly.  
#
# $Id: CDMSsite.sh,v 1.4 2010/08/31 16:48:57 kelsey Exp $
# 20100729  Michael Kelsey
# 20100825  Use |dnsdomainname|, not NIS; define alias if command not found
# 20100830  Need to quote back-quoted strings in test
# 20100831  Allow for second "=" in value string

# Get path to this script to find matching config.* files

cdmspath=`dirname $0`

# If no site name has been defined, figure it out

if [ -z "`which dnsdomainname 2>/dev/null`" ]; then
   dnsdomainname() {
     nslookup `hostname -s`|awk '/Name:/{print $2}'|cut -f2- -d.
   }
fi

if [ -z "$CDMS_SITE" ]; then
   try_site=`dnsdomainname`
   if [ -z "$try_site" -o ! -r "$cdmspath/config.$try_site" ]; then
      try_site=`hostname -s`
      if [ -z "$try_site" -o ! -r "$cdmspath/config.$try_site" ]; then
         try_site=""
      fi
   fi
   [ -n "$try_site" ] && export CDMS_SITE=$try_site
fi

# If no site name can be defined, just quit

[ -z "$CDMS_SITE" ] && exit

# Now just echo out the minimal assignment statements (skipping comments)

sitefile=$cdmspath/config.$CDMS_SITE

echo "export CDMS_SITE=$CDMS_SITE ;"
sed -n 's/\(^[^#][^=]*=[^#]*\).*/export \1 ;/p' $sitefile
