# $Id: GNUmakefile.rel,v 1.2 2010/08/31 00:07:39 kelsey Exp $
#
# RMC top-level GNUmakefile -- controls building of package libraries
# and executables, and self-configures working/release area.
#
# 20100722  Michael Kelsey

ifdef VERBOSE
  MAKEVBS := -w
else
  MAKEVBS := -s
endif

export NOTFIRST
export VERBOSE
export MAKEVBS

# Do a recursive MAKE up front to set search path for includes
ifndef NOTFIRST
.PHONY : $(MAKECMDGOALS)
$(MAKECMDGOALS) :
	@echo `date` : RMC GNUmakefile starting ...
	@$(MAKE) $(MAKEVBS) -I RMCbuild NOTFIRST=yes $@
else	# NOTFIRST
# See targets.gmk for top-level targets, and packages.gmk for within-packages
include arch.gmk
include commands.gmk
include packages.gmk
include targets.gmk
include help.gmk

endif	# NOTFIRST
