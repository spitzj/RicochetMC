#  $Id: GNUmakefile,v 1.9 2010/12/23 16:51:56 kelsey Exp $
#
#  Top-level makefile for building CDMS general geometry simulation
#  ----------------------------------------------------------------

# Library is build from sources in several places
SUBDIRS := labs detectors interface

# Binary is just a unit test
BINS := CDMS_GGSim

# Executable depends on ROOT through CDMSactions/CDMSSimRoot class
ROOTLIBS := $(shell $(ROOTSYS)/bin/root-config --glibs) -lMinuit -lHtml
ROOTLIBS := $(filter-out -lNew,$(ROOTLIBS))
ROOTLIBS := $(filter-out -lpthread,$(ROOTLIBS))

CDMSLIBS := -lCDMSactions -lCDMSgeometry -lCDMSsources -lCDMSg4base

LDLIBS += $(ROOTLIBS)

include build.gmk
