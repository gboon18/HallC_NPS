# $Id: GNUmakefile,v 1.1 1999-01-07 16:05:40 gunter Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := DVCS
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
#  G4INSTALL = /pbs/software/centos-7-x86_64/geant4/10.03.p03/share/Geant4-10.3.3/geant4make
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

##################################### ROOT FLAGS and Libraries ##############################

ifdef ROOTSYS
ifndef G4UI_USE_ROOT
   ROOTCPPFLAGS   = $(shell $(ROOTSYS)/5.34.36/bin/root-config --cflags)
   Gen_Sources    = "-I/sps/clas/hosanko/NPS/DVCS/dvcs_gen" #"-I/pbs/throng/clas/hosanko/DVCS/dvcs_gen"

   Gen_Dir 	  = /sps/clas/hosanko/NPS/DVCS/dvcs_gen/ #/pbs/throng/clas/hosanko/DVCS/dvcs_gen/
   DVCS_DIR	  = ${NPS_SOFT}
#   ANALYZER_DIR	  = ${ANALYZER}

   DVCS_SRC	  = "-I/pbs/throng/clas/hosanko/NPS/NPS_SOFT" #"-I/pbs/throng/clas/alex/soft/soft-frederic"
#   ALALYZER_SRC	  = -I$(ANALYZER_DIR)/src
#   ANALYZER_SRC  += -I$(ANALYZER_DIR)/hana_scaler
#   ANALYZER_SRC  += -I$(ANALYZER_DIR)/hana_decode
#   INCPATH        = -I/usr/include/Qt #-I/usr/include/qt5 #-I/usr/include/Qt

   #ROOTCPPFLAGS  := $(filter-out -D_REENTRANT,$(ROOTCPPFLAGS))
   #ROOTCPPFLAGS  := $(filter-out -pthread,$(ROOTCPPFLAGS))
   CPPFLAGS      += -DG4ANALYSIS_USE_ROOT $(ROOTCPPFLAGS)
   CPPFLAGS      += $(Gen_Sources)
   CPPFLAGS	 += $(DVCS_SRC)
#   CPPFLAGS      += $(ANALYZER_SRC)
#   CPPFLAGS      += $(INCPATH)

   ROOTLIBS       = $(shell $(ROOTSYS)/5.34.36/bin/root-config --nonew --glibs)
   ROOTLIBS      := $(filter-out -lNew,$(ROOTLIBS))
   TGenBaseSO     = -lTGenBase
   TGenGeoSO      = -lTGenGeo 
   TGenDVCSSO     = -lTGenDVCS
   TGVKellySO	  = -lTGVKelly
   SO_LIBS =  -L$(Gen_Dir) $(TGenBaseSO) $(TGenGeoSO) $(TGenDVCSSO) $(TGVKellySO)
   DVCS_Libs = -L$(DVCS_DIR) -lDVCS
#   ANALYZER_LIBS = -L$(ANALYZER_DIR) -lHallA
#   ANALYZER_LIBS += -L$(ANALYZER_DIR)/hana_scaler -lscaler
#   ANALYZER_LIBS += -L$(ANALYZER_DIR)/hana_decode -ldc

   #ROOTLIBS      := $(filter-out -lThread,$(ROOTLIBS))
   #ROOTLIBS      := $(filter-out -lpthread,$(ROOTLIBS))
   #ROOTLIBS      := $(filter-out -pthread,$(ROOTLIBS))
#   ROOTLIBS      += -lMinuit
 #  ROOTLIBS      += -lHallA
 #  ROOTLIBS      += -lDVCS
# LDLIBS        += $(ROOTLIBS)
EXTRALIBS	+= $(ROOTLIBS)
EXTRALIBS	+= $(SO_LIBS)
EXTRALIBS	+= $(DVCS_Libs)
#EXTRALIBS	+= $(ANALYZER_LIBS)
EXTRALIBS	+= /pbs/throng/clas/carlos/dvcs/soft/hrse_acceptance_maps2.o

 endif
 endif


#ROOTCFLAGS	= $(shell$(ROOTSYS)/bin/root-config --cflags)
#CPPFLAGS	+= $(ROOTCFLAGS)
#	echo $(CPPFLAGS)
#ROOTLIBS	= $(shell $(ROOTSYS)/bin/root-config  --glibs)
#EXTRALIBS	+= $(ROOTLIBS)
include $(G4INSTALL)/config/binmake.gmk


#include $(G4INSTALL)/config/binmake.gmk
