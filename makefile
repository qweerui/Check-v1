#------------------------------------------------------------------------------
#             Tan-Tracker Carbon cycle Data Assimilation System               #
#------------------------------------------------------------------------------
#BOP
#
#   !Description : Tan-Tracker Source Code TTsrc dir
#
#   Record of revision:
#              Date:        Programmer:         Description of change:
#        ------------     --------------       ------------------------
#        2018.01.03             Hanrui            NV2.0 original code
#EOP
#------------------------------------------------------------------------------
#BOC
#
#--- include header file ---
include makefile.header

all:
	cd $(TTsrc); make
	cd $(TTMain); make
	cd $(TimeCTL); make
	cd $(Sampling); make
	cd $(TimeBack); make
clean:
	cd $(TTsrc); make clean
	cd $(TTMain); make clean
	cd $(TimeCTL); make clean
	cd $(Sampling); make clean
	cd $(TimeBack); make clean
	rm $(bin)/*

#EOC
