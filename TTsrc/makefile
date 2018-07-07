#------------------------------------------------------------------------------
#             Tan-Tracker Carbon cycle Data Assimilation System               #
#------------------------------------------------------------------------------
#BOP
#
#	!Description : Tan-Tracker Source Code TTsrc dir
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
include ../makefile.header

#---  List of source files: everything ending in .F ---
SOURCES		:=$(wildcard *.f90) 
OBJECTS		:=$(SOURCES:.f90=.o)
MODULES		:=$(OBJECTS:.o=.mod)

#--- Library file ---
LIBRARY   	:=libTTsrc.a

#--- compile depedency ---
all: 		$(OBJECTS) $(LIBRARY)

$(LIBRARY)	:	$(OBJECTS)
				ar rcs $(LIBRARY) $(OBJECTS)

#--- modules' dependency ---
#	level 1
inputoption.o	:		inputoption.f90
						$(FC) -c $< $(FLAG) 

common.o		:		common.f90
						$(FC) -c $< $(FLAG) 

mathmatic.o 	:		mathmatic.f90
						$(FC) -c $< $(FLAG)

statistic.o 	:		statistic.f90
						$(FC) -c $< $(FLAG)

reshapeTT.o 	:		reshapeTT.f90
						$(FC) -c $< $(FLAG)

#	level 2
readwrite_nc.o	:		readwrite_nc.f90	\
						common.o
						$(FC) -c $< $(FLAG) -I$(NETCDFINC)

extract_bpch.o	:		extract_bpch.f90	\
						inputoption.o 		\
						common.o
						$(FC) -c $< $(FLAG) 

random.o 		:		random.f90 			\
						mathmatic.o
						$(FC) -c $< $(FLAG)

interpolation.o :		interpolation.f90	\
						common.o
						$(FC) -c $< $(FLAG)

timecounter.o 	:		timecounter.f90 	\
						common.o 			\
						inputoption.o
						$(FC) -c $< $(FLAG)

podcore.o 		:		podcore.f90 		\
						common.o 			\
						inputoption.o 		\
						mathmatic.o 		\
						reshapeTT.o
						$(FC) -c $< $(FLAG)

#	level 3
dainout.o 	:			dainout.f90 		\
						common.o 			\
						inputoption.o 		\
						statistic.o 		\
						timecounter.o 		\
						readwrite_nc.o 		\
						extract_bpch.o 		\
						reshapeTT.o
						$(FC) -c $< $(FLAG)

operator_ocogst.o	:	operator_ocogst.f90  \
						common.o 			\
						extract_bpch.o 		\
						interpolation.o 	\
						statistic.o 		\
						inputoption.o
						$(FC) -c $< $(FLAG) -I$(NETCDFINC)

#	level 4
operator_tansat.o	:	operator_tansat.f90  \
						common.o 			\
						operator_ocogst.o 	\
						extract_bpch.o 		\
						interpolation.o 	\
						inputoption.o
						$(FC) -c $< $(FLAG) -I$(NETCDFINC)

fast_loc.o 		:		fast_loc.f90 		\
						common.o 			\
						inputoption.o 		\
						readwrite_nc.o 		\
						operator_ocogst.o 	\
						mathmatic.o 		\
						reshapeTT.o
						$(FC) -c $< $(FLAG)

#	level 5
nlscore.o 		:		nlscore.f90 		\
						common.o 			\
						inputoption.o 		\
						mathmatic.o 		\
						fast_loc.o 			\
						reshapeTT.o
						$(FC) -c $< $(FLAG)

clean:		
	rm -f *.mod *.o *.a

#EOC
