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
include ../makefile.header


#--- compile depedency ---
all:	TTDriver.exe
	cp TTDriver.exe $(bin)

TTDriver.exe	:	main.o
			$(FC) $< $(FLAG) -L$(TTsrc) -lTTsrc -o TTDriver.exe

#--- modules' dependency ---
main.o 		:	main.f90
			$(FC) -c $(FLAG) -I$(TTsrc) -I$(NETCDFINC) $<

clean:
	rm *.o *.exe

#EOC
