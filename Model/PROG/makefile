NETCDFHOME="/opt/cray/pe/netcdf/4.8.1.3/GNU/8.2"
NETCDFINC="-I$(NETCDFHOME)/include"
NETCDFLIBS="-L$(NETCDFHOME)/lib -lnetcdff -lnetcdf"

all: MPM

MPM: MPM_Utils.F90  MPM_Procedures.F90  MPM_SSA.F90
	elmerf90 $(NETCDFINC) $^ -o $@ $(NETCDFLIBS)

clean:
	rm -f MPM
	rm -f *.o
	rm -f *.mod
