FC = mpif90-openmpi-mp
FFLAGS = -O3 -fdefault-real-8 -fdefault-double-8 
CODEHOME = ~/Dropbox/papers/UltimaThule/fortran

main: $(CODEHOME)/yoshikozai.f90
	$(FC) $(FFLAGS) $(CODEHOME)/yoshikozai.f90 -o yoshikozai.x

clean:
	rm -f yoshikozai.x