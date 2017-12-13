include ../../make.inc

LIBXC_SRC = \
libxc_funcs.f90 libxc.f90 

OBJ = $(LIBXC_SRC:.f90=.o) 

#
# Suffix rule for Fortran 90
#
%.o : %.f90
	$(F90) $(F90_OPTS) -c -o $(*F).o $<


# Targets

# we only need the *mod files
lib: $(OBJ)
	cp *.o *.mod ../../

# does not delete *.x files
clean:
	rm -rf *.o *.mod 




