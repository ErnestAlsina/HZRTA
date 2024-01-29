# Makefile for building the MCMC F90 code for meta

# Select compiler uncommenting the COMPILER and OPTIONS variables accordingly

# Intel Fortran compiler
#COMPILER = f95 
#OPTIONS = -c 
# DEBUG = -g -traceback -check noarg_temp_created

# GFortran
 COMPILER = gfortran
  OPTIONS = -c -frecord-marker=4 -O2 -funroll-loops

OBJECTS = hz_rta.o rt_formal.o matrmod.o redis.o abs_prof.o math.o profile.o parameters.o global_vars.o

metalevel: $(OBJECTS)
	$(COMPILER) $(OBJECTS) -o hz_rta
clean:
	find . -maxdepth 1 -name "*.o" -delete ; find . -maxdepth 1 -name "*.mod" -delete

global_vars.o: global_vars.f90
	$(COMPILER) $(OPTIONS) global_vars.f90

parameters.o: parameters.f90
	$(COMPILER) $(OPTIONS) parameters.f90

math.o: math.f90 parameters.o
	$(COMPILER) $(OPTIONS) math.f90

profile.o: profile.f90 parameters.o
	$(COMPILER) $(OPTIONS) profile.f90

abs_prof.o: abs_prof.f90 parameters.o global_vars.o profile.o math.o
	$(COMPILER) $(OPTIONS) abs_prof.f90

redis.o: redis.f90 parameters.o global_vars.o profile.o math.o abs_prof.o
	$(COMPILER) $(OPTIONS) redis.f90

rt_formal.o: rt_formal.f90 parameters.o global_vars.o profile.o math.o abs_prof.o redis.o
	$(COMPILER) $(OPTIONS) rt_formal.f90

matrmod.o: matrmod.f90 parameters.o global_vars.o profile.o math.o abs_prof.o redis.o rt_formal.o
	$(COMPILER) $(OPTIONS) matrmod.f90

hz_rta.o: hz_rta.f90 parameters.o global_vars.o profile.o abs_prof.o math.o redis.o matrmod.o rt_formal.o
	$(COMPILER) $(OPTIONS) hz_rta.f90
