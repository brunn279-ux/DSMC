### Simple Makefile for MGDS_lite post-processor ###

##### Compiler options
FC    = gfortran
FOPTS = -fPIC -O3 -ffixed-line-length-132
#FOPTS = -fPIC -g -O0 -fbounds-check -fbacktrace -ffixed-line-length-132

# Need modules from main code (and from precompiled viz_struct/io_viz now in ../src)
INCLUDE = -I ../src -I .

############################ Targets ##########################
all: check_src code link

# Extremely compact check that the precompiled objects exist
check_src:
	@test -f ../src/viz_struct.o -a -f ../src/io_viz.o -a -f ../src/main_struct.o || \
	( echo "ERROR: Missing precompiled objects in ../src (need viz_struct.o, io_viz.o, main_struct.o)"; exit 1 )

code:
	$(FC) -c $(FOPTS) $(INCLUDE) -J . amr.f90
	$(FC) -c $(FOPTS) $(INCLUDE) -J . surface_out.f90
	$(FC) -c $(FOPTS) $(INCLUDE) -J . flow_out.f90
	$(FC) -c $(FOPTS) $(INCLUDE) -J . main_viz.f90

link:
	$(FC) amr.o surface_out.o flow_out.o main_viz.o \
	      ../src/viz_struct.o ../src/io_viz.o ../src/main_struct.o \
	      -o viz_mgds_lite

clean:
	rm -f *.o *.mod viz_mgds_lite
