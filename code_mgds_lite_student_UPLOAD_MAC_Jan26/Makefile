### Student Makefile (assumes precompiled src/*.o and src/*.mod exist) ###

# Compiler (MUST BE GFORTRAN)
FC    = gfortran

# A few different compiler options (uncomment the one you want):
## Fast, optimized compilation
FOPTS = -fPIC -O3 -ffixed-line-length-132
## Slower, but has good debugging flags
#FOPTS = -fPIC -g -O0 -fbounds-check -fbacktrace -ffixed-line-length-132

SRC_DIR  = src
DSMC_DIR = DSMC

# DSMC needs to find .mod files from the precompiled code
INCLUDE = -I$(SRC_DIR)
DSMC_MODS   = -I$(DSMC_DIR) -J$(DSMC_DIR)

# Precompiled objects (Should already be built!)
MGDS_OBJS = \
  $(SRC_DIR)/main.o \
  $(SRC_DIR)/readin.o \
  $(SRC_DIR)/init.o \
  $(SRC_DIR)/geometry.o \
  $(SRC_DIR)/random.o \
  $(SRC_DIR)/io.o \
  $(SRC_DIR)/particle.o \
  $(SRC_DIR)/main_struct.o \
  $(SRC_DIR)/surface.o \
  $(SRC_DIR)/dsmc_iteration.o \
  $(SRC_DIR)/trisort.o \
  $(SRC_DIR)/NASA_cutcell.o

# Student-editable DSMC sources
DSMC_SRCS = \
  $(DSMC_DIR)/user_struct.f90 \
  $(DSMC_DIR)/set_params.f90 \
  $(DSMC_DIR)/BCs.f90 \
  $(DSMC_DIR)/collide.f90 \
  $(DSMC_DIR)/extra_utils.f90

DSMC_OBJS = $(DSMC_SRCS:.f90=.o)
DSMC_LIB  = $(DSMC_DIR)/libMGDSuser.a

.PHONY: all clean check_src

all: mgds_lite

# Prebuilt check for source files
check_src:
	@missing=""; for f in $(MGDS_OBJS); do [ -f $$f ] || missing="$$missing $$f"; done; \
	if [ -n "$$missing" ]; then \
	  echo "ERROR: Missing precompiled object(s) in $(SRC_DIR):$$missing"; \
	  echo "       Restore correct source files in $(SRC_DIR) (.o/.mod files)."; \
	  exit 1; \
	fi

# Build DSMC library from DSMC sources
$(DSMC_LIB): $(DSMC_OBJS)
	ar -r $@ $^

# Compile DSMC objects in-place (DSMC/*.o, DSMC/*.mod)
$(DSMC_DIR)/%.o: $(DSMC_DIR)/%.f90
		$(FC) -c $(FOPTS) $(INCLUDE) $(DSMC_MODS) $< -o $@

# Final link (single flags set)
mgds_lite: check_src $(DSMC_LIB)
	$(FC) $(FOPTS) $(MGDS_OBJS) $(DSMC_LIB) -o $@

clean:
	rm -f mgds_lite
	rm -f $(DSMC_DIR)/*.o $(DSMC_DIR)/*.mod $(DSMC_LIB)

