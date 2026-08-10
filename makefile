# simple make file
SRCS=main.f90
OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

# Compiler/Linker settings
FC = gfortran
export LIBRARY_PATH += :$(HOME)/.rsg_libs
INCLUDE_PATH := $(HOME)/.rsg_libs/include
FCFLAGS =  -I $(INCLUDE_PATH) -O2 -c -cpp -DDEBUG -Wall -DDEBUG -DASSERTS -Wextra -Wconversion  -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds -Wno-tabs  -Wall #-fmax-errors=5
FLFLAGS =  -l2CO2_PES_Aug26
PROGRAM = 2CO2_PES.out
PRG_OBJ = $(PROGRAM).o


default:
	$(MAKE) -C src all
	$(MAKE) programme

programme :$(PROGRAM)

all:
	$(MAKE) -C src all
	$(MAKE)

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) -o $@ $^ $(FLFLAGS)

# If something doesn't work right, have a 'make debug' to
# show what each variable contains.
debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "MODS = $(MODS)"
	@echo "MOD_OBJS = $(MOD_OBJS)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRG_OBJ = $(PRG_OBJ)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(patsubst %.o,%.mod,$(MOD_OBJS)) *.mod
	$(RM) $(PROGRAM)
	$(MAKE) -C src clean
	$(MAKE) -C tests clean

.PHONY: debug default clean
