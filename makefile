# simple make file
SRCS=2CO2_PES_GP_Symm.f90 testing_2CO2.f90  PES_GP_Symm.f90 


OBJS=$(patsubst %.f90,%.o,$(SRCS))

# Ditto for mods (They will be in both lists)
MODS=$(wildcard mod*.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

# Compiler/Linker settings
FC = gfortran
FLFLAGS = 
FCFLAGS =  -c -cpp -DDEBUG -Wall -DDEBUG -Wextra -Wconversion -fcheck=bounds -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow  #-fmax-errors=5
PROGRAM = 2CO2_PES.out
PRG_OBJ = $(PROGRAM).o

# make without parameters will make first target found.
default : $(PROGRAM)

# Compiler steps for all objects
$(OBJS) : %.o : %.f90
	$(FC) $(FCFLAGS) -o $@ $<

# Linker
$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^


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

.PHONY: debug default clean

# Dependencies

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)

# Blocks and allocations depends on shared
mod_blocks.o mod_allocations.o : mod_shared.o
