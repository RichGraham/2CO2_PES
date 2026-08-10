FC := gfortran
SRCS := 2CO2_PES_GP_Symm.f90 testing_2CO2.f90
OBJS := $(SRCS:%.f90=%.o)
DEPS := $(OBJS:.o=.d)

MODS=$(wildcard mod*.f90)
MOD_OBJS=$(patsubst %.f90,%.o,$(MODS))

INCLUDE_PATH := $(HOME)/.rsg_libs/include
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))

FCFLAGS =  -I $(INCLUDE_PATH) -fPIC -O2 -c -cpp -DASSERTS -Wall -pedantic-errors -Wsurprising -std=f2018 -Wextra -Wconversion  -ffpe-trap=invalid -ffpe-trap=zero,overflow,underflow -fbacktrace -fdump-core -fcheck=bounds  -Wno-tabs -fall-intrinsics #-Werror #-fmax-errors=5

FLFLAGS =

all: lib2CO2_PES_Aug26.a

lib2CO2_PES_Aug26.a: $(OBJS)
	$(AR) $(FLFLAGS) -r $@ $?
	ln -s -f $(ROOT_DIR)$@ $(HOME)/.rsg_libs
	ln -s -f $(ROOT_DIR)*.mod $(HOME)/.rsg_libs/include/

%.mod %.o : %.f90
	$(FC) -c -MMD -MP $(FCFLAGS)  $< -o $@

%.o : %.f90
	$(FC) $(FCFLAGS)  -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a *.d

