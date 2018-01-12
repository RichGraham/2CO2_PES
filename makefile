# simple make file
SOURCES=2CO2_PES_GP_Symm.f90 PES_GP_Symm.f90
PRODUCT=2CO2_PES.out


all: $(PRODUCT)

$(PRODUCT) : $(SOURCES)
	gfortran -o $(PRODUCT) $(SOURCES)
