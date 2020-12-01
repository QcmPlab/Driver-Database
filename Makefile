#TO BE CHANGED BY USER:
#$ DRIVER NAME without .f90 extension
#$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
#$ PLATFORM: supported platform are intel, gnu
#$ EXECUTABLE TARGET DIRECTORY (default if $HOME/.bin in the PATH)
EXE=ed_bhz_2d
FC=mpif90
PLAT=gnu
DIREXE=$(HOME)/.bin

define colorecho	
	@tput setaf 6
	@tput bold
	@echo $1
	@tput sgr0
endef


#NO NEED TO CHANGE DOWN HERE, only expert mode.
#########################################################################
GLOB_INC=$(shell pkg-config --cflags dmft_tools scifor)
GLOB_LIB=$(shell pkg-config --libs   dmft_tools scifor)

LANC_INC=$(shell pkg-config --cflags lanc_ed)
LANC_LIB=$(shell pkg-config --libs   lanc_ed)

DMFT_INC=$(shell pkg-config --cflags dmft_ed)
DMFT_LIB=$(shell pkg-config --libs   dmft_ed)


ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created
FPPFLAG =-fpp
endif
ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPFLAG =-cpp
endif


##$ REVISION SOFTWARE VARIABLES
REV=$(shell git rev-parse HEAD)
REV=$(shell git rev-parse HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90

all: FLAG:=${FFLAG}
all: MYINC:=$(DMFT_INC) $(GLOB_INC)
all: MYLIB:=$(DMFT_LIB) $(GLOB_LIB)
all:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${MYINC} ${MYLIB}
	@echo "Done"

debug: FLAG:=${DFLAG}
debug: MYINC:=$(DMFT_INC) $(GLOB_INC)
debug: MYLIB:=$(DMFT_LIB) $(GLOB_LIB)
debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${MYINC} ${MYLIB}
	@echo "Done"



lanc: FLAG:=${FFLAG}
lanc: MYINC:=$(LANC_INC) $(GLOB_INC)
lanc: MYLIB:=$(LANC_LIB) $(GLOB_LIB)
lanc:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${MYINC} ${MYLIB}
	@echo "Done"

lanc_debug: FLAG:=${DFLAG}
lanc_debug: MYINC:=$(LANC_INC) $(GLOB_INC)
lanc_debug: MYLIB:=$(LANC_LIB) $(GLOB_LIB)
lanc_debug:
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ")
	@echo ""
	$(FC) $(FLAG) $(EXE).f90 -o $(DIREXE)/$(EXE) ${MYINC} ${MYLIB}
	@echo "Done"


clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~
	@rm -fv  $(DIREXE)/$(EXE)



#########################################################################
