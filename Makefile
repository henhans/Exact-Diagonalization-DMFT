#=========================================================================
include sfmake.inc
FC=$(SFMPI)/mpif90
#FPP=MPI
#=========================================================================
#--> HUBBARD MODELS:
#EXE=ed_hm_bethe
#EXE=ed_hm_2dsquare
#EXE=ed_hm_bethe_mixgf

#--> PERIODIC ANDERSON & P-D MODELS
#EXE=ed_pam_1b
#EXE=ed_pam_2b
#EXE=ed_lda1b
#EXE=ed_lda

#--> B-H-Z MODELS
#EXE=ed_2x2bhz
EXE=ed_bhz

DIR =drivers
DIREXE=$(HOME)/.bin

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

#COMPILATION:
OBJS= MATRIX_SPARSE.o ED_EIGENSPACE.o ED_VARS_GLOBAL.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_AUX_FUNX.o ED_BATH.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o ED_CHI2FIT.o ED_DIAG.o DMFT_ED.o

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD) -fpp -D_$(FPP) #-openmp
all: ARGS=$(SFLIBS)
all:compile

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT) -fpp -D_$(FPP) #-openmp
opt: ARGS=$(SFLIBS)
opt:compile

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB) -fpp -D_$(FPP)
debug: ARGS=$(SFLIBS_DEB)
debug:compile

compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH)$(FPP) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)$(FPP)

.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) 

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)
