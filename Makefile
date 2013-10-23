#=========================================================================
include sfmake.inc
#FC=$(SFMPI)/mpif90
#=========================================================================
#EXE=ed_lda1b
#EXE=ed_pam_test
EXE=ed_hm_bethe
DIR =drivers
DIREXE=$(HOME)/.bin

BRANCH=$(shell git rev-parse --abbrev-ref HEAD)

#COMPILATION:
OBJS= MATRIX_SPARSE.o ARPACK_LANCZOS.o PLAIN_LANCZOS.o ED_EIGENSPACE.o ED_VARS_GLOBAL.o ED_AUX_FUNX.o ED_BATH.o ED_HAMILTONIAN.o ED_GREENS_FUNCTIONS.o ED_OBSERVABLES.o ED_CHI2FIT.o ED_DIAG.o DMFT_ED.o

#=================STANDARD COMPILATION====================================
all: FLAG=$(STD) -openmp
all: ARGS=$(SFLIBS)
all:compile

#================OPTIMIZED COMPILATION====================================
opt: FLAG=$(OPT) -openmp
opt: ARGS=$(SFLIBS)
opt:compile

#================DEBUGGIN COMPILATION=====================================
debug: FLAG=$(DEB)
debug: ARGS=$(SFLIBS_DEB)
debug:compile


compile: version $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

.f90.o:	
	$(FC) $(FLAG) -c $< $(SFINCLUDE) 

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)
