EXE =fulled_AFpam_bethe
DIR =models_interface
DIREXE=$(HOME)/.bin

#########################################################################
include $(SFDIR)/etc/lib.mk
#########################################################################

#COMPILATION:
OBJS     = \
CGFIT.o \
ED_VARS_GLOBAL.o \
ED_AUX_FUNX.o \
ED_GETH.o \
ED_GETGF.o \
ED_GETOBS.o \
ED_CHI2FIT.o \
ED_DIAG.o \
DMFT_FULLED.o

OBJS_OPT = \
CGFIT_OPT.o \
ED_VARS_GLOBAL_OPT.o \
ED_AUX_FUNX_OPT.o \
ED_GETH_OPT.o \
ED_GETGF_OPT.o \
ED_GETOBS_OPT.o \
ED_CHI2FIT_OPT.o \
ED_DIAG_OPT.o \
DMFT_FULLED_OPT.o

OBJS_DEB = \
CGFIT_DEB.o \
ED_VARS_GLOBAL_DEB.o \
ED_AUX_FUNX_DEB.o \
ED_GETH_DEB.o \
ED_GETGF_DEB.o \
ED_GETOBS_DEB.o \
ED_CHI2FIT_DEB.o \
ED_DIAG_DEB.o \
DMFT_FULLED_DEB.o

#=================STANDARD COMPILATION====================================
all: 	version $(OBJS)
	@echo " ........... compile: normal ........... "
	$(FC) $(STD) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(SFMODS) $(SFLIBS)
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)
	@echo " ...................... done .............................. "


#=================OPTIMIZED COMPILATION====================================
opt: 	version $(OBJS_OPT)
	@echo " ........... compile: optimized ........... "
	$(FC) $(OPT) $(OBJS_OPT) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) $(SFMODS) $(SFLIBS)
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)
	@echo " ...................... done .............................. "


#================DEBUGGIN COMPILATION=====================================
debug:	version $(OBJS_DEB)	
	@echo " ........... compile : debug   ........... "
	${FC} ${DEB} ${OBJS_DEB} $(DIR)/${EXE}.f90 -o ${DIREXE}/${EXE} ${SFMODS} ${SFLIBS}
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)
	@echo " ...................... done .............................. "



ED_VARS_GLOBAL.o: ED_VARS_GLOBAL.f90
	$(FC) $(STD) -c ED_VARS_GLOBAL.f90 $(SFMODS)
ED_AUX_FUNX.o: ED_AUX_FUNX.f90
	$(FC) $(STD) -c ED_AUX_FUNX.f90 $(SFMODS)
ED_DIAG.o: ED_DIAG.f90
	$(FC) $(STD) -c ED_DIAG.f90 $(SFMODS)
ED_GETH.o: ED_GETH.f90
	$(FC) $(STD) -c ED_GETH.f90 $(SFMODS)
ED_GETGF.o: ED_GETGF.f90
	$(FC) $(STD)  -c ED_GETGF.f90 $(SFMODS)
ED_GETOBS.o: ED_GETOBS.f90
	$(FC) $(STD) -c ED_GETOBS.f90 $(SFMODS)
ED_CHI2FIT.o: ED_CHI2FIT.f90
	$(FC) $(STD) -c ED_CHI2FIT.f90 $(SFMODS)
CGFIT.o: CGFIT.f90
	$(FC) $(STD) -c CGFIT.f90 $(SFMODS)
DMFT_FULLED.o: DMFT_FULLED.f90
	$(FC) $(STD) -c DMFT_FULLED.f90 $(SFMODS)

ED_VARS_GLOBAL_OPT.o: ED_VARS_GLOBAL.f90
	$(FC) $(OPT) -c ED_VARS_GLOBAL.f90 $(SFMODS) -o ED_VARS_GLOBAL_OPT.o
ED_AUX_FUNX_OPT.o: ED_AUX_FUNX.f90
	$(FC) $(OPT) -c ED_AUX_FUNX.f90 $(SFMODS) -o ED_AUX_FUNX_OPT.o
ED_DIAG_OPT.o: ED_DIAG.f90
	$(FC) $(OPT) -c ED_DIAG.f90 $(SFMODS) -o ED_DIAG_OPT.o
ED_GETH_OPT.o: ED_GETH.f90
	$(FC) $(OPT) -c ED_GETH.f90 $(SFMODS) -o ED_GETH_OPT.o
ED_GETGF_OPT.o: ED_GETGF.f90
	$(FC) $(OPT)  -c ED_GETGF.f90 $(SFMODS) -o ED_GETGF_OPT.o
ED_GETOBS_OPT.o: ED_GETOBS.f90
	$(FC) $(OPT) -c ED_GETOBS.f90 $(SFMODS) -o ED_GETOBS_OPT.o
ED_CHI2FIT_OPT.o: ED_CHI2FIT.f90
	$(FC) $(OPT) -c ED_CHI2FIT.f90 $(SFMODS) -o ED_CHI2FIT_OPT.o
CGFIT_OPT.o: CGFIT.f90
	$(FC) $(OPT) -c CGFIT.f90 $(SFMODS) -o CGFIT_OPT.o
DMFT_FULLED_OPT.o: DMFT_FULLED.f90
	$(FC) $(OPT) -c DMFT_FULLED.f90 $(SFMODS) -o DMFT_FULLED_OPT.o

ED_VARS_GLOBAL_DEB.o: ED_VARS_GLOBAL.f90
	$(FC) $(DEB) -c ED_VARS_GLOBAL.f90 $(SFMODS) -o ED_VARS_GLOBAL_DEB.o
ED_AUX_FUNX_DEB.o: ED_AUX_FUNX.f90
	$(FC) $(DEB) -c ED_AUX_FUNX.f90 $(SFMODS) -o ED_AUX_FUNX_DEB.o
ED_DIAG_DEB.o: ED_DIAG.f90
	$(FC) $(DEB) -c ED_DIAG.f90 $(SFMODS) -o ED_DIAG_DEB.o
ED_GETH_DEB.o: ED_GETH.f90
	$(FC) $(DEB) -c ED_GETH.f90 $(SFMODS) -o ED_GETH_DEB.o
ED_GETGF_DEB.o: ED_GETGF.f90
	$(FC) $(DEB)  -c ED_GETGF.f90 $(SFMODS) -o ED_GETGF_DEB.o
ED_GETOBS_DEB.o: ED_GETOBS.f90
	$(FC) $(DEB) -c ED_GETOBS.f90 $(SFMODS) -o ED_GETOBS_DEB.o
ED_CHI2FIT_DEB.o: ED_CHI2FIT.f90
	$(FC) $(DEB) -c ED_CHI2FIT.f90 $(SFMODS) -o ED_CHI2FIT_DEB.o
CGFIT_DEB.o: CGFIT.f90
	$(FC) $(DEB) -c CGFIT.f90 $(SFMODS) -o CGFIT_DEB.o
DMFT_FULLED_DEB.o: DMFT_FULLED.f90
	$(FC) $(DEB) -c DMFT_FULLED.f90 $(SFMODS) -o DMFT_FULLED_DEB.o


clean: 
	@echo 'removing *.mod *.o *~'
	@rm -f *.mod *.o *~ revision.inc

#########################################################################
include  $(SFDIR)/etc/version.mk
#########################################################################
