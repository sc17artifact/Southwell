#!/bin/bash

REMOVE = rm -f
SRC_DIR = src/
METLIB = libmetis.a
SOURCES = ${SRC_DIR}Misc.c \
	  ${SRC_DIR}FileUtils.c \
	  ${SRC_DIR}Laplace.c \
	  ${SRC_DIR}MatrixUtils.c \
 	  ${SRC_DIR}Multicolor.c \
	  ${SRC_DIR}pseudo_peripheral_node.c \

SEQ_SOURCES = ${SRC_DIR}SEQ_Setup.c \
	      ${SRC_DIR}SEQ_Solve.c \

SMEM_SOURCES = ${SRC_DIR}SMEM_Misc.c \
               ${SRC_DIR}SMEM_Setup.c \
               ${SRC_DIR}SMEM_Solve.c \

DMEM_SOURCES = ${SRC_DIR}DMEM_Misc.c \
               ${SRC_DIR}DMEM_Setup.c \
               ${SRC_DIR}DMEM_Solve.c \
	       ${SRC_DIR}DMEM_SolveUtils.c \
	       ${SRC_DIR}DMEM_PSouthUtils.c \
	       ${SRC_DIR}DMEM_DSouthUtils.c \
	       $(SRC_DIR)DMEM_MCGSutils.c \


foMPI_DIR  = ./foMPI/foMPI-0.2.1/
CASPER_DIR = ./casper/lib/
METIS_DIR  = ./metis/lib/
MKL_LINK =

foMPI_LINK = ${foMPI_DIR}fompi.ar -ldmapp ${foMPI_DIR}mpitypes/install/lib/libmpitypes.a ${foMPI_DIR}libtopodisc/libtopodisc.a


MPI_INTEL_COMPILE = mpiicpc -std=c++11 -O3 -qopenmp -mkl
MPI_CORI_COMPILE = CC -std=c++11 -O3 -qopenmp -mkl #-hfp3
MPI_GPP_COMPILE =  mpicxx -fopenmp -std=c++0x -O3

MPI_CORI_CASPER_COMPILE = ${MPI_CORI_COMPILE}
MPI_CORI_foMPI_COMPILE  = ${MPI_CORI_COMPILE} -DUSE_FOMPI -I${foMPI_DIR}

DMEM_COMPILE     = ${MPI_CORI_CASPER_COMPILE}
DMEM_EXTERN_LIBS = ${METIS_DIR}libmetis.a ${foMPI_LINK} ${CASPER_DIR}libcasper.a

OMP_GPP_COMPILE = g++ -fopenmp -std=c++0x -O3 #-qopt-report=5 -qopt-report-phase=vec
SMEM_COMPILE   = $(MPI_CORI_COMPILE)


DMEM_HOST: DMEM_CLEANHOST DMEM_Southwell

DMEM_MIC: DMEM_CLEANMIC DMEM_Southwell_mic

SMEM_HOST: cleanhost SMEM_Southwell

SMEM_MIC: cleanmic SMEM_SouthwellMic PrintVecMic VaryThreadsMic ResMic SweepMic

DMEM_Southwell: ${SRC_DIR}DMEM_Southwell.c
	$(DMEM_COMPILE) $(DMEM_SOURCES) $(SOURCES) ${SRC_DIR}DMEM_Southwell.c -o DMEM_Southwell ${DMEM_EXTERN_LIBS}

DMEM_Southwell_mic: ${SRC_DIR}DMEM_Southwell.c
	$(DMEM_COMPILE) -mmic $(DMEM_SOURCES) $(SOURCES) ${SRC_DIR}DMEM_Southwell.c -o DMEM_Southwell_mic ${MPI_MIC_LIBS_POSTFIX}

SMEM_Southwell: ${SRC_DIR}SMEM_Southwell.c
	$(SMEM_COMPILE) $(SMEM_SOURCES) $(SOURCES) ${SRC_DIR}SMEM_Southwell.c -o SMEM_Southwell libmetis.a

SMEM_SouthwellMic: ${SRC_DIR}SMEM_Southwell.c
	$(SMEM_COMPILE) -mmic $(SMEM_SOURCES) $(SOURCES) ${SRC_DIR}SMEM_Southwell.c -o SMEM_SouthwellMic libmetis_mic.a

SEQ_Southwell: ${SRC_DIR}SEQ_Southwell.c
	$(SMEM_COMPILE) $(SEQ_SOURCES) $(SOURCES) ${SRC_DIR}SEQ_Southwell.c -o SEQ_Southwell $(MKL_LINK)



PrintVec: PrintVec.c
	$(SMEM_COMPILE) PrintVec.c -o PrintVec

PrintVecMic: PrintVec.c
	$(SMEM_COMPILE) -mmic PrintVec.c -o PrintVecMic

VaryThreads: VaryThreads.c
	$(SMEM_COMPILE) $(SOURCES) VaryThreads.c -o VaryThreads libmetis.a

VaryThreadsMic: VaryThreads.c
	$(SMEM_COMPILE) -mmic $(SOURCES) VaryThreads.c -o VaryThreadsMic libmetis_mic.a

Res: Res.c
	$(SMEM_COMPILE) $(SOURCES) Res.c -o Res libmetis.a

ResMic: Res.c
	$(SMEM_COMPILE) -mmic $(SOURCES) Res.c -o ResMic libmetis_mic.a

Sweep: Sweep.c
	$(SMEM_COMPILE) $(SOURCES) Sweep.c -o Sweep libmetis.a

SweepMic: Sweep.c
	$(SMEM_COMPILE) -mmic $(SOURCES) Sweep.c -o SweepMic libmetis_mic.a

TextToBin: ${SRC_DIR}TextToBin.c
	$(SMEM_COMPILE) ${SRC_DIR}TextToBin.c -o TextToBin

DMEM_CLEANHOST:
	$(REMOVE) DMEM_Southwell

DMEM_CLEANMIC:
	$(REMOVE) DMEM_Southwell_mic


cleanhost :
	$(REMOVE) SMEM_Southwell PrintVec VaryThreads Res Sweep TextToBin

cleanmic :
	$(REMOVE) SMEM_SouthwellMic PrintVecMic VaryThreadsMic ResMic SweepMic


SEQ_CLEAN:
	$(REMOVE) SEQ_Southwell
