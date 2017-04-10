#!/bin/bash

METLIB         = libmetis.a
SMEM_SOURCES   = Misc.c \
		 FileUtils.c \
		 Laplace.c \
		 MatrixUtils.c \
		 Multicolor.c \
		 SEQ_Solve.c \
		 SMEM_Misc.c \
		 SMEM_Setup.c \
		 SMEM_Solve.c 
COMPILE        = g++ -fopenmp -std=c++0x -O3 #-qopt-report=5 -qopt-report-phase=vec 
REMOVE         = rm -f

host: cleanhost SMEM_Southwell

mic: cleanmic SMEM_SouthwellMic PrintVecMic VaryThreadsMic ResMic SweepMic

SMEM_Southwell: SMEM_Southwell.c
	$(COMPILE) $(SMEM_SOURCES) SMEM_Southwell.c -o SMEM_Southwell libmetis.a

SMEM_SouthwellMic: SMEM_Southwell.c
	$(COMPILE) -mmic $(SMEM_SOURCES) SMEM_Southwell.c -o SMEM_SouthwellMic libmetis_mic.a

PrintVec: PrintVec.c
	$(COMPILE) PrintVec.c -o PrintVec

PrintVecMic: PrintVec.c
	$(COMPILE) -mmic PrintVec.c -o PrintVecMic

VaryThreads: VaryThreads.c
	$(COMPILE) $(SMEM_SOURCES) VaryThreads.c -o VaryThreads libmetis.a

VaryThreadsMic: VaryThreads.c
	$(COMPILE) -mmic $(SMEM_SOURCES) VaryThreads.c -o VaryThreadsMic libmetis_mic.a

Res: Res.c
	$(COMPILE) $(SMEM_SOURCES) Res.c -o Res libmetis.a

ResMic: Res.c
	$(COMPILE) -mmic $(SMEM_SOURCES) Res.c -o ResMic libmetis_mic.a

Sweep: Sweep.c
	$(COMPILE) $(SMEM_SOURCES) Sweep.c -o Sweep libmetis.a

SweepMic: Sweep.c
	$(COMPILE) -mmic $(SMEM_SOURCES) Sweep.c -o SweepMic libmetis_mic.a

TextToBin: TextToBin.c
	$(COMPILE) TextToBin.c -o TextToBin

cleanhost :
	$(REMOVE) SMEM_Southwell PrintVec VaryThreads Res Sweep TextToBin

cleanmic :
	$(REMOVE) SMEM_SouthwellMic PrintVecMic VaryThreadsMic ResMic SweepMic
