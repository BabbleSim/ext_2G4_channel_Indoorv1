# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

BSIM_BASE_PATH?=$(abspath ../ )
include ${BSIM_BASE_PATH}/common/pre.make.inc

2G4_libPhyComv1_COMP_PATH?=$(abspath ${BSIM_COMPONENTS_PATH}/ext_2G4_libPhyComv1)
2G4_phy_v1_COMP_PATH?=$(abspath ${BSIM_COMPONENTS_PATH}/ext_2G4_phy_v1)

SRCS:=src/channel_2G4Indoorv1.c \
		  src/channel_2G4Indoorv1_argparse.c \
		  src/channel_2G4Indoorv1_pathloss.c

INCLUDES:= -I${libUtilv1_COMP_PATH}/src/ \
           -I${libRandv2_COMP_PATH}/src/ \
           -I${libPhyComv1_COMP_PATH}/src/ \
           -I${2G4_phy_v1_COMP_PATH}/src/ \
           -I${2G4_libPhyComv1_COMP_PATH}/src

LIB_NAME:=lib_2G4Channel_Indoorv1
A_LIBS:=
SO_LIBS:=

DEBUG:=-g
OPT:=
ARCH:=
WARNINGS:=-Wall -pedantic
COVERAGE:=
CFLAGS:=${ARCH} ${DEBUG} ${OPT} ${WARNINGS} -MMD -MP -std=c99 -fPIC ${INCLUDES}
LDFLAGS:=${ARCH} ${COVERAGE} -lm -lfftw3
CPPFLAGS:=-D_XOPEN_SOURCE=700

include ${BSIM_BASE_PATH}/common/make.lib_so.inc
