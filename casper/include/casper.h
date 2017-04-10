/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 * (C) 2015 by Argonne National Laboratory.
 *     See COPYRIGHT in top-level directory.
 */

#ifndef CASPER_H_INCLUDED
#define CASPER_H_INCLUDED

/* CASPER_VERSION is the version string. CASPER_NUMVERSION is the
 * numeric version that can be used in numeric comparisons.
 *
 * CASPER_VERSION uses the following format:
 * Version: [MAJ].[MIN].[REV][EXT][EXT_NUMBER]
 * Example: 1.0.7rc1 has
 *          MAJ = 1
 *          MIN = 0
 *          REV = 7
 *          EXT = rc
 *          EXT_NUMBER = 1
 *
 * CASPER_NUMVERSION will convert EXT to a format number:
 *          ALPHA (a) = 0
 *          BETA (b)  = 1
 *          RC (rc)   = 2
 *          PATCH (p) = 3
 * Regular releases are treated as patch 0
 *
 * Numeric version will have 1 digit for MAJ, 2 digits for MIN, 2
 * digits for REV, 1 digit for EXT and 2 digits for EXT_NUMBER. So,
 * 1.0.7rc1 will have the numeric version 10007201.
 */

#define CASPER_VERSION "1.0b1"
#define CASPER_NUMVERSION 10000101
#define CASPER_RELEASE_DATE "2016-11-12 23:49:21 -0600"
#define CASPER_BUILD_INFO "--prefix=/global/u2/j/jwolfson/GitSouthwell/extern_libs/cray-mpich/casper/ CC=cc CC=cc CFLAGS=  LDFLAGS= LIBS= CPPFLAGS=    LT_SYS_LIBRARY_PATH= CPP=cc -E"

/* Get the number of ghost processes. */
int CSP_ghost_size(int *ng);

#endif /* CASPER_H_INCLUDED */
