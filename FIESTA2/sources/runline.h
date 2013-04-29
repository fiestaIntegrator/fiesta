/*
    Copyright (C) Alexander Smirnov and Mikhail Tentyukov. 
    This file is part of the program CIntegrate.
    The program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License version 2 as
    published by the Free Software Foundation.

    The program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/

#ifndef RUNLINE_H
#define RUNLINE_H 1

#include "comdef.h"
#include "scanner.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef COM_INLINE
#define RL_INLINE COM_INLINE
#else
#define RL_INLINE inline
#endif
FLOAT runline(rt_triad_t *rt,INTERNAL_FLOAT *a);
FLOAT runExpr(FLOAT *x,scan_t *theScan);
#ifdef __cplusplus
}
#endif
#endif
