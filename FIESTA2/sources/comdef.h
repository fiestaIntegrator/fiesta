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

#ifndef COMDEF_H
#define COMDEF_H

#define VEGAS vegas_
#define VEGAS1 vegas1_
#define VEGAS2 vegas2_
#define VEGAS3 vegas3_
#define RESULT result_
#define COM_INLINE inline

#define  NO_FILES 1
#define  FOR_MATHEMATICA 1

#define WITH_OPTIMIZATION 1


/*
#ifdef WITH_OPTIMIZATION
#define STATISTICS_OUT 1
#endif
*/

/*#define COM_INT int*/
#define COM_INT long int


/*IEEE FP type, used for input-output:*/
#define FLOAT double

/*Could'n handle smaller values!*/
#define ABS_MONOM_MIN 1.0E-308

/*Initial precision for internal FP, in bits:*/
#ifndef PRECISION
#define PRECISION 64
#endif 


#define NATIVE 1
#define MPFR 4

/*The macro MIXED_ARITHMETIC forces FIESTA to compile two variants of 
  rt triples: the usual one with INTERNAL_FLOAT internal numbers and
  the second one with native FLOAT numbers. If one of x-variable is 
  smaller than some threshold, the function is evaluated using 
  INTERNAL_FLOAT, if all x-vars are big enough, the function is 
  evaluated using native FLOAT arithmetic.
  MIXED_ARITHMETIC requires ARITHMETIC != NATIVE.
  Note, when ARITHMETIC != NATIVE, all results are stored directly in 
  rt triples.*/
#ifndef MIXED_ARITHMETIC
#define MIXED_ARITHMETIC 1
#else

#if MIXED_ARITHMETIC == 0
#undef MIXED_ARITHMETIC
#endif

#endif

/*Generally speaking, this should come from the command line:*/
#ifndef ARITHMETIC
#define ARITHMETIC MPFR
#endif

#if ARITHMETIC==NATIVE
#undef PRECISION
#ifdef MIXED_ARITHMETIC
#undef MIXED_ARITHMETIC
#endif
#endif

/*FP type used for internal representation:*/
#if ARITHMETIC == NATIVE
#define INTERNAL_FLOAT FLOAT
#elif ARITHMETIC == MPFR
#define INTERNAL_FLOAT mpfr_t
#else
#error "Unknown ARITHMETIC"
#endif

#if ARITHMETIC != NATIVE
extern int g_default_precision;
extern int g_default_precisionChanged;

#endif
#ifdef MIXED_ARITHMETIC

#define MPSMALLX 0.001
extern FLOAT g_mpsmallX;
#define MPTHRESHOLD 1E-9
extern FLOAT g_mpthreshold;
#define MPMIN 1E-48
extern FLOAT g_mpmin;
extern int g_mpminChanged;
#define MPPRECISIONSHIFT 38
extern int g_mpPrecisionShift;
#endif
#endif
