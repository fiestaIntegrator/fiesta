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

/*
  This file contains tools for manipulation by several types of 
  re-allocatable arrays ("vectors"):
  collect_t       -- array of generic void*;
  collectStr_t    -- a string;
  collectInt_t    -- array of integers;
  collectFloat_t -- array of INTERNAL_FLOATs.
  collectNativeFloat_t -- array of FLOATs, defined only when 
   the macro MIXED_ARITHMETIC  is defined. 
  All inline functions and public prototypes are in this file,
  linkable functions are placed to the file collectstr.c

  Note, this is not a universal toolkit, it is written specially
  for the program CIntegrate, and all the methods are not
  structured well.  For each vector we have: a constructor 
  (initCollect...),  a reallocator (reallocCollect...) and a 
  method for adding an element. For collect_t and collectFloat_t
  there is also a method for extracting the last added element.
  There are no any destructors since they are just 
  "free(buffer)".
*/

#ifndef COLLECTSTR_H
#define COLLECTSTR_H 1

#include "comdef.h"

#if ARITHMETIC == MPFR
#include <mpfr.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef COM_INLINE
#define CS_INLINE COM_INLINE
#else
#define CS_INLINE inline
#endif

#ifdef COM_INT
#define CS_INT COM_INT
#else
#define CS_INT int
#endif

/*Initial sizes:*/
#define INIT_STR_SIZE 128
#define INIT_INT_SIZE 128
#define INIT_FLOAT_SIZE 256

#ifndef FLOAT
#define FLOAT double
#endif

#include <stdlib.h>
#include <stdio.h>

struct collect_s{
    void **buf;
    size_t top;
    size_t fill;
};

struct collectStr_s{
    char *buf;
    size_t top;
    size_t fill;
};

struct collectInt_s{
    CS_INT *buf;
    size_t top;
    size_t fill;
};

struct collectFloat_s{
    INTERNAL_FLOAT *buf;
    size_t top;
    size_t fill;
};

#ifdef MIXED_ARITHMETIC
struct collectNativeFloat_s{
    FLOAT *buf;
    size_t top;
    size_t fill;
};

typedef struct collectNativeFloat_s collectNativeFloat_t;
#endif

typedef struct collect_s collect_t;
typedef struct collectStr_s collectStr_t;
typedef struct collectInt_s collectInt_t;
typedef struct collectFloat_s collectFloat_t;

/************ array of generic void*: ******************/
void reallocCollect(collect_t *theCollection );
collect_t *initCollect( collect_t *theCollection);

static CS_INLINE collect_t *pushCell( collect_t *theCollection, void *theCell )
{
   if( theCollection->fill >= theCollection->top )
     reallocCollect(theCollection);
   theCollection->buf[theCollection->fill++]=theCell;
   return theCollection;
}/*pushCell*/

static CS_INLINE void *popCell( collect_t *theCollection)
{
   if(theCollection->fill < 1)
      return NULL;/*stack underflow*/
   return theCollection->buf[--theCollection->fill];
}/*popCell*/
/************ :array of generic void* ******************/
/************ string (array of char): ******************/
void reallocCollectStr(collectStr_t *theSTring );
collectStr_t *initCollectStr( collectStr_t *theSTring );

static CS_INLINE collectStr_t *addChar( collectStr_t *theSTring, char c )
{
   if( theSTring->fill >= theSTring->top )
     reallocCollectStr(theSTring);
   theSTring->buf[theSTring->fill++]=c;
   return theSTring;
}/*addChar*/
/************ :string (array of char) ******************/
/************ array of int: ****************************/
void reallocCollectInt(collectInt_t *theInt );
collectInt_t *initCollectInt( collectInt_t *theInt );
static CS_INLINE collectInt_t *addInt( collectInt_t *theInt, CS_INT c )
{
   if( theInt->fill >= theInt->top )
     reallocCollectInt(theInt);
   theInt->buf[theInt->fill++]=c;
   return theInt;
}/*addInt*/
static CS_INLINE CS_INT popInt( collectInt_t *theInt)
{
   if(theInt->fill < 1)
      return 0;/*stack underflow*/
   return theInt->buf[--theInt->fill];
}/*popInt*/
/************ :array of int *****************************/

/************ array of INTERNAL_FLOAT:****************************/

static CS_INLINE void initIFvar(INTERNAL_FLOAT *f)
{
#if ARITHMETIC == NATIVE
   return;
#elif ARITHMETIC == MPFR
   mpfr_init(*f);
#endif
}


#if ARITHMETIC != NATIVE
static CS_INLINE void clearIFvar(INTERNAL_FLOAT *f)
{
#if ARITHMETIC == MPFR
   mpfr_clear(*f);
#endif
}
#endif



#if ARITHMETIC != NATIVE
static CS_INLINE void int2IFloat(INTERNAL_FLOAT *i, int c)
{
#if ARITHMETIC == MPFR
   mpfr_set_si(*i,c,GMP_RNDN);
#endif
}/*int2IFloat*/

static CS_INLINE void float2IFloat(INTERNAL_FLOAT *i, 
                                   FLOAT *c)
{
#if ARITHMETIC == MPFR
   mpfr_set_d(*i,*c,GMP_RNDN);
#endif

}/*float2IFloat*/
#endif

void reallocCollectFloat(collectFloat_t *theFloat );
collectFloat_t *initCollectFloat( size_t iniSize, collectFloat_t *theFloat );
static CS_INLINE collectFloat_t *addFloat( collectFloat_t *theFloat, FLOAT c)

{
   if( theFloat->fill >= theFloat->top )
     reallocCollectFloat(theFloat);
#if ARITHMETIC == NATIVE
      theFloat->buf[theFloat->fill++]=c;
#else
   initIFvar(theFloat->buf+theFloat->fill);
   float2IFloat(theFloat->buf+theFloat->fill,&c);
   theFloat->fill++;
#endif
   return theFloat;
}/*addFloat*/
/************ :array of INTERNAL_FLOAT ****************************/

/************ array of FLOAT:****************************/
#ifdef MIXED_ARITHMETIC

void reallocCollectNativeFloat(collectNativeFloat_t *theFloat );
collectNativeFloat_t *initCollectNativeFloat( size_t iniSize, collectNativeFloat_t *theFloat );
static CS_INLINE collectNativeFloat_t *addNativeFloat( collectNativeFloat_t *theFloat, FLOAT c)
{
   if( theFloat->fill >= theFloat->top )
     reallocCollectNativeFloat(theFloat);
   theFloat->buf[theFloat->fill++]=c;
   return theFloat;
}/*addNativeFloat*/
#endif
/************ :array of FLOAT ****************************/



#ifdef __cplusplus
}
#endif

#endif
