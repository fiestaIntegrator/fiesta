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

#ifndef SCANNER_H
#define SCANNER_H 1


#include "comdef.h"

#if ARITHMETIC == MPFR
#include <stdio.h>
#include <mpfr.h>
#define PREC_MIN MPFR_PREC_MIN
/*Just arbitrary restriction to fit signed int on 32-bit systems, 2^31-1:*/
#define PREC_MAX 2147483647
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define CONST_BUF 256
/* Parameters affectin the performance: */
/*If the number of triads is more then this value, the result 
will be stored by address rather than by value:*/
#define INDIRECT_ADDRESSING_THRESHOLD 300000

#define OFTUSED 2048
/*128/8=16, one we do not use, so it is reasonable to use 
 NREGS= n*16 - 1, n=0,1,2...:*/
#define NREGS 31


/*size of a sigle tline pool:*/
#define TLINE_POOL 2047

/*Numbers without fractions from 0 to MAX_TAB are tabbed: */

#define MAX_TAB 4

/* :Parameters affectin the performance */

/*pow(a,1)...pow(a,MAX_POW_INLINE-1) are inlined:*/
#define MAX_POW_INLINE 8

#define INBUFSIZE 4096

/*Size of a buffer - 1 (since last must be '\0'):*/
#define MULTILINESIZE 16383


#ifdef COM_INLINE
#define SC_INLINE COM_INLINE
#else
#define SC_INLINE inline
#endif

#ifdef COM_INT
#define SC_INT COM_INT
#else
#define SC_INT int
#endif

void halt(int retval,char *fmt, ...);

#include "triesi.h"
#include "collectstr.h"
#include "hash.h"

/*compile-time operations:*/
/*unary without result:*/
#define OP_JMP 1
/*unary with result:*/
/*Inline powers:*/
#define OP_IPOW2 2
#define OP_IPOW3 3
#define OP_IPOW4 4
#define OP_IPOW5 5
#define OP_IPOW6 6

#define OP_NEG MAX_POW_INLINE
#define OP_INV (MAX_POW_INLINE+1)
#define OP_LOG (MAX_POW_INLINE +2)
/* -- attention! Must be the biggest unary!:*/
#define OP_CPY (MAX_POW_INLINE +10)
/*binrary without result:*/
#define OP_LEQ (MAX_POW_INLINE +13)
#define OP_L (MAX_POW_INLINE +14)
#define OP_EQ (MAX_POW_INLINE +15)
#define OP_GEQ (MAX_POW_INLINE +16)
#define OP_G (MAX_POW_INLINE +17)
#define OP_NEQ (MAX_POW_INLINE +18)
/*binary with result:*/
#define OP_CPY2 (MAX_POW_INLINE +20)
#define OP_IPOW (MAX_POW_INLINE +21)
#define OP_POW (MAX_POW_INLINE +22)
#define OP_MINUS (MAX_POW_INLINE +23)
/* -- attention! Must be the biggest one!:*/
#define OP_DIV (MAX_POW_INLINE +24)
/*binary commuting:*/
#define OP_PLUS (-1)
#define OP_MUL (-2)

   /*
     Commuting (all binary):
      -1 -- OP_PLUS
      -2 - OP_MUL

     Non-commuting:

       unary: only first agrument is relevant:
       without result:
       1 -- OP_JMP
       with result: (000 100 001 101)
       MAX_POW_INLINE -- OP_NEG
       MAX_POW_INLINE+1 -- OP_INV
       MAX_POW_INLINE+2 -- OP_LOG
       (just copy the first operand to result:)
       MAX_POW_INLINE +10 -- OP_CPY 

       binary:
       without result: (000 100 010 110)
       MAX_POW_INLINE +13 -- <=
       MAX_POW_INLINE +14 -- <
       MAX_POW_INLINE +15 -- =
       MAX_POW_INLINE +16 -- >=
       MAX_POW_INLINE +17 -- >
       MAX_POW_INLINE +18 -- !=
       with result: (000 001 010 011 100 101 110 111)
       (copy the first operand to the second one and to the result:)
       MAX_POW_INLINE +20 -- OP_CPY2
       (ipow may only be with integer second operator,010 011 110 111:)
       MAX_POW_INLINE +21 -- OP_IPOW -- integer power
       MAX_POW_INLINE +22 -- OP_POW  -- non-integer power
       MAX_POW_INLINE +23 -- OP_MINUS
       MAX_POW_INLINE +24 -- OP_DIV

    */


/* 
Run-time operations.  The result may be stored directly in a triad,
or in some extrnal array.
*/
#define ROP_CPY2_P 1
#define ROP_CPY_P 2
#define ROP_CPY2_F 3
#define ROP_CPY_F 4
#define ROP_LOG_P 5
#define ROP_LOG_F 6
#define ROP_NEG_P 7
#define ROP_NEG_F 8
#define ROP_INV_P 9
#define ROP_INV_F 10
#define ROP_POW_P 11
#define ROP_POW_F 12
#define ROP_IPOW_P 13
#define ROP_IPOW_F 14
#define ROP_IPOW2_P 15
#define ROP_IPOW2_F 16
#define ROP_IPOW3_P 17
#define ROP_IPOW3_F 18
#define ROP_IPOW4_P 19
#define ROP_IPOW4_F 20
#define ROP_IPOW5_P 21
#define ROP_IPOW5_F 22
#define ROP_IPOW6_P 23
#define ROP_IPOW6_F 24
#define ROP_MINUS_P 25
#define ROP_MINUS_F 26
#define ROP_DIV_P 27
#define ROP_DIV_F 28
#define ROP_PLUS_P 29
#define ROP_PLUS_F 30
#define ROP_MUL_P 31
#define ROP_MUL_F 32
#define ROP_JMP 33
#define ROP_LEQ 34
#define ROP_L 35
#define ROP_EQ 36
#define ROP_GEQ 37
#define ROP_G 38
#define ROP_NEQ 39

struct scan_struct;
/*Compile Time triads:*/
/*This triad is allocated using mmap:*/
struct ct_triad_struct {
   /*Dynamically growing arrays:*/
   char *operation;  /*<= 0 -- commuting, > 0 -- non-commuting*/

   /*A triad type: 0 -- main cache, >0 -- one of an if-chain;
     -1 -- non- optimizable, like (the first) "if", "jpm", cpy,...;:*/
   SC_INT *triadType; 
   /* triadType is reused building rt-trads as ctNotUsed1: if >0, current triad is
   the last one which uses the result of thiad pointed out by this field*/

   SC_INT *lastUsing;/*last triad wich uses the result of this one*/

  /*How many other triads refer to the result of this one:*/
   SC_INT *refCounter;
   /* refCounter is reused building rt-trads as ctNotUsed2: if >0, current triad is
   the last one which uses the result of thiad pointed out by this field*/
   /* Operands:
      0 -- NULL, rezult of NOP, ignored
      >0 -- triads
      -nx..-1 -- x
      -nf-nx ... -nx-1 -- f
      < -nf-nx -- fline args
    */
   SC_INT *firstOperand; 
   SC_INT *secondOperand;

   SC_INT free;
   SC_INT max;

   struct scan_struct *theScan;
   hash_table_t hashTable;
};
typedef struct ct_triad_struct ct_triad_t;

/*Run time triads:*/
typedef struct rt_triad_addrP_struct {
   /*Not arrays! Just pointers:*/
   INTERNAL_FLOAT *firstOperand;
   INTERNAL_FLOAT *secondOperand;
   INTERNAL_FLOAT *result;
}rt_triadaddrP_t;

typedef struct rt_triad_addrIP_struct {
   INTERNAL_FLOAT *firstOperand;
   SC_INT secondOperand;
   INTERNAL_FLOAT *result;
}rt_triadaddrIP_t;

typedef struct rt_triad_addrF_struct {
   /*Not arrays! Just pointers:*/
   INTERNAL_FLOAT *firstOperand;
   INTERNAL_FLOAT *secondOperand;
   INTERNAL_FLOAT result;
}rt_triadaddrF_t;

#ifdef MIXED_ARITHMETIC
/*about MIXED_ARITHMETIC see comments in comdef.h */

/*The same as rt_triadaddrF_t but with native arithmetic:*/
typedef struct rt_triad_addrFN_struct {
   /*Not arrays! Just pointers:*/
   FLOAT *firstOperand;
   FLOAT *secondOperand;
   FLOAT result;
}rt_triadaddrFN_t;
/*The same as rt_triadaddrIF_t but with native arithmetic:*/
typedef struct rt_triad_addrIFN_struct {
   FLOAT *firstOperand;
   SC_INT secondOperand;
   FLOAT result;
}rt_triadaddrIFN_t;
#endif

typedef struct rt_triad_addrIF_struct {
   INTERNAL_FLOAT *firstOperand;
   SC_INT secondOperand;
   INTERNAL_FLOAT result;
}rt_triadaddrIF_t;

typedef   union rt_triad_addr_union {
   rt_triadaddrP_t  aP;
   rt_triadaddrIP_t aIP;
   rt_triadaddrF_t  aF;
#ifdef MIXED_ARITHMETIC
   rt_triadaddrFN_t  aFN;
   rt_triadaddrIFN_t  aIFN;
#endif
   rt_triadaddrIF_t aIF;
   SC_INT           aJ;
}rt_triadaddr_t;

/*This triad is allocated using malloc:*/
typedef struct rt_triad_struct {
   /*Allocated arrays:*/
   char *operations;
   rt_triadaddr_t *operands;
#ifdef MIXED_ARITHMETIC
   rt_triadaddr_t *nativeOperands;
#endif
   SC_INT length;
   int status; /*just a multipurpose flag*/
}rt_triad_t;

typedef struct multiLine_struct{
   collect_t mline;
   int nbuf;/*Translation counter*/
   int free;
}multiLine_t;

typedef struct scan_struct{
#ifndef NO_FILES
   char buf[INBUFSIZE+1];
#endif
   char *theChar;
   multiLine_t *multiLine;
   int fd;
   ct_triad_t ctTriad;/*Compile Time triads*/
   rt_triad_t rtTriad;/*RunTime triads*/

   int nx;/*dimesion of x*/
   int nf;/*number of f;*/
#ifdef MIXED_ARITHMETIC
   int *maxX;/* maxX[i] in a maximal negative power of a bare x[i]*/
#endif

   /*
      0..MAX_TAB-1 -tabbed constants
      MAX_TAB -- undefined (reserved?)
      MAX_TAB+1 ... 2*MAX_TAB-1 -- inverse tabbed (i.e., 1/n is in n+MAX_TAB cell)
      2*MAX_TAB ... 2*MAX_TAB+NCONST-1 -- predefined constants (like Pi)
      2*MAX_TAB+NCONST ... fline->fill-1 -- other constants
    */

   collectFloat_t *fline;
#ifdef MIXED_ARITHMETIC
   collectNativeFloat_t *fNativeLine;
#endif

   INTERNAL_FLOAT *x;
#ifdef MIXED_ARITHMETIC
   FLOAT *nativeX;
   /*Collects addresses of all activated MP variables:*/
   collect_t *allMPvariables;
#endif
   SC_INT *f;

   collect_t *tline;
   /*
     structure of tline:
     [0] the 128-bytes aligned array of NREGS+2 INTERNAL_FLOAT's
     [1] the address of the allocated block with this array
     [2] the array of OFTUSED+1 INTERNAL_FLOAT's
     [3] and further -- arrays of TLINE_POOL INTERNAL_FLOAT's.
    */

   FLOAT *ep;
   int *ep_i;
   int wasCut;
   /*
      pstack is a data stack.
      >=0 -- triads
      -nx..-1 -- x
      -nf-nx ... -nx -- f
      < -nf-nx -- fline args
    */
   collectInt_t *pstack;

   /*Symbolic constants like PolyLog and numerical strings:*/
   trieRoot_t newConstantsTrie;
#ifdef MIXED_ARITHMETIC
   /*indices of newConstantsTrie->mpool.pool[]:*/
   collectInt_t *constStrings;
   mpool_t newConstantsPool;
#endif
   int flags;

   /*if-related stuff:*/
   /* +/- conditions, + for the first "if" brunch, "-" -- for the second one:*/
   collectInt_t *condChain;
   /*Chain of all coundition leading to the current status:*/
   collectInt_t *allCondChain;
   /*Array  of stored chains:*/
   collect_t *allocatedCondChains;
}scan_t;

#ifndef NO_FILES
scan_t *newScanner(char *fname);
#endif

scan_t *newScannerMultiStr(multiLine_t *multiLine);
int initMultiLine(multiLine_t *multiLine);
int addToMultiLine(multiLine_t *multiLine, char *str);
void destroyMultiLine(multiLine_t *multiLine);

scan_t *newScannerFromStr(char *str);
void destroyScanner( scan_t *theScan);
int masterScan(scan_t *theScan);

#ifdef MIXED_ARITHMETIC

static SC_INLINE FLOAT ipow(FLOAT x,SC_INT y)
{
   FLOAT z, u;
#else
static SC_INLINE INTERNAL_FLOAT ipow(INTERNAL_FLOAT x,SC_INT y)
{
   INTERNAL_FLOAT z, u;
#endif
/*
   determines x to the power y
*/
   if ( y == 2 ) u = x*x;
   else {
      if ( ( y & 1 ) != 0 ) u = x;
      else u = 1.0L;
      z = x;
      y >>= 1;
      while ( y ) {
         z = z*z;
         if ( ( y & 1 ) != 0 ) u *= z;
         y >>= 1;
      }
   }

   return u;
}/*ipow*/

#ifdef __cplusplus
}
#endif
#endif
