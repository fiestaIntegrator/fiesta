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
In this file one can find an interpretor of the array of quadruples
("rt_triad_t", for historical reasons). The array is built by the
function buildRtTriad from the file scanner.c.

The interface function runExpr() makes some initializations and
invokes the interpreter runline().
*/

#include "scanner.h"
#include "runline.h"

#include <errno.h>
#include <string.h>
#include <math.h>

/*#define DEBUG 1*/
/*#define DEBUG_PREC 1*/

#if ARITHMETIC != NATIVE
static INTERNAL_FLOAT *l_one;
#endif

#if ARITHMETIC == NATIVE

#define CHC1(a)
#define CHC2(a,b)
#define CHC3(a,b,c) 

#define CPY(a,b) (a) = (b)
#define LOG(a,b) (a) = log(b)
#define NEG(a,b) (a) = -(b)
#define INV(a,b) (a) = 1.0 / (b)
#define POW(a,b,c) (a) = pow((b),(c))
#define IPOW(a,b,c) (a) = ipow((b),(c))
#define MUL(a,b,c) (a) = (b) * (c)
#define POW3(a,b) (a) = (b)*(b)*(b)
#define POW4(r,b) (*a)=(b)*(b);(r)=(*a)*(*a)
#define POW5(r,b) (*a)=(b)*(b);(r)=(*a)*(*a)*(b)
#define POW6(r,b) (*a)=(b)*(b);(r)=(*a)*(*a)*(*a)
#define MINUS(a,b,c) (a)=(b)-(c)
#define DIV(a,b,c) (a)=(b)/(c)
#define PLUS(a,b,c) (a)=(b)+(c)
#define LEQ(a,b) (a)<=(b)
#define LLL(a,b) (a)<(b)
#define EQ(a,b) (a) == (b)
#define GEQ(a,b) (a)>=(b)
#define GGG(a,b) (a) > (b)
#define NEQ(a,b) (a) != (b)

static RL_INLINE FLOAT IFloat2Float(INTERNAL_FLOAT *f)
{
return *f;
}/*IFloat2Float*/

#elif ARITHMETIC == MPFR

#ifdef DEBUG_PREC
#define CHC1(a) \
   if(mpfr_get_prec(a)<g_default_precision)fprintf(stderr,"%ld<%d\n",\
   mpfr_get_prec(a),g_default_precision);

#define CHC2(a,b) \
   if(mpfr_get_prec(a)<g_default_precision)fprintf(stderr,"%ld<%d\n",\
   mpfr_get_prec(a),g_default_precision);\
   if(mpfr_get_prec(b)<g_default_precision)fprintf(stderr,"%ld<%d\n",\
   mpfr_get_prec(b),g_default_precision);

#define CHC3(a,b,c) \
   if(mpfr_get_prec(a)<g_default_precision)fprintf(stderr,"%ld<%d\n",\
   mpfr_get_prec(a),g_default_precision);\
   if(mpfr_get_prec(b)<g_default_precision)fprintf(stderr,"%ld<%d\n",\
   mpfr_get_prec(b),g_default_precision);\
   if(mpfr_get_prec(c)<g_default_precision)fprintf(stderr,"%ld<%d\n",\
   mpfr_get_prec(c),g_default_precision);
#else
#define CHC1(a)
#define CHC2(a,b)
#define CHC3(a,b,c) 
#endif

#define CPY(a,b)  { CHC2(a,b) mpfr_set((a),(b),GMP_RNDN); }
#define LOG(a,b) { CHC1(a) if (mpfr_cmp_d((b),1.0E-100) <=0)\
                      mpfr_set_d((a),-230.258509299405,GMP_RNDN);\
                   else mpfr_log((a),(b),GMP_RNDN); }
#define NEG(a,b) { CHC2(a,b) mpfr_neg((a),(b),GMP_RNDN);}
#define INV(a,b) { CHC2(a,b) mpfr_div((a),*l_one,(b),GMP_RNDN);}
#define POW(a,b,c) {CHC3(a,b,c) mpfr_pow((a),(b),(c),GMP_RNDN);}
#define IPOW(a,b,c) {CHC2(a,b) mpfr_pow_si((a),(b),(c),GMP_RNDN);}
#define MUL(a,b,c) { CHC3(a,b,c) mpfr_mul((a),(b),(c),GMP_RNDN);}
#define POW3(a,b) { CHC1(a) mpfr_pow_ui((a),(b),3,GMP_RNDN);}
#define POW4(r,b) { CHC1(r) mpfr_pow_ui((r),(b),4,GMP_RNDN);}
#define POW5(r,b) { CHC1(r) mpfr_pow_ui((r),(b),5,GMP_RNDN);}
#define POW6(r,b) { CHC1(r) mpfr_pow_ui((r),(b),6,GMP_RNDN);}
#define MINUS(a,b,c) { CHC3(a,b,c) mpfr_sub((a),(b),(c),GMP_RNDN);}
#define DIV(a,b,c) { CHC3(a,b,c) mpfr_div((a),(b),(c),GMP_RNDN);}
#define PLUS(a,b,c) { CHC3(a,b,c) mpfr_add((a),(b),(c),GMP_RNDN);}
#define LEQ(a,b) (mpfr_cmp((a),(b))<=0)
#define LLL(a,b) (mpfr_cmp((a),(b))<0)
#define EQ(a,b) (mpfr_cmp((a),(b))==0)
#define GEQ(a,b) (mpfr_cmp((a),(b))>=0)
#define GGG(a,b) (mpfr_cmp((a),(b))>0)
#define NEQ(a,b) (mpfr_cmp((a),(b)))


static RL_INLINE FLOAT IFloat2Float(INTERNAL_FLOAT *f)
{
/*@@@*/ /*mpfr_fprintf(stderr,"%RE\n",*f);*/
return mpfr_get_d(*f,GMP_RNDN);
}/*IFloat2Float*/

#endif



/*Parametr INTERNAL_FLOAT a is just auxiliary variable, for non-native
  arithmetic it is used for the final result:*/
FLOAT runline(rt_triad_t *rt,INTERNAL_FLOAT *a)
{
char *rto=rt->operations;
char *rtoStop=rt->operations+rt->length-1;
rt_triadaddr_t *rta=rt->operands;

   for(;rto<rtoStop;rto++,rta++){
      switch(*rto){
         case ROP_CPY2_P :
            CPY(*(rta->aP.secondOperand),*(rta->aP.firstOperand));
            /*no break*/
         case ROP_CPY_P :
            CPY(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_CPY2_F :
            CPY(*(rta->aF.secondOperand),*(rta->aF.firstOperand));
            /*no break*/
         case ROP_CPY_F :
               CPY(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:CPY*/
         /*LOG:*/
         case ROP_LOG_P:
            LOG(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_LOG_F:
            LOG(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:LOG*/
         /*NEG:*/
         case ROP_NEG_P :
            NEG(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_NEG_F :
            NEG(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:NEG*/
         /*INV:*/
         case ROP_INV_P :
            INV(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_INV_F :
            INV(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:INV*/
         /*POW:*/
         case ROP_POW_P :
            POW(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_POW_F :
            POW(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
         /*:POW*/
         /*IPOW:*/
         /*Integer power, a special case, the second argument is used just as it is:*/
         case ROP_IPOW_P:
            IPOW(*(rta->aIP.result),*(rta->aIP.firstOperand),rta->aIP.secondOperand);
            break;
         case ROP_IPOW_F :
            IPOW(rta->aIF.result,*(rta->aIF.firstOperand),rta->aIF.secondOperand);
            break;
         /*:IPOW*/
         /*IPOW2:*/
         /*Integer power 2*/
         case ROP_IPOW2_P:
            MUL(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.firstOperand));
            break;
         case ROP_IPOW2_F :
            MUL(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.firstOperand));
            break;
         /*:IPOW2*/
         /*Integer power 3*/
         case ROP_IPOW3_P:
            POW3(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_IPOW3_F:
            POW3(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:IPOW3*/
         /*Integer power 4*/
         case ROP_IPOW4_P:
            POW4(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_IPOW4_F:
            POW4(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:IPOW4*/
         /*Integer power 5*/
         case ROP_IPOW5_P:
            POW5(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_IPOW5_F:
            POW5(rta->aF.result,*(rta->aF.firstOperand));
            break;
         /*:IPOW5*/
         /*Integer power 6*/
         case ROP_IPOW6_F:
           POW6(rta->aF.result,*(rta->aF.firstOperand));
            break;
         case ROP_IPOW6_P:
           POW6(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         /*:IPOW6*/
         /*MINUS:*/
         case ROP_MINUS_P :
            MINUS(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_MINUS_F :
            MINUS(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
         /*:MINUS*/
         /*DIV:*/
         case ROP_DIV_P :
            DIV(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_DIV_F :
            DIV(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
         /*:DIV*/
         /*PLUS:*/
         case ROP_PLUS_P :
            PLUS(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_PLUS_F :
            PLUS(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
         /*:PLUS*/
         /*MUL:*/
         case ROP_MUL_P :
            MUL(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_MUL_F :
            MUL(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
        /*:MUL*/
        /*JMP:*/
         case ROP_JMP:
               rto+=rta->aJ;
               rta+=rta->aJ;
            break;
         /*:JMP*/
         /*LEQ:*/
         case ROP_LEQ:
            CHC2(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            if( LEQ(*(rta->aF.firstOperand), *(rta->aF.secondOperand)) ){
               rto++;rta++;
            }
            break;
         /*:LEQ*/
         /*L:*/
         case ROP_L:
            CHC2(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            if( LLL(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
               rto++;rta++;
            }
            break;
         /*:L*/
         /*EQ:*/
         case ROP_EQ:
            CHC2(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            if( EQ(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
               rto++;rta++;
            }
            break;
         /*:EQ*/
         /*GEQ:*/
         case ROP_GEQ:
            CHC2(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            if( GEQ(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
               rto++;rta++;
            }
            break;
         /*:GEQ*/
         /*G:*/
         case ROP_G:
            CHC2(*(rta->aF.firstOperand),*(rta->aF.secondOperand));

            if( GGG(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
               rto++;rta++;
            }
            break;
         /*:G*/
         /*NEQ:*/
         case ROP_NEQ:
            CHC2(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            if( NEQ(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
               rto++;rta++;
            }
            break;
         /*:NEQ*/
      }/*switch(*rto)*/
   }/*for(;rto<rtoStop;rto++,rta++)*/

   switch(*rto){
      case ROP_CPY2_P :
      case ROP_CPY_P :
         return IFloat2Float(rta->aP.firstOperand);
      case ROP_CPY2_F :
      case ROP_CPY_F :
         return IFloat2Float(rta->aF.firstOperand);
      case ROP_LOG_P:
         LOG(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_LOG_F:
         LOG(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_NEG_P :
         NEG(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_NEG_F :
         NEG(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_INV_P :
         INV(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_INV_F :
         INV(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_POW_P :
         POW(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
      case ROP_POW_F :
         POW(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
      case ROP_IPOW_P :
         IPOW(*a,*(rta->aIP.firstOperand),rta->aIP.secondOperand);
         return IFloat2Float(a);
      case ROP_IPOW_F :
         IPOW(*a,*(rta->aIF.firstOperand),rta->aIF.secondOperand);
         return IFloat2Float(a);
      case ROP_IPOW2_P:
         MUL(*a,*(rta->aP.firstOperand),*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW2_F :
         MUL(*a,*(rta->aF.firstOperand),*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW3_P:
         POW3(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW3_F:
         POW3(*a, *(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW4_P:
         POW4(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW4_F:
         POW4(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW5_P:
         POW5(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW5_F:
         POW6(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW6_F:
         POW6(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
      case ROP_IPOW6_P:
        POW6(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
      case ROP_MINUS_P :
         MINUS(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
      case ROP_MINUS_F :
         MINUS(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
      case ROP_DIV_P :
         DIV(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
      case ROP_DIV_F :
         DIV(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
      case ROP_PLUS_P :
         PLUS(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
      case ROP_PLUS_F :
         PLUS(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
      case ROP_MUL_P :
         MUL(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
      case ROP_MUL_F :
         MUL(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
   }/*switch(*rto)*/
   
   /*An error!*/
   halt (-10, "Runline error %d\n",*rto);
   return 0.0;
}/*runline*/

#ifdef MIXED_ARITHMETIC
FLOAT runlineNative(rt_triad_t *rt)
{
char *rto=rt->operations;
char *rtoStop=rt->operations+rt->length-1;
rt_triadaddr_t *rta=rt->nativeOperands;
FLOAT a;
   for(;rto<rtoStop;rto++,rta++){
      switch(*rto){
         case ROP_CPY2_F :
            *(rta->aFN.secondOperand)=*(rta->aFN.firstOperand);
            /*no break*/
         case ROP_CPY_F :
            rta->aFN.result = *(rta->aFN.firstOperand);
            break;
         /*:CPY*/
         /*LOG:*/
         case ROP_LOG_F:
            if(*(rta->aFN.firstOperand) <= 1.0E-100l )
               rta->aFN.result = -230.258509299405l;
            else
               rta->aFN.result = log(*(rta->aFN.firstOperand));

            /*rta->aFN.result = log(*(rta->aFN.firstOperand));*/
            break;
         /*:LOG*/
         /*NEG:*/
         case ROP_NEG_F :
            rta->aFN.result = - *(rta->aFN.firstOperand);
            break;
         /*:NEG*/
         /*INV:*/
         case ROP_INV_F :
            rta->aFN.result = 1.0 / *(rta->aFN.firstOperand);
            break;
         /*:INV*/
         /*POW:*/
         case ROP_POW_F :
            rta->aFN.result = pow(*(rta->aFN.firstOperand),*(rta->aFN.secondOperand));
            break;
         /*:POW*/
         /*IPOW:*/
         /*Integer power, a special case, the second argument is used just as it is:*/
         case ROP_IPOW_F :
            rta->aIFN.result = ipow(*(rta->aIFN.firstOperand),rta->aIFN.secondOperand);
            break;
         /*:IPOW*/
         /*IPOW2:*/
         /*Integer power 2*/
         case ROP_IPOW2_F :
            rta->aFN.result = *(rta->aFN.firstOperand) *  *(rta->aFN.firstOperand);
            break;
         /*:IPOW2*/
         /*Integer power 3*/
         case ROP_IPOW3_F:
            rta->aFN.result = *(rta->aFN.firstOperand) *  
                                  *(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
            break;
         /*:IPOW3*/
         /*Integer power 4*/
         case ROP_IPOW4_F:
            a=*(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
            rta->aFN.result = a*a;
            break;
         /*:IPOW4*/
         /*Integer power 5*/
         case ROP_IPOW5_F:
            a=*(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
            rta->aFN.result = a * a * *(rta->aFN.firstOperand);
            break;
         /*:IPOW5*/
         case ROP_IPOW6_F:
           a=*(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
           rta->aFN.result = a*a*a;
            break;
         /*:IPOW6*/
         /*MINUS:*/
         case ROP_MINUS_F :
            rta->aFN.result = *(rta->aFN.firstOperand)-*(rta->aFN.secondOperand);
            break;
         /*:MINUS*/
         /*DIV:*/
         case ROP_DIV_F :
            rta->aFN.result = *(rta->aFN.firstOperand) / *(rta->aFN.secondOperand);
            break;
         /*:DIV*/
         /*PLUS:*/
         case ROP_PLUS_F :
            rta->aFN.result = *(rta->aFN.firstOperand)+*(rta->aFN.secondOperand);
            break;
         /*:PLUS*/
         /*MUL:*/
         case ROP_MUL_F :
            rta->aFN.result = *(rta->aFN.firstOperand) * *(rta->aFN.secondOperand);
            break;
        /*:MUL*/
        /*JMP:*/
         case ROP_JMP:
               rto+=rta->aJ;
               rta+=rta->aJ;
            break;
         /*:JMP*/
         /*LEQ:*/
         case ROP_LEQ:
            if( *(rta->aFN.firstOperand) <= *(rta->aFN.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:LEQ*/
         /*L:*/
         case ROP_L:
            if( *(rta->aFN.firstOperand) < *(rta->aFN.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:L*/
         /*EQ:*/
         case ROP_EQ:
            if( *(rta->aFN.firstOperand) == *(rta->aFN.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:EQ*/
         /*GEQ:*/
         case ROP_GEQ:
            if( *(rta->aFN.firstOperand) >= *(rta->aFN.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:GEQ*/
         /*G:*/
         case ROP_G:
            if( *(rta->aFN.firstOperand) > *(rta->aFN.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:G*/
         /*NEQ:*/
         case ROP_NEQ:
            if( *(rta->aFN.firstOperand) != *(rta->aFN.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:NEQ*/
      }/*switch(*rto)*/
   }/*for(;rto<rtoStop;rto++,rta++)*/

   switch(*rto){
      case ROP_CPY2_F :
      case ROP_CPY_F :
         return *(rta->aFN.firstOperand);
      case ROP_LOG_F:
            if(*(rta->aFN.firstOperand) <= 1.0E-100l )
               return -230.258509299405l;
            else
               return log(*(rta->aFN.firstOperand));
      case ROP_NEG_F :
         return - *(rta->aFN.firstOperand);
      case ROP_INV_F :
         return 1.0 / *(rta->aFN.firstOperand);
      case ROP_POW_F :
         return  pow(*(rta->aFN.firstOperand),*(rta->aFN.secondOperand));
      case ROP_IPOW_F :
         return ipow(*(rta->aIFN.firstOperand),rta->aIFN.secondOperand);
      case ROP_IPOW2_F :
         return *(rta->aFN.firstOperand) *  *(rta->aFN.firstOperand);
      case ROP_IPOW3_F:
         return *(rta->aFN.firstOperand) *  
                *(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
      case ROP_IPOW4_F:
         a=*(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
         return a*a;
      case ROP_IPOW5_F:
         a=*(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
         return a * a * *(rta->aFN.firstOperand);
      case ROP_IPOW6_F:
         a=*(rta->aFN.firstOperand) * *(rta->aFN.firstOperand);
         return a*a*a;
      case ROP_MINUS_F :
         return *(rta->aFN.firstOperand)-*(rta->aFN.secondOperand);
      case ROP_DIV_F :
         return *(rta->aFN.firstOperand) / *(rta->aFN.secondOperand);
      case ROP_PLUS_F :
         return *(rta->aFN.firstOperand)+*(rta->aFN.secondOperand);
      case ROP_MUL_F :
         return *(rta->aFN.firstOperand) * *(rta->aFN.secondOperand);
   }/*switch(*rto)*/
   
   /*An error!*/
   halt (-10, "Runline error %d\n",*rto);
   return 0.0;
}/*runlineNative*/

FLOAT l_prod;   

static RL_INLINE int isArithmeticNative(FLOAT *x,scan_t *theScan)
{
int i;
   l_prod=1.0;
   if(theScan->maxX[0]>0){
      for(i=1;i<theScan->nx; i++){
         if(theScan->maxX[i]>0) l_prod*=ipow(x[i],theScan->maxX[i]);
      }
      if(l_prod>=g_mpthreshold)
         return 1;
     return 0;
   }/*if(theScan->maxX[0]>0)*/
   return 1;
}/*isArithmeticNative*/

static RL_INLINE void createTruncatedX(FLOAT *x,scan_t *theScan)
{
FLOAT prod, mpsmallX=g_mpsmallX/10.0l,ep=ABS_MONOM_MIN/10.0l;
int i,n=0;

   if(mpsmallX < ABS_MONOM_MIN)
      mpsmallX = ABS_MONOM_MIN;

   for (i=1;i<theScan->nx;i++)
      if(x[i]<mpsmallX)
         theScan->nativeX[i]=(mpsmallX-x[i])/2.0;
      else
         theScan->nativeX[i]=0.0l;

   do{
      n++;
      prod=1.0;
      for(i=1;i<theScan->nx; i++)
         prod*=ipow(x[i]+theScan->nativeX[i],theScan->maxX[i]);
      if(prod < ABS_MONOM_MIN)
         for(i=1;i<theScan->nx; i++)
            theScan->nativeX[i]+= theScan->nativeX[i]/2.0l;
      else
         for(i=1;i<theScan->nx; i++)
            theScan->nativeX[i]-= theScan->nativeX[i]/2.0l;
      if(n==128)
         break;
   }while( (prod<ABS_MONOM_MIN)||(prod-ABS_MONOM_MIN > ep) );

   for(i=1;i<theScan->nx; i++){
      FLOAT tmp=x[i]+theScan->nativeX[i];
      float2IFloat(theScan->x+i,&tmp);
   }
   l_prod=prod;
}/*createTruncatedX*/

static RL_INLINE void resetConstants(scan_t *theScan)
{
   int i,j;
   for(j=i=0; i <theScan->fNativeLine->fill; i++){
      char *str=(char*) (theScan->newConstantsPool.pool + theScan->constStrings->buf[i]);
#if ARITHMETIC == MPFR
      mpfr_set_str(theScan->fline->buf[j++],str,10,GMP_RNDN);
#endif
      if(theScan->constStrings->buf[i+1] == -1){
         /*Store also */
         /*fline->buf[1] is 1, fline->buf[fline->fill-1] is 'val', we
           store 1/val into fline->buf[fline->fill]:*/
#if ARITHMETIC == MPFR
         mpfr_div(theScan->fline->buf[j],
               theScan->fline->buf[1],
               theScan->fline->buf[j-1],GMP_RNDN);
         j++;
#endif
         i++;/*Skip next "-1" in theScan->constStrings*/
      }/*if(theScan->constStrings->buf[i+1] == -1)*/
   }/*for(i=0; i <theScan->fNativeline->fill; i++)*/
}/*resetConstants*/

static RL_INLINE int resetPrecision(scan_t *theScan)
{
int i;
   mpfr_set_default_prec(g_default_precision);
   for(i=theScan->allMPvariables->fill-1; i>=0; i--)
      mpfr_prec_round( *((mpfr_t *)(theScan->allMPvariables->buf[i])),g_default_precision,GMP_RNDN);
   resetConstants(theScan);
   return 0;
}/*resetPrecision*/

static RL_INLINE int intlog2(FLOAT x){return ceil(log(x)*1.44269504088896l);}


FLOAT runExpr(FLOAT *x,scan_t *theScan)
{
FLOAT res;
int i;
int useNativeArithmetic;
static int l_default_precision_mem=0;;
   x--;/*The patch, for historical reasons the rest of the function 
         assumes that the array indices are started from 1, not 0*/
   useNativeArithmetic=isArithmeticNative(x,theScan);
   /*copy x to the scratch:*/
   if(useNativeArithmetic)
      memcpy(theScan->nativeX+1,x+1,sizeof(FLOAT)*(theScan->nx - 1));
   else
   {
      if( l_prod < ABS_MONOM_MIN)
         createTruncatedX(x,theScan);
         /*now l_prod is re-calculated*/

      if(l_prod>=g_mpmin){
         if(l_default_precision_mem){
            /*Restore the default precision*/
            g_default_precision=l_default_precision_mem;
            l_default_precision_mem=0;
            resetPrecision(theScan);
         }
      }/*if(l_prod>=g_mpmin)*/
      else
      {
         int newPrec=intlog2(1.0l/l_prod)+g_mpPrecisionShift;
         if(l_default_precision_mem){
            if(g_default_precision<newPrec){
               g_default_precision=newPrec;
               resetPrecision(theScan);
            }
         }/*if(l_default_precision_mem)*/
         else{
            l_default_precision_mem=g_default_precision;
            g_default_precision=newPrec;
            resetPrecision(theScan);
         }/*if(l_default_precision_mem)...else*/
      }/*if(l_prod>=g_mpmin)...else*/
   }/*if(useNativeArithmetic)...else*/

   for(i=1;i<theScan->nx;i++)
       float2IFloat(theScan->x+i,x+i);

   if(theScan->wasCut){
      if(theScan->ep[0] > 0.0l){/*"Global" cut*/
         for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0]){
            if(useNativeArithmetic)
               theScan->nativeX[i]=theScan->ep[0];
            /*else -- all cuts for MP mode are done*/

         }/*for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0])*/
      }else /*"local" cut, imopse it individually for each x:*/
         for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i]){
            if(x[i]<theScan->ep[i]){
               if(useNativeArithmetic)
                  theScan->nativeX[i]=theScan->ep[i];
               /*else -- all cuts for MP mode are done*/
            }
         }/*for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i])*/
   }/*if(theScan->wasCut)*/

      if(useNativeArithmetic){
         res=runlineNative(&theScan->rtTriad);
      }
      else{
         INTERNAL_FLOAT a;
         initIFvar(&a);
         l_one=theScan->fline->buf+1;
         res=runline(&theScan->rtTriad,&a);
      }
   return res;
}/*runExpr*/

#else
/*No MIXED_ARITHMETIC*/
FLOAT runExpr(FLOAT *x,scan_t *theScan)
{
FLOAT res;
int i;
   x--;/*The patch, for historical reasons the rest of the function 
         assumes that the array indices are started from 1, not 0*/
   /*copy x to the scratch:*/
#if ARITHMETIC == NATIVE
   memcpy(theScan->x+1,x+1,sizeof(FLOAT)*(theScan->nx - 1));
#else
   {/*block*/
         for(i=1;i<theScan->nx;i++)
            float2IFloat(theScan->x+i,x+i);
   }/*block*/
#endif

   if(theScan->wasCut){
      if(theScan->ep[0] > 0.0){/*"Global" cut*/
         for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0]){
#if ARITHMETIC == NATIVE
            theScan->x[i]=theScan->ep[0];
#else
            float2IFloat(theScan->x+i,theScan->ep);
#endif
         }/*for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0])*/
      }else /*"local" cut, imopse it individually for each x:*/
         for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i]){
            if(x[i]<theScan->ep[i]){
#if ARITHMETIC == NATIVE
               theScan->x[i]=theScan->ep[i];
#else
               float2IFloat(theScan->x+i,theScan->ep+i);
#endif
            }
         }/*for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i])*/
   }/*if(theScan->wasCut)*/

   {/*block*/
      INTERNAL_FLOAT a;
#if ARITHMETIC != NATIVE
      initIFvar(&a);
      l_one=theScan->fline->buf+1;
#endif
      res=runline(&theScan->rtTriad,&a);
#if ARITHMETIC != NATIVE
      clearIFvar(&a);
#endif
   }/*block*/

   return res;
}/*runExpr*/

#endif
