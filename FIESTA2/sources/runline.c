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

#if ARITHMETIC != NATIVE
static INTERNAL_FLOAT *l_one;
#endif

#if ARITHMETIC == NATIVE

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

#define CPY(a,b) mpfr_set((a),(b),GMP_RNDN)
#define LOG(a,b) { if (mpfr_cmp_d((b),1.0E-100) <=0)\
                      mpfr_set_d((a),-230.258509299405,GMP_RNDN);\
                   else mpfr_log((a),(b),GMP_RNDN); }
#define NEG(a,b) mpfr_neg((a),(b),GMP_RNDN)
#define INV(a,b) mpfr_div((a),*l_one,(b),GMP_RNDN)
#define POW(a,b,c) mpfr_pow((a),(b),(c),GMP_RNDN)
#define IPOW(a,b,c) mpfr_pow_si((a),(b),(c),GMP_RNDN)
#define MUL(a,b,c) mpfr_mul((a),(b),(c),GMP_RNDN)
#define POW3(a,b) mpfr_pow_ui((a),(b),3,GMP_RNDN)
#define POW4(r,b) mpfr_pow_ui((r),(b),4,GMP_RNDN)
#define POW5(r,b) mpfr_pow_ui((r),(b),5,GMP_RNDN)
#define POW6(r,b) mpfr_pow_ui((r),(b),6,GMP_RNDN)
#define MINUS(a,b,c) mpfr_sub((a),(b),(c),GMP_RNDN)
#define DIV(a,b,c) mpfr_div((a),(b),(c),GMP_RNDN)
#define PLUS(a,b,c) mpfr_add((a),(b),(c),GMP_RNDN)
#define LEQ(a,b) (mpfr_cmp((a),(b))<=0)
#define LLL(a,b) (mpfr_cmp((a),(b))<0)
#define EQ(a,b) (mpfr_cmp((a),(b))==0)
#define GEQ(a,b) (mpfr_cmp((a),(b))>=0)
#define GGG(a,b) (mpfr_cmp((a),(b))>0)
#define NEQ(a,b) mpfr_cmp((a),(b))

static RL_INLINE FLOAT IFloat2Float(INTERNAL_FLOAT *f)
{
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
//            *(rta->aP.secondOperand)=*(rta->aP.firstOperand);
            CPY(*(rta->aP.secondOperand),*(rta->aP.firstOperand));
            /*no break*/
         case ROP_CPY_P :
//               *(rta->aP.result) = *(rta->aP.firstOperand);
            CPY(*(rta->aP.result),*(rta->aP.firstOperand));
            break;
         case ROP_CPY2_F :
            CPY(*(rta->aF.secondOperand),*(rta->aF.firstOperand));
//            *(rta->aF.secondOperand)=*(rta->aF.firstOperand);
            /*no break*/
         case ROP_CPY_F :
               CPY(rta->aF.result,*(rta->aF.firstOperand));
//               rta->aF.result = *(rta->aF.firstOperand);
            break;
         /*:CPY*/
         /*LOG:*/
         case ROP_LOG_P:
            LOG(*(rta->aP.result),*(rta->aP.firstOperand));
//            *(rta->aP.result) = log(*(rta->aP.firstOperand));
            break;
         case ROP_LOG_F:
            LOG(rta->aF.result,*(rta->aF.firstOperand));
//            rta->aF.result = log(*(rta->aF.firstOperand));
            break;
         /*:LOG*/
         /*NEG:*/
         case ROP_NEG_P :
            NEG(*(rta->aP.result),*(rta->aP.firstOperand));
//            *(rta->aP.result) = -*(rta->aP.firstOperand);
            break;
         case ROP_NEG_F :
            NEG(rta->aF.result,*(rta->aF.firstOperand));
//            rta->aF.result = - *(rta->aF.firstOperand);
            break;
         /*:NEG*/
         /*INV:*/
         case ROP_INV_P :
            INV(*(rta->aP.result),*(rta->aP.firstOperand));
//            *(rta->aP.result) = 1.0 / *(rta->aP.firstOperand);
            break;
         case ROP_INV_F :
            INV(rta->aF.result,*(rta->aF.firstOperand));
//            rta->aF.result = 1.0 / *(rta->aF.firstOperand);
            break;
         /*:INV*/
         /*POW:*/
         case ROP_POW_P :
            POW(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
//            *(rta->aP.result) = pow(*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_POW_F :
            POW(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
//            rta->aF.result = pow(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
         /*:POW*/
         /*IPOW:*/
         /*Integer power, a special case, the second argument is used just as it is:*/
         case ROP_IPOW_P:
            IPOW(*(rta->aIP.result),*(rta->aIP.firstOperand),rta->aIP.secondOperand);
//            *(rta->aIP.result) = ipow(*(rta->aIP.firstOperand),rta->aIP.secondOperand);
            break;
         case ROP_IPOW_F :
            IPOW(rta->aIF.result,*(rta->aIF.firstOperand),rta->aIF.secondOperand);
//            rta->aIF.result = ipow(*(rta->aIF.firstOperand),rta->aIF.secondOperand);
            break;
         /*:IPOW*/
         /*IPOW2:*/
         /*Integer power 2*/
         case ROP_IPOW2_P:
            MUL(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.firstOperand));
//            *(rta->aP.result) = *(rta->aP.firstOperand) *  *(rta->aP.firstOperand);
            break;
         case ROP_IPOW2_F :
            MUL(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.firstOperand));
//            rta->aF.result = *(rta->aF.firstOperand) *  *(rta->aF.firstOperand);
            break;
         /*:IPOW2*/
         /*Integer power 3*/
         case ROP_IPOW3_P:
            POW3(*(rta->aP.result),*(rta->aP.firstOperand));
//            *(rta->aP.result) = *(rta->aP.firstOperand) *
//                                  *(rta->aP.firstOperand) * *(rta->aP.firstOperand);
            break;
         case ROP_IPOW3_F:
            POW3(rta->aF.result,*(rta->aF.firstOperand));
//            rta->aF.result = *(rta->aF.firstOperand) *  
//                                  *(rta->aF.firstOperand) * *(rta->aF.firstOperand);
            break;
         /*:IPOW3*/
         /*Integer power 4*/
         case ROP_IPOW4_P:
            POW4(*(rta->aP.result),*(rta->aP.firstOperand));
//            a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
//            *(rta->aP.result) = a*a;
            break;
         case ROP_IPOW4_F:
            POW4(rta->aF.result,*(rta->aF.firstOperand));
//            a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
//            rta->aF.result = a*a;
            break;
         /*:IPOW4*/
         /*Integer power 5*/
         case ROP_IPOW5_P:
            POW5(*(rta->aP.result),*(rta->aP.firstOperand));
//            a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
//            *(rta->aP.result) = a * a * *(rta->aP.firstOperand);
            break;
         case ROP_IPOW5_F:
            POW5(rta->aF.result,*(rta->aF.firstOperand));
//            a=*(rta->aF.firstOperand) * *(rta->aP.firstOperand);
//            rta->aF.result = a * a * *(rta->aF.firstOperand);
            break;
         /*:IPOW5*/
         /*Integer power 6*/
         case ROP_IPOW6_F:
           POW6(rta->aF.result,*(rta->aF.firstOperand));
//           a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
//           rta->aF.result = a*a*a;
            break;
         case ROP_IPOW6_P:
           POW6(*(rta->aP.result),*(rta->aP.firstOperand));
//           a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
//           *(rta->aP.result) = a*a*a;
            break;
         /*:IPOW6*/
         /*MINUS:*/
         case ROP_MINUS_P :
            MINUS(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
//            *(rta->aP.result) = *(rta->aP.firstOperand)-*(rta->aP.secondOperand);
            break;
         case ROP_MINUS_F :
            MINUS(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
//            rta->aF.result = *(rta->aF.firstOperand)-*(rta->aF.secondOperand);
            break;
         /*:MINUS*/
         /*DIV:*/
         case ROP_DIV_P :
            DIV(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
//            *(rta->aP.result) = *(rta->aP.firstOperand) / *(rta->aP.secondOperand);
            break;
         case ROP_DIV_F :
            DIV(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
//            rta->aF.result = *(rta->aF.firstOperand) / *(rta->aF.secondOperand);
            break;
         /*:DIV*/
         /*PLUS:*/
         case ROP_PLUS_P :
            PLUS(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
//            *(rta->aP.result) = *(rta->aP.firstOperand)+*(rta->aP.secondOperand);
            break;
         case ROP_PLUS_F :
            PLUS(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
//            rta->aF.result = *(rta->aF.firstOperand)+*(rta->aF.secondOperand);
            break;
         /*:PLUS*/
         /*MUL:*/
         case ROP_MUL_P :
            MUL(*(rta->aP.result),*(rta->aP.firstOperand),*(rta->aP.secondOperand));
//            *(rta->aP.result) = *(rta->aP.firstOperand) * *(rta->aP.secondOperand);
            break;
         case ROP_MUL_F :
            MUL(rta->aF.result,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
//            rta->aF.result = *(rta->aF.firstOperand) * *(rta->aF.secondOperand);
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
            if( LEQ(*(rta->aF.firstOperand), *(rta->aF.secondOperand)) ){
//            if( *(rta->aP.firstOperand) <= *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:LEQ*/
         /*L:*/
         case ROP_L:
            if( LLL(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
//            if( *(rta->aP.firstOperand) < *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:L*/
         /*EQ:*/
         case ROP_EQ:
            if( EQ(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
//            if( *(rta->aP.firstOperand) == *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:EQ*/
         /*GEQ:*/
         case ROP_GEQ:
            if( GEQ(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
//            if( *(rta->aP.firstOperand) >= *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:GEQ*/
         /*G:*/
         case ROP_G:
            if( GGG(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
//            if( *(rta->aP.firstOperand) > *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:G*/
         /*NEQ:*/
         case ROP_NEQ:
            if( NEQ(*(rta->aF.firstOperand),*(rta->aF.secondOperand)) ){
//            if( *(rta->aP.firstOperand) != *(rta->aP.secondOperand) ){
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
//         return *(rta->aP.firstOperand);
      case ROP_CPY2_F :
      case ROP_CPY_F :
         return IFloat2Float(rta->aF.firstOperand);
//         return *(rta->aF.firstOperand);
      case ROP_LOG_P:
         LOG(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         return log(*(rta->aP.firstOperand));
      case ROP_LOG_F:
         LOG(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         return log(*(rta->aF.firstOperand));
      case ROP_NEG_P :
         NEG(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         return -*(rta->aP.firstOperand);
      case ROP_NEG_F :
         NEG(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         return - *(rta->aF.firstOperand);
      case ROP_INV_P :
         INV(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         return 1.0 / *(rta->aP.firstOperand);
      case ROP_INV_F :
         INV(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         return 1.0 / *(rta->aF.firstOperand);
      case ROP_POW_P :
         POW(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
//         return pow(*(rta->aP.firstOperand),*(rta->aP.secondOperand));
      case ROP_POW_F :
         POW(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
//         return  pow(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
      case ROP_IPOW_P :
         IPOW(*a,*(rta->aIP.firstOperand),rta->aIP.secondOperand);
         return IFloat2Float(a);
//         return ipow(*(rta->aIP.firstOperand),rta->aIP.secondOperand);
      case ROP_IPOW_F :
         IPOW(*a,*(rta->aIF.firstOperand),rta->aIF.secondOperand);
         return IFloat2Float(a);
//         return ipow(*(rta->aIF.firstOperand),rta->aIF.secondOperand);
      case ROP_IPOW2_P:
         MUL(*a,*(rta->aP.firstOperand),*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         return  *(rta->aP.firstOperand) *  *(rta->aP.firstOperand);
      case ROP_IPOW2_F :
         MUL(*a,*(rta->aF.firstOperand),*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         return *(rta->aF.firstOperand) *  *(rta->aF.firstOperand);
      case ROP_IPOW3_P:
         POW3(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         return *(rta->aP.firstOperand) *  
//                *(rta->aP.firstOperand) * *(rta->aP.firstOperand);
      case ROP_IPOW3_F:
         POW3(*a, *(rta->aF.firstOperand));
         return IFloat2Float(a);
//         return *(rta->aF.firstOperand) *  
//                *(rta->aF.firstOperand) * *(rta->aF.firstOperand);
      case ROP_IPOW4_P:
         POW4(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
//         return a*a;
      case ROP_IPOW4_F:
         POW4(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
//         return a*a;
      case ROP_IPOW5_P:
         POW5(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//         a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
//         return a * a * *(rta->aP.firstOperand);
      case ROP_IPOW5_F:
         POW6(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         a=*(rta->aF.firstOperand) * *(rta->aP.firstOperand);
//         return a * a * *(rta->aF.firstOperand);
      case ROP_IPOW6_F:
         POW6(*a,*(rta->aF.firstOperand));
         return IFloat2Float(a);
//         a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
//         return a*a*a;
      case ROP_IPOW6_P:
        POW6(*a,*(rta->aP.firstOperand));
         return IFloat2Float(a);
//        a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
//        return a*a*a;
      case ROP_MINUS_P :
         MINUS(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aP.firstOperand)-*(rta->aP.secondOperand);
      case ROP_MINUS_F :
         MINUS(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aF.firstOperand)-*(rta->aF.secondOperand);
      case ROP_DIV_P :
         DIV(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aP.firstOperand) / *(rta->aP.secondOperand);
      case ROP_DIV_F :
         DIV(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aF.firstOperand) / *(rta->aF.secondOperand);
      case ROP_PLUS_P :
         PLUS(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aP.firstOperand)+*(rta->aP.secondOperand);
      case ROP_PLUS_F :
         PLUS(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aF.firstOperand)+*(rta->aF.secondOperand);
      case ROP_MUL_P :
         MUL(*a,*(rta->aP.firstOperand),*(rta->aP.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aP.firstOperand) * *(rta->aP.secondOperand);
      case ROP_MUL_F :
         MUL(*a,*(rta->aF.firstOperand),*(rta->aF.secondOperand));
         return IFloat2Float(a);
//         return *(rta->aF.firstOperand) * *(rta->aF.secondOperand);
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

            if(*(rta->aFN.firstOperand) <= 1.0E-10l )
               rta->aFN.result =-23.0258509299405l;
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

            if(*(rta->aFN.firstOperand) <= 1.0E-10l )
               return -23.0258509299405l;
            else
               return log(*(rta->aFN.firstOperand));
//         return log(*(rta->aFN.firstOperand));
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
      l_prod=1.0;
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
FLOAT prod, ep=g_mpmin/10.0;
int i,n=0;

   for (i=1;i<theScan->nx;i++)
      if(x[i]<g_mpsmallX)
         theScan->nativeX[i]=(g_mpsmallX-x[i])/2.0;
      else
         theScan->nativeX[i]=0.0;

   do{
      n++;
      prod=1.0;
      for(i=1;i<theScan->nx; i++)
         prod*=ipow(x[i]+theScan->nativeX[i],theScan->maxX[i]);
      if(prod < g_mpmin)
         for(i=1;i<theScan->nx; i++)
            theScan->nativeX[i]+= theScan->nativeX[i]/2.0;
      else
         for(i=1;i<theScan->nx; i++)
            theScan->nativeX[i]-= theScan->nativeX[i]/2.0;
      if(n==128)
         break;
   }while( (prod<g_mpmin)||(prod-g_mpmin > ep) );

   for(i=1;i<theScan->nx; i++){
      FLOAT tmp=x[i]+theScan->nativeX[i];
      float2IFloat(theScan->x+i,&tmp);
   }
#if 0
/*@@@*/
fprintf(stderr,"g_mpmin=%.15lf g_mpsmallX=%.15lf ",g_mpmin,g_mpsmallX);
for(i=1;i<theScan->nx; i++){
      FLOAT tmp=x[i]+theScan->nativeX[i];
      fprintf(stderr,"x[%d]=%.15lf rez[%d]=%.15lf ;",i,x[i],i,tmp);
}
fprintf(stderr,"\ncount=%d\n",n);
#endif
}/*createTruncatedX*/

#endif

FLOAT runExpr(FLOAT *x,scan_t *theScan)
{
FLOAT res;
#ifdef MIXED_ARITHMETIC
  int useNativeArithmetic;
#endif
   x--;/*The patch, for historical reasons the rest of the function 
         assumes that the array indices are started from 1, not 0*/
#ifdef MIXED_ARITHMETIC
   useNativeArithmetic=isArithmeticNative(x,theScan);
#endif
   /*copy x to the scratch:*/
#if ARITHMETIC == NATIVE
   memcpy(theScan->x+1,x+1,sizeof(FLOAT)*(theScan->nx - 1));
#else
#ifdef MIXED_ARITHMETIC
   if(useNativeArithmetic)
      memcpy(theScan->nativeX+1,x+1,sizeof(FLOAT)*(theScan->nx - 1));
   else
#endif
   {/*block*/
      int i;
#ifdef MIXED_ARITHMETIC
      if(l_prod>=g_mpmin)
#endif
         for(i=1;i<theScan->nx;i++)
            float2IFloat(theScan->x+i,x+i);
#ifdef MIXED_ARITHMETIC
      else
         createTruncatedX(x,theScan);
#endif
   }/*block*/
#endif

   if(theScan->wasCut){
      int i;
      if(theScan->ep[0] > 0.0){/*"Global" cut*/
         for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0]){
#if ARITHMETIC == NATIVE
            theScan->x[i]=theScan->ep[0];
#else
#ifdef MIXED_ARITHMETIC
            if(useNativeArithmetic)
               theScan->nativeX[i]=theScan->ep[0];
            /*else -- all cuts for MP mode are done*/
#else
            float2IFloat(theScan->x+i,theScan->ep);
#endif
#endif
         }/*for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0])*/
      }else /*"local" cut, imopse it individually for each x:*/
         for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i]){
            if(x[i]<theScan->ep[i]){
#if ARITHMETIC == NATIVE
               theScan->x[i]=theScan->ep[i];
#else
#ifdef MIXED_ARITHMETIC
               if(useNativeArithmetic)
                  theScan->nativeX[i]=theScan->ep[i];
               /*else -- all cuts for MP mode are done*/
#else
               float2IFloat(theScan->x+i,theScan->ep+i);
#endif
#endif
            }
         }/*for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i])*/
   }/*if(theScan->wasCut)*/

#ifdef MIXED_ARITHMETIC
      if(useNativeArithmetic)
         res=runlineNative(&theScan->rtTriad);
      else
#endif   
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
#if 0
/*@@@*/
{
int i;
   for(i=1;i<theScan->nx; i++)
      fprintf(stderr,"%.12lf ",x[i]);
fprintf(stderr,"res=%.12lf\n",res);
}
#endif
   return res;
}/*runExpr*/
