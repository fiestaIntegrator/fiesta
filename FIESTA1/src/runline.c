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
static RL_INLINE FLOAT ipow(FLOAT x,SC_INT y)
{
   FLOAT z, u;

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


FLOAT runline(rt_triad_t *rt)
{
FLOAT a;
char *rto=rt->operations;
char *rtoStop=rt->operations+rt->length-1;
rt_triadaddr_t *rta=rt->operands;
   for(;rto<rtoStop;rto++,rta++){
      switch(*rto){
         case ROP_CPY2_P :
            *(rta->aP.secondOperand)=*(rta->aP.firstOperand);
            /*no break*/
         case ROP_CPY_P :
               *(rta->aP.result) = *(rta->aP.firstOperand);
            break;
         case ROP_CPY2_F :
            *(rta->aF.secondOperand)=*(rta->aF.firstOperand);
            /*no break*/
         case ROP_CPY_F :
               rta->aF.result = *(rta->aF.firstOperand);
            break;
         /*:CPY*/
         /*LOG:*/
         case ROP_LOG_P:
            *(rta->aP.result) = log(*(rta->aP.firstOperand));
            break;
         case ROP_LOG_F:
            rta->aF.result = log(*(rta->aF.firstOperand));
            break;
         /*:LOG*/
         /*NEG:*/
         case ROP_NEG_P :
            *(rta->aP.result) = -*(rta->aP.firstOperand);
            break;
         case ROP_NEG_F :
            rta->aF.result = - *(rta->aF.firstOperand);
            break;
         /*:NEG*/
         /*INV:*/
         case ROP_INV_P :
            *(rta->aP.result) = 1.0 / *(rta->aP.firstOperand);
            break;
         case ROP_INV_F :
            rta->aF.result = 1.0 / *(rta->aF.firstOperand);
            break;
         /*:INV*/
         /*POW:*/
         case ROP_POW_P :
            *(rta->aP.result) = pow(*(rta->aP.firstOperand),*(rta->aP.secondOperand));
            break;
         case ROP_POW_F :
            rta->aF.result = pow(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
            break;
         /*:POW*/
         /*IPOW:*/
         /*Integer power, a special case, the second argument is used just as it is:*/
         case ROP_IPOW_P:
            *(rta->aIP.result) = ipow(*(rta->aIP.firstOperand),rta->aIP.secondOperand);
            break;
         case ROP_IPOW_F :
            rta->aIF.result = ipow(*(rta->aIF.firstOperand),rta->aIF.secondOperand);
            break;
         /*:IPOW*/
         /*IPOW2:*/
         /*Integer power 2*/
         case ROP_IPOW2_P:
            *(rta->aP.result) = *(rta->aP.firstOperand) *  *(rta->aP.firstOperand);
            break;
         case ROP_IPOW2_F :
            rta->aF.result = *(rta->aF.firstOperand) *  *(rta->aF.firstOperand);
            break;
         /*:IPOW2*/
         /*Integer power 3*/
         case ROP_IPOW3_P:
            *(rta->aP.result) = *(rta->aP.firstOperand) *  
                                  *(rta->aP.firstOperand) * *(rta->aP.firstOperand);
            break;
         case ROP_IPOW3_F:
            rta->aF.result = *(rta->aF.firstOperand) *  
                                  *(rta->aF.firstOperand) * *(rta->aF.firstOperand);
            break;
         /*:IPOW3*/
         /*Integer power 4*/
         case ROP_IPOW4_P:
            a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
            *(rta->aP.result) = a*a;
            break;
         case ROP_IPOW4_F:
            a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
            rta->aF.result = a*a;
            break;
         /*:IPOW4*/
         /*Integer power 5*/
         case ROP_IPOW5_P:
            a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
            *(rta->aP.result) = a * a * *(rta->aP.firstOperand);
            break;
         case ROP_IPOW5_F:
            a=*(rta->aF.firstOperand) * *(rta->aP.firstOperand);
            rta->aF.result = a * a * *(rta->aF.firstOperand);
            break;
         /*:IPOW5*/
         /*Integer power 6*/
         case ROP_IPOW6_F:
           a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
           rta->aF.result = a*a*a;
            break;
         case ROP_IPOW6_P:
           a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
           *(rta->aP.result) = a*a*a;
            break;
         /*:IPOW6*/
         /*MINUS:*/
         case ROP_MINUS_P :
            *(rta->aP.result) = *(rta->aP.firstOperand)-*(rta->aP.secondOperand);
            break;
         case ROP_MINUS_F :
            rta->aF.result = *(rta->aF.firstOperand)-*(rta->aF.secondOperand);
            break;
         /*:MINUS*/
         /*DIV:*/
         case ROP_DIV_P :
            *(rta->aP.result) = *(rta->aP.firstOperand) / *(rta->aP.secondOperand);
            break;
         case ROP_DIV_F :
            rta->aF.result = *(rta->aF.firstOperand) / *(rta->aF.secondOperand);
            break;
         /*:DIV*/
         /*PLUS:*/
         case ROP_PLUS_P :
            *(rta->aP.result) = *(rta->aP.firstOperand)+*(rta->aP.secondOperand);
            break;
         case ROP_PLUS_F :
            rta->aF.result = *(rta->aF.firstOperand)+*(rta->aF.secondOperand);
            break;
         /*:PLUS*/
         /*MUL:*/
         case ROP_MUL_P :
            *(rta->aP.result) = *(rta->aP.firstOperand) * *(rta->aP.secondOperand);
            break;
         case ROP_MUL_F :
            rta->aF.result = *(rta->aF.firstOperand) * *(rta->aF.secondOperand);
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
            if( *(rta->aP.firstOperand) <= *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:LEQ*/

         /*L:*/
         case ROP_L:
            if( *(rta->aP.firstOperand) < *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
         /*:L*/
         /*EQ:*/
         case ROP_EQ:
            if( *(rta->aP.firstOperand) == *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:EQ*/
         /*GEQ:*/
         case ROP_GEQ:
            if( *(rta->aP.firstOperand) >= *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:GEQ*/
         /*G:*/
         case ROP_G:
            if( *(rta->aP.firstOperand) > *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:G*/
         /*NEQ:*/
         case ROP_NEQ:
            if( *(rta->aP.firstOperand) != *(rta->aP.secondOperand) ){
               rto++;rta++;
            }
            break;
         /*:NEQ*/
      }/*switch(*rto)*/
   }/*for(;rto<rtoStop;rto++,rta++)*/
   switch(*rto){
      case ROP_CPY2_P :
      case ROP_CPY_P :
         return *(rta->aP.firstOperand);
      case ROP_CPY2_F :
      case ROP_CPY_F :
         return *(rta->aF.firstOperand);
      case ROP_LOG_P:
         return log(*(rta->aP.firstOperand));
      case ROP_LOG_F:
         return log(*(rta->aF.firstOperand));
      case ROP_NEG_P :
         return -*(rta->aP.firstOperand);
      case ROP_NEG_F :
         return - *(rta->aF.firstOperand);
      case ROP_INV_P :
         return 1.0 / *(rta->aP.firstOperand);
      case ROP_INV_F :
         return 1.0 / *(rta->aF.firstOperand);
      case ROP_POW_P :
         return pow(*(rta->aP.firstOperand),*(rta->aP.secondOperand));
      case ROP_POW_F :
         return  pow(*(rta->aF.firstOperand),*(rta->aF.secondOperand));
      case ROP_IPOW_P :
         return ipow(*(rta->aIP.firstOperand),rta->aIP.secondOperand);
      case ROP_IPOW_F :
         return ipow(*(rta->aIF.firstOperand),rta->aIF.secondOperand);
      case ROP_IPOW2_P:
         return  *(rta->aP.firstOperand) *  *(rta->aP.firstOperand);
      case ROP_IPOW2_F :
         return *(rta->aF.firstOperand) *  *(rta->aF.firstOperand);
      case ROP_IPOW3_P:
         return *(rta->aP.firstOperand) *  
                *(rta->aP.firstOperand) * *(rta->aP.firstOperand);
      case ROP_IPOW3_F:
         return *(rta->aF.firstOperand) *  
                *(rta->aF.firstOperand) * *(rta->aF.firstOperand);
      case ROP_IPOW4_P:
         a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
         return a*a;
      case ROP_IPOW4_F:
         a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
         return a*a;
      case ROP_IPOW5_P:
         a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
         return a * a * *(rta->aP.firstOperand);
      case ROP_IPOW5_F:
         a=*(rta->aF.firstOperand) * *(rta->aP.firstOperand);
         return a * a * *(rta->aF.firstOperand);
      case ROP_IPOW6_F:
         a=*(rta->aF.firstOperand) * *(rta->aF.firstOperand);
         return a*a*a;
      case ROP_IPOW6_P:
        a=*(rta->aP.firstOperand) * *(rta->aP.firstOperand);
        return a*a*a;
      case ROP_MINUS_P :
         return *(rta->aP.firstOperand)-*(rta->aP.secondOperand);
      case ROP_MINUS_F :
         return *(rta->aF.firstOperand)-*(rta->aF.secondOperand);
      case ROP_DIV_P :
         return *(rta->aP.firstOperand) / *(rta->aP.secondOperand);
      case ROP_DIV_F :
         return *(rta->aF.firstOperand) / *(rta->aF.secondOperand);
      case ROP_PLUS_P :
         return *(rta->aP.firstOperand)+*(rta->aP.secondOperand);
      case ROP_PLUS_F :
         return *(rta->aF.firstOperand)+*(rta->aF.secondOperand);
      case ROP_MUL_P :
         return *(rta->aP.firstOperand) * *(rta->aP.secondOperand);
      case ROP_MUL_F :
         return *(rta->aF.firstOperand) * *(rta->aF.secondOperand);
   }/*switch(*rto)*/
   
   /*An error!*/
   halt (-10, "Runline error %d\n",*rto);
   return 0.0;
}/*runline*/

FLOAT runExpr(FLOAT *x,scan_t *theScan)
{
   /*copy x to the scratch:*/
   memcpy(theScan->x+1,x+1,sizeof(FLOAT)*(theScan->nx - 1));
   
   if(theScan->wasCut){
      int i;
      if(theScan->ep[0] > 0.0){/*"Global" cut*/
         for(i=1;i<theScan->nx; i++)if(x[i]<theScan->ep[0])
            theScan->x[i]=theScan->ep[0];
      }else /*"local" cut, imopse it individually for each x:*/
         for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i]){
            if(x[i]<theScan->ep[i]){
               theScan->x[i]=theScan->ep[i];
            }
         }/*for(i=1;i<theScan->nx; i++)if(theScan->ep_i[i])*/
   }/*if(theScan->wasCut)*/
   return runline(&theScan->rtTriad);
}/*runExpr*/
