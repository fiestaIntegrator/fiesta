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
The main file for the CIntegrate program.

The program is not supposed to start directly, it is
launched from Mathematica via the Mathlink protocol.

The routine AddString() builds the integrand step by step collecting
the incoming lines. Every moment the process may be canceled by
invoking the function ClearString().  After (almost) all of the lines
are collected, the routine Integrate() adds the rest and invokes the
translator masterScan(), file "scanner.c".  After the incoming
expression is translated into the tetrad array the routine Integrate()
invokes the VEGAS routine for numberOfIterations1 repetitions with
numberOfCalls1 sampling points in order to ``warmup'' the
grid. Finally, the VEGAS routine is invoked by the Integrate() routine
for numberOfIterations2 repetitions with numberOfCalls2 sampling
points for the real integration. Parameters numberOfCalls1,
numberOfIterations1, numberOfCalls2 and numberOfIterations2 could be
changed by the function setpoints().

The integrand for the VEGAS routine, the function theIntegrand(),
invokes the interpreter runExpr(), see the file "runline.c".

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mathlink.h>
#include <math.h>

#include "scanner.h"
#include "runline.h"

#include "vegas.h"

#define CUT_MULTIPLIER 10.0

extern void Integrate(char* s);

#ifdef STATISTICS_OUT
extern long int allT;
extern long int newT;
#endif

int setpoints(int i1,int i2,int i11,int i22);

void PutErrorMessage(char* function,char* message) {
			fprintf(stderr, "%s: %s\n",function,message);
			MLPutFunction(stdlink,"CompoundExpression",2);
			MLPutFunction(stdlink,"Message",2);
			MLPutFunction(stdlink,"MessageName",2);
			MLPutSymbol(stdlink, function);
			MLPutString(stdlink, "failed");
			MLPutString(stdlink, message);
			MLPutReal(stdlink,0.0);
}

static scan_t *fscan;

static multiLine_t ml;


FLOAT theIntegrand( FLOAT *x, FLOAT *wgt)
{
   return runExpr(--x,fscan);
}/*theIntegrand*/

static int numberOfCalls1 = 10000;
static int numberOfIterations1 = 5;
static int numberOfCalls2 = 100000;
static int numberOfIterations2 = 15;

int setpoints(int i1,int i2,int i11,int i22)
{
   numberOfCalls1=i1;
   numberOfIterations1=i2;
   numberOfCalls2=i11;
   numberOfIterations2=i22;

   return 0;
}/*setpoints*/

FLOAT *xCuts=NULL;
int nCuts=0;
int wasCut=0;

int setcut(int x, FLOAT ep)
{
   if(ep >=1.0){
      PutErrorMessage("SetCut","Ep >=1!");
      return 1;
   }
   if(x>=nCuts){
      FLOAT *tmp=calloc((x+1),sizeof(FLOAT));
      *tmp=-1.0;
      if(nCuts>0){
          memcpy(tmp,xCuts,nCuts*sizeof(FLOAT));
          free(xCuts);
      }
      xCuts=tmp;
      nCuts=x+1;
   }/*if(x>=nCuts)*/
   xCuts[x]=ep;
   if(ep < 0.0){/*clear x cut*/
      int i;
      wasCut=1;
      for(i=0;i<nCuts;i++)
         if(xCuts[i]>0.0)
            return 0;
      wasCut=0;
      return 0;
   }/*if(ep < 0.0)*/
   wasCut=1;
   return 0;
}/*setcut*/

int doAddString (char* s)
{
  if( ml.mline.buf == NULL)
     if(initMultiLine(&ml)){
        destroyMultiLine(&ml);
        return 1;
     }/*if(initMultiLine(&ml))*/

  if(addToMultiLine(&ml,s)){
     destroyMultiLine(&ml);
     return 1;
  }/*if(addToMultiLine(&ml,s))*/
  return 0;
}/*doAddString*/

void AddString (char* s)
{
   if(doAddString(s))
      MLPutSymbol(stdlink, "False");
   else
      MLPutSymbol(stdlink, "True");
}/*AddString*/


void Integrate(char* s) {
int i;

   if( (s!=NULL)&&( *s != '\0' ) )
      if(doAddString(s))
         return;
   
   fscan=newScannerMultiStr(&ml);

   if(fscan==NULL){
      PutErrorMessage("CIntegrate","Can't initialize scanner");
      return;
   }

   if(masterScan(fscan))
      return ;
   
#ifdef STATISTICS_OUT
   printf("Tetrads: total %d, left after optimization %d\n",allT,newT);
#endif

   if(wasCut){
      if(xCuts[0]>0.0){/*"Global cut*/
         *(fscan->ep=malloc(sizeof(FLOAT)))=xCuts[0];
         fscan->wasCut=1;
      }else{/*"Local" cut, 
         individually for each x:*/
         fscan->wasCut=0;
         for(i=1; (i<nCuts) && (i<fscan->nx);i++)
            if( xCuts[i]>0.0)
               ( fscan->wasCut)++;
         if(fscan->wasCut){
            fscan->ep=calloc(fscan->nx,sizeof(FLOAT));
            fscan->ep_i=calloc(fscan->nx,sizeof(int));
            fscan->ep[0]=xCuts[0];
            for(i=1; (i<nCuts) && (i<fscan->nx);i++)
               if(xCuts[i]>0.0){
                  fscan->ep_i[i]=1;
                  fscan->ep[i]=xCuts[i];
               }/*if(xCuts[i]>0.0)*/
         }/*if(fscan->wasCut)*/
      }/*else*/
   }/*if(wasCut)*/
   /*Now in fscan->rtTriade is the translated expression*/
   {/*Block*/
     double acc=1.e-5;
     int dim=fscan->nx-1, 
         ncall = numberOfCalls1,
         itmx = numberOfIterations1, 
         nprn = 1,
         doof=0;

         VEGAS(theIntegrand,&acc,&dim,&ncall,&itmx,&nprn,&doof);
         if(numberOfCalls2>0){
            ncall = numberOfCalls2;
            itmx = numberOfIterations2;
            VEGAS1(theIntegrand,&acc,&dim,&ncall,&itmx,&nprn,&doof);
         }/*if(numberOfCalls2>0)*/
   }/*Block*/

   destroyScanner(fscan);

   //RESULT.s1;
	MLPutFunction(stdlink,"List",2);
	MLPutReal(stdlink,RESULT.s1);
	MLPutReal(stdlink,RESULT.s2);
}/*Integrate*/
 
void ClearString (void)
{
  destroyMultiLine(&ml);
  MLPutSymbol(stdlink, "True");
}/*ClearString*/

int main(int argc, char *argv[])
{
	return MLMain(argc, argv);

}
