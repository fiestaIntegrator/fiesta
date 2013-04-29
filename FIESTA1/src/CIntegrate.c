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

#include "integrators.h"

#define CUT_MULTIPLIER 10.0

integrator_t g_integratorArray[MAX_INTEGRATOR];
char *g_integratorNamesArray[MAX_INTEGRATOR];
int g_currentIntegrator=0;
int g_topIntegrator=0;

scan_t *g_fscan=NULL;

#ifdef STATISTICS_OUT
extern long int allT;
extern long int newT;
#endif

void getCurrentIntegratorParameters(void)
{
   /*Att! Very dangerous! No overfolow checkup!*/
   char res[2048];
   if(getAllPars(g_currentIntegrator, res))
      MLPutString(stdlink,"Fault");
   else
      MLPutString(stdlink,res);
}/*getCurrentIntegratorParameters*/

void setCurrentIntegratorParameter(char *name, char *value)
{
   int ret=setPar(g_currentIntegrator,name,value);
   if(ret==0)
      MLPutSymbol(stdlink, "True");
   else
      MLPutSymbol(stdlink, "False");
}/*setCurrentIntegratorParameter*/

void setIntegrator(char *name)
{
int i;
   /*If te integrator array is already initialized, the function does nothing:*/
   initIntegrators();/*fixme: check the returned code!*/
   for(i=0;i<g_topIntegrator;i++)
      if( strcmp(name,g_integratorNamesArray[i])==0 ){
         g_currentIntegrator=i;
         MLPutSymbol(stdlink, "True");
         return;
      }
   MLPutSymbol(stdlink, "False");
}/*setIntegrator*/

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

static multiLine_t ml;

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
result_t results;

   if( (s!=NULL)&&( *s != '\0' ) )
      if(doAddString(s))
         return;
   g_fscan=newScannerMultiStr(&ml);

   if(g_fscan==NULL){
      PutErrorMessage("CIntegrate","Can't initialize scanner");
      return;
   }

   if(masterScan(g_fscan))
      return ;
   
#ifdef STATISTICS_OUT
   printf("Tetrads: total %d, left after optimization %d\n",allT,newT);
#endif

   if(wasCut){
      if(xCuts[0]>0.0){/*"Global cut*/
         *(g_fscan->ep=malloc(sizeof(FLOAT)))=xCuts[0];
         g_fscan->wasCut=1;
      }else{/*"Local" cut, 
         individually for each x:*/
         g_fscan->wasCut=0;
         for(i=1; (i<nCuts) && (i<g_fscan->nx);i++)
            if( xCuts[i]>0.0)
               ( g_fscan->wasCut)++;
         if(g_fscan->wasCut){
            g_fscan->ep=calloc(g_fscan->nx,sizeof(FLOAT));
            g_fscan->ep_i=calloc(g_fscan->nx,sizeof(int));
            g_fscan->ep[0]=xCuts[0];
            for(i=1; (i<nCuts) && (i<g_fscan->nx);i++)
               if(xCuts[i]>0.0){
                  g_fscan->ep_i[i]=1;
                  g_fscan->ep[i]=xCuts[i];
               }/*if(xCuts[i]>0.0)*/
         }/*if(g_fscan->wasCut)*/
      }/*else*/
   }/*if(wasCut)*/
   /*Now in g_fscan->rtTriade is the translated expression*/

   /*If te integrator array is already initialized, the function does nothing:*/
   initIntegrators();/*fixme: check the returned code!*/

   g_integratorArray[g_currentIntegrator](&results);
   destroyScanner(g_fscan);
	MLPutFunction(stdlink,"List",2);
	MLPutReal(stdlink,results.s1);
	MLPutReal(stdlink,results.s2);
}/*Integrate*/
 
void ClearString (void)
{

  destroyMultiLine(&ml);
  MLPutSymbol(stdlink, "True");
}/*ClearString*/

int main(int argc, char *argv[])
{
/*
volatile FILE *zzz=NULL;
while(zzz == NULL){
   sleep(1);
}
*/
	return MLMain(argc, argv);

}
