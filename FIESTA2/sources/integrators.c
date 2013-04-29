#include "integrators.h"

#include <string.h>

typedef struct paramsArrayStruct{
   char *name;
   char *fmt;
   void *addr;
}paramsArrayStruct_t;

#define FMT_F "%lf"
#define FMT_E "%lE"
#define FMT_I "%d"

static paramsArrayStruct_t *integratorParams[MAX_INTEGRATOR];

int registerIntegrator(integrator_t theIntegrator, char *name, paramsArrayStruct_t *params)
{
   if(g_topIntegrator >= MAX_INTEGRATOR)
      return -1;
   g_integratorNamesArray[g_topIntegrator]=malloc(strlen(name)+1);
   strcpy(g_integratorNamesArray[g_topIntegrator],name);
   integratorParams[g_topIntegrator]=params;
   g_integratorArray[g_topIntegrator++]=theIntegrator;
   return 0;
}/*registerIntegrator*/

/*********** CUBA: ***************************/
/*note, #include<cuba.h> is contained in "integrators.h"*/

static void IntegrandCuba(const int *ndim, const double xx[],
                      const int *ncomp, double ff[])
{
ff[0]=runExpr((double*)xx,g_fscan);
}/*IntegrandCuba*/

#define EPSREL 1e-5
#define EPSABS 1e-12
#define VERBOSE 0
#define LAST 4
#define MINEVAL 0
#define MAXEVAL 50000

#define NSTART 1000
#define NINCREASE 500

#define NNEW 1000
#define FLATNESS 25.

#define KEY1 0
#define KEY2 0
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25

#define NGIVEN 0

#define KEY 0

static int key=KEY;
static int key1=KEY1;
static int key2=KEY2;
static int key3=KEY3;

static FLOAT epsrel=EPSREL;
static FLOAT epsabs=EPSABS;
static int mineval=MINEVAL;
static int maxeval=MAXEVAL;
static int nstart=NSTART;
static int nincrease=NINCREASE;
static int nnew=NNEW;
static FLOAT flatness=FLATNESS; 
static int maxpass=MAXPASS;
static FLOAT border=BORDER;
static FLOAT maxchisq=MAXCHISQ;
static FLOAT mindeviation=MINDEVIATION;
/*
static int ngiven=NGIVEN;
*/
static paramsArrayStruct_t vegasCubaParameters[]={
   {"epsrel",FMT_E,&epsrel},
   {"epsabs",FMT_E,&epsabs},
   {"mineval",FMT_I,&mineval},
   {"maxeval",FMT_I,&maxeval},
   {"nstart",FMT_I,&nstart},
   {"nincrease",FMT_I,&nincrease},
   {NULL,NULL,NULL}
};

static int vegasCuba(result_t *result)
{
  int  neval, fail;
  double integral[1], error[1], prob[1];
  Vegas(g_fscan->nx-1, 1, IntegrandCuba,
    epsrel, epsabs, VERBOSE, mineval, maxeval,
    nstart, nincrease,
    &neval, &fail, integral, error, prob);
    result->s1=integral[0];
    result->s2=error[0];
    return (fail>=0);
}/*vegasCuba*/


static paramsArrayStruct_t suaveCubaParameters[]={
   {"epsrel",FMT_E,&epsrel},
   {"epsabs",FMT_E,&epsabs},
   {"mineval",FMT_I,&mineval},
   {"maxeval",FMT_I,&maxeval},
   {"nnew",FMT_I,&nnew},
   {"flatness",FMT_F,&flatness},
   {NULL,NULL,NULL}
};

static int suaveCuba(result_t *result)
{
  int  nregions, neval, fail;
  double integral[1], error[1], prob[1];

  Suave(g_fscan->nx-1,1, IntegrandCuba,
    epsrel, epsabs, VERBOSE | LAST, mineval, maxeval,
    nnew, flatness,
    &nregions, &neval, &fail, integral, error, prob);
  result->s1=integral[0];
  result->s2=error[0];
  return (fail>=0);
}/*suaveCuba*/


static paramsArrayStruct_t divonneCubaParameters[]={
   {"epsrel",FMT_E,&epsrel},
   {"epsabs",FMT_E,&epsabs},
   {"mineval",FMT_I,&mineval},
   {"maxeval",FMT_I,&maxeval},
   {"key1",FMT_I,&key1},
   {"key2",FMT_I,&key2},
   {"key3",FMT_I,&key3},
   {"maxpass",FMT_I,&maxpass},
   {"border",FMT_F,&border},
   {"maxchisq",FMT_F,&maxchisq},
   {"mindeviation",FMT_F,&mindeviation},
/*   {"ngiven",FMT_I,&ngiven},*/
   {NULL,NULL,NULL}
};


static int divonneCuba(result_t *result)
{
  int  nregions, neval, fail;
  double integral[1], error[1], prob[1];

  Divonne(g_fscan->nx-1,1, IntegrandCuba,
    epsrel, epsabs, VERBOSE, mineval, maxeval,
    key1, key2, key3, maxpass, border, maxchisq, mindeviation,
    NGIVEN, g_fscan->nx-1, NULL, 0, NULL,
    &nregions, &neval, &fail, integral, error, prob);
  result->s1=integral[0];
  result->s2=error[0];
  return (fail>=0);
}/*divonneCuba*/

static paramsArrayStruct_t cuhreCubaParameters[]={
   {"epsrel",FMT_E,&epsrel},
   {"epsabs",FMT_E,&epsabs},
   {"mineval",FMT_I,&mineval},
   {"maxeval",FMT_I,&maxeval},
   {"key",FMT_I,&key},
   {NULL,NULL,NULL}
};

static int cuhreCuba(result_t *result)
{
  int nregions, neval, fail;
  double integral[1], error[1], prob[1];
  Cuhre(g_fscan->nx-1, 1, IntegrandCuba,
    epsrel, epsabs, VERBOSE | LAST, mineval, maxeval,
    key,
    &nregions, &neval, &fail, integral, error, prob);
    result->s1=integral[0];
    result->s2=error[0];
    return (fail>=0);
}/*cuhreCuba*/

/*********** :CUBA ***************************/

int setpoints(int i1,int i2,int i11,int i22)
{
return 0;
}
#if 0
/********** The Fortran VEGAS: ***************/

static FLOAT accF=1.e-5;
static int numberOfCalls1 = 10000;
static int numberOfIterations1 = 5;
static int numberOfCalls2 = 100000;
static int numberOfIterations2 = 15;

static paramsArrayStruct_t vegasFParameters[]={
   {"acc",FMT_E,&accF},
   {"numberOfCalls1",FMT_I,&numberOfCalls1},
   {"numberOfIterations1",FMT_I,&numberOfIterations1},
   {"numberOfCalls2",FMT_I,&numberOfCalls2},
   {"numberOfIterations2",FMT_I,&numberOfIterations2},
   {NULL,NULL,NULL}
};

int setpoints(int i1,int i2,int i11,int i22)
{
   numberOfCalls1=i1;
   numberOfIterations1=i2;
   numberOfCalls2=i11;
   numberOfIterations2=i22;

   return 0;
}/*setpoints*/

static FLOAT theIntegrandF( FLOAT *x, FLOAT *wgt)
{
   return runExpr(x,g_fscan);
}/*theIntegrand*/


static int vegasF(result_t *result) 
{
int dim=g_fscan->nx-1, 
    ncall = numberOfCalls1,
    itmx = numberOfIterations1, 
    nprn = 0,
    doof=0;
    VEGAS(theIntegrandF,&accF,&dim,&ncall,&itmx,&nprn,&doof);
    if(numberOfCalls2>0){
       ncall = numberOfCalls2;
       itmx = numberOfIterations2;
       VEGAS1(theIntegrandF,&accF,&dim,&ncall,&itmx,&nprn,&doof);
    }/*if(numberOfCalls2>0)*/
    result->s1=RESULT.s1;
    result->s2=RESULT.s2;
    return 0;
}/*vegasF*/
/********** :The Fortran VEGAS ***************/
#endif

/*Returns 0 on success:*/
int setPar(int nIntegrator, char *name, char *val)
{
   paramsArrayStruct_t *itr;

   for(itr=integratorParams[nIntegrator];itr->name!=NULL;itr++)
      if(strcmp(itr->name,name)==0){
         if(sscanf(val,itr->fmt,itr->addr)!=1)
            return 2;
         return 0;
      }
   return 1;
}/*setPar*/

/*If the output of getAllPars is >= STR_WIDTH it breaks the line:*/
#define STR_WIDTH 50

/*Att! Very dangerous! No overfolow checkup!*/
int getAllPars(int nIntegrator, char *result)
{
   paramsArrayStruct_t *itr;
   char *str=result;
   char *lastStop=result;
   int wasComma=0;
   *str++=' ';*str++=' ';*str++=' ';
   *str++='{';
   for(itr=integratorParams[nIntegrator];itr->name!=NULL;itr++){
      sprintf(str,"{\"%s\",\"",itr->name);
      while((*str)!='\0')str++;
      if(strcmp(FMT_I,itr->fmt)==0){
         if(sprintf(str,itr->fmt,*((int*)itr->addr))<=0)
            return -2;
      }
      else if(strcmp(FMT_F,itr->fmt)==0){
         if(sprintf(str,itr->fmt,*((FLOAT*)itr->addr))<=0)
            return -2;
      }
      else if(strcmp(FMT_E,itr->fmt)==0){
         if(sprintf(str,itr->fmt,*((FLOAT*)itr->addr))<=0)
            return -2;
      }
      else /*Not implemented format*/
         return -3;

      while((*str)!='\0')str++;

      *str++='"';
      *str++='}';
      *str++=',';wasComma=1;

      if(str-lastStop >= STR_WIDTH){
        *str++='\n';*str++=' ';*str++=' ';*str++=' ';
        lastStop=str;
      }
   }/*for(itr=integratorParams[nIntegrator];itr->name!=NULL;itr++)*/
   if(wasComma)
      while(*str != ',')str--;
   *str++='}';*str='\0';
   return 0;
}/*getAllPars*/

static int justEvaluate(result_t *result) {
double xx[10];
float res;
/*xx[0]=1;
xx[1]=1;
xx[2]=0.000000000000;
xx[3]=1;
xx[4]=1;
xx[5]=1;   */ 

xx[0]=0.000;
xx[1]=1;
xx[2]=1;
xx[3]=1;
/*
xx[1]=1;
xx[2]=1; 
xx[3]=1;
xx[4]=1;
xx[5]=1;
xx[6]=1;
xx[7]=1;
xx[8]=1;*/

res=runExpr((double*)xx,g_fscan);
result->s1=res;
result->s2=1;
return 1;
};

static paramsArrayStruct_t justEvaluateParameters[]={
  {NULL,NULL,NULL}  
};



int initIntegrators(void)
{
   if(g_topIntegrator>0)
      return 0;
#if 0
   if(registerIntegrator(&vegasF,"vegasf",vegasFParameters))
      return -1;
#endif
   if(registerIntegrator(&vegasCuba,"vegasCuba",vegasCubaParameters))
      return -1;
   if(registerIntegrator(&suaveCuba,"suaveCuba",suaveCubaParameters))
      return -1;
   if(registerIntegrator(&divonneCuba,"divonneCuba",divonneCubaParameters))
      return -1;
   if(registerIntegrator(&cuhreCuba,"cuhreCuba",cuhreCubaParameters))
      return -1;
   if(registerIntegrator(&justEvaluate,"justEvaluate",justEvaluateParameters))
      return -1;
   return 0;
}/*initIntegrators*/
