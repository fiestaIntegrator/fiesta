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
In this file one can find a translator of an incoming expression into
triples (ctTriad) and quadruples (rtTriad, for historical reasons).

masterScan() perform initialisations and invokes the translator
simpleScan() for all functions and then for the expression.

simpleScan() is a recursive translator, it invokes scanTerm() for each
term and then it sums up all the terms in some "native" order.

For subexpressions scanTerm() invokes the next copy of simpleScan().
scanTerm() scans each factor to a "commuting diad":
Let's consider a term. Let's assume it has an implicit coefficient 1 in
front.  We consider every factor as a commuting diad of a form
"operation", "operand", e.g. 3*x[1]*x[3]/x[2] will be
considered as a set of the following commuting diads:
'*', 3
'*', x[1]
'*', x[3]
'/',x[2]
Now we can re-order all these diads without affecting the result.
We try to evaluate such a sequence in some native order, see the function
processDiads().

After the sequence of the triples (ctTriad) is ready, the function
buildRtTriad() translates the sequence to the array of quadruples rtTriad
of a form "operation", "operand 1", "operand 2" and "result". The
fields "operand 1" and "operand 2" are pointers to some memory
addresses.  The "result" field is either the result of the quadruple
evaluation (double type), or the address of the memory cell in which
the result must be stored. The point is that the size of the address
field on 64 bit platforms is the same as the size of double, and the
usage of indirect references is rather slow, so it seems to be
reasonable to use a new memory cell for each intermediate result.

However, on a modern CPU such an approach appears to be completely wrong.
Modern processors have a big write-back cache memory (megabytes), and
provided the same memory address is permanently updated it is never
written to RAM but stays in the cache. On the other hand, if we write intermediate
results every time to a new address, all this data is sooner or later
written to RAM (after the cache exhausted). That is why
for "small" jobs, when all intermediate results fit into the
processor cache, the interpreter stores results directly to the
"result" field of a quadruple, while for "large-scale" jobs the
translator creates some buffers and provides the quadruples with re-usable
addresses of elements of these buffers.

The algorithm switches from the "small" to "large" model when the
number of generated triples reaches some threshold which strongly
depends on the size of the processor L2 (L3) cache per active core and
it is a subject for tuning for every specific architecture. The
corresponding value is hard-coded, it can be changed in the file
scanner.h, the macro INDIRECT_ADDRESSING_THRESHOLD.

After quadruple array is ready it can be evaluated by the interpreter
runExpr(), file runline.c
*/

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <errno.h>
#include <string.h>

#ifndef NO_MATH
#include <mathlink.h>
#endif

#include "scanner.h"
#include "tries.h"
#include "constants.h"

#define O_B '['
#define C_B ']'

#define O_B_S "["
#define C_B_S "]"

/*mmap/munmap:*/
#include <unistd.h>
#include <sys/mman.h>


/*Attention!
The following three macros must be defined before including queue.h, otherwise
queues will use malloc() instead of mmap:*/
#define qFLAlloc(theSize) mmap(0,theSize,PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
#define qFLFree(theScratch,theLength) munmap(theScratch,theLength);
#define qFLAllocFail MAP_FAILED
#include "queue.h"


#ifndef NO_FILES
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#endif

#include <stdarg.h>

#define MAX_FLOAT_LENGTH 128

/*FL_something are bitwise flags:*/
#define FL_WITHOPTIMIZATION 1
#define FL_NEVERTLINE 2

#ifdef FOR_MATHEMATICA

void PutErrorMessage(char*,char*);

int stoprun=0;

void halt(int retval,char *fmt, ...)
{
va_list arg_ptr;
char buf[256];

   va_start (arg_ptr, fmt);
   vsnprintf(buf,255,fmt,arg_ptr);
   va_end (arg_ptr);
   if(stoprun==0)
      PutErrorMessage("CIntegrate",buf);
   stoprun=1;
   if(retval<0)
      exit(-retval);
}/*halt*/

#else

void halt(int retval,char *fmt, ...)
{
va_list arg_ptr;
   va_start (arg_ptr, fmt);
   vfprintf(stderr,fmt,arg_ptr);
   va_end (arg_ptr);
   exit(retval);

}/*halt*/

#endif

#if ARITHMETIC != NATIVE
static  void initIFArray(INTERNAL_FLOAT *f,SC_INT l)
{
   SC_INT i;
   for(i=0;i<l;i++)
      initIFvar(f++);
}/*initIFArray*/
#endif

/*
 The quickselect algorith, returns the value of the `k'-th largest
 element from the array `v' of dimension `dim'. `k'-th largest element
 means that in the array there are not more than `k' elements which are
 bigger than the found one, and at least `k' elements which are not
 smaller than the found one.
 Implementation by M.Tentyukov based on the original C.A.R. Hoare ideas.
*/
static SC_INLINE SC_INT selectK(SC_INT *v, SC_INT dim, SC_INT k) 
{ SC_INT
left=0,right=dim-1; SC_INT m,l,r;

   k--;
   while(left<right){
      /*0..left-1 are greatest elements,
        right+1..dim-1 -- small elements.*/
     
      /*some optimization:*/
      if(k==left){/*just one element is missing, find 
                    and return the maximum:*/
         m=v[left];
         for(r=right; r>left;r--)if(v[r]>m) m=v[r];
         return m;
      }
      /*partitioning: left elements are _greater_ (or 
        equal) than the pivot, right elements are _less_ (or
        equal) than the pivot. This is the OPPOSITE w.r.t.
        the usual quicksort!
       */
      l=left;r=right;m=(l+r)>>1;
      do{
         while(v[l]>v[m])l++;
         while(v[r]<v[m])r--;
restart:
         if(l<r){
            SC_INT tmp=v[l]; v[l]=v[r];v[r]=tmp;
            if(m==l){
               m=r;
               l++;
               while(v[l]>v[m])l++;
               goto restart;
            }
            if(m==r){
               m=l;
               r--;
               while(v[r]<v[m])r--;
               goto restart;
            }
            l++;r--;
         }/*if(l<r)*/
      }while(l<r);
      /*Here if l>r, then v[l]==v[r]*/
      /*Now from 0 to r all elements greater or equal v[m],
        from l to dim all elements less or equal v[m]*/
      if(k<r)/*Too many big elements, 
               drop the right partition and repeat:*/
         right=r-1;
      else if (k>l)/*Not enough big elements,try get 
                     more from the right partition:*/
         left=l+1;
      else/*found*/
         return v[r];
   }/*while(left<right)*/
   return v[right];
}/*selectK*/

/*Set in masterScan():*/
SC_INT totalInputLength=-1;

SC_INT estimateHashSize(scan_t *theScan){
return (totalInputLength / 4) + 5;
}/*estimateHashSize*/

#ifndef NO_FILES
scan_t *newScanner(char *fname)
{
   scan_t *res;
   int fd=0;
#if ARITHMETIC == MPFR
      mpfr_set_default_prec(g_default_precision);
#endif

   if(fname!=NULL){
      fd=open(fname,O_RDONLY);
      if(fd <0)
         return NULL;
   }/*if(fname!=NULL)*/
   /*else -- fd=0 by initialization*/

   res=malloc(sizeof(scan_t));
   res->multiLine=NULL;

   res->buf[0]='\0';
   res->theChar=res->buf;
   res->fd=fd;
   res->fline=NULL;
#ifdef MIXED_ARITHMETIC
   res->maxX=NULL;
   res->fNativeLine=NULL;
   res->nativeX=NULL;
   res->allMPvariables=NULL;
#endif

   res->x=NULL;
   res->f=NULL;
   res->ep=NULL;
   res->ep_i=NULL;
   res->wasCut=0;
   res->pstack=NULL;
   res->flags=0;
#ifdef WITH_OPTIMIZATION
   res->flags|=FL_WITHOPTIMIZATION;
#endif
   res->tline=NULL;
   res->condChain=NULL;
   res->allCondChain=NULL;   
   res->rtTriad.length=0;
   res->ctTriad.max=0;
   res->ctTriad.free=0;
   res->ctTriad.theScan=res;
   return res;
}/*newScanner*/
#endif

scan_t *newScannerFromStr(char *str)
{
   scan_t *res;
#if ARITHMETIC == MPFR
   mpfr_set_default_prec(g_default_precision);
#endif
   res=malloc(sizeof(scan_t));
   res->multiLine=NULL;

   res->theChar=str;
   res->fd=-1;
   res->fline=NULL;
#ifdef MIXED_ARITHMETIC
   res->fNativeLine=NULL;
   res->nativeX=NULL;
   res->maxX=NULL;
   res->allMPvariables=NULL;
#endif
   res->x=NULL;
   res->f=NULL;

   res->ep=NULL;
   res->ep_i=NULL;
   res->wasCut=0;
   res->pstack=NULL;
   res->flags=0;

#ifdef WITH_OPTIMIZATION
   res->flags|=FL_WITHOPTIMIZATION;
#endif
   res->tline=NULL;
   res->condChain=NULL;
   res->allCondChain=NULL;
   res->allocatedCondChains=NULL;
#ifdef MIXED_ARITHMETIC
   res->constStrings=NULL;
   res->newConstantsPool.full=0;
#endif
   res->rtTriad.length=0;
   res->ctTriad.max=0;
   res->ctTriad.free=0;
   res->ctTriad.theScan=res;

   return res;
}/*newScannerFromStr*/

scan_t *newScannerMultiStr(multiLine_t *multiLine)
{
   scan_t *res;
#if ARITHMETIC == MPFR
   mpfr_set_default_prec(g_default_precision);
#endif
   res=malloc(sizeof(scan_t));
   res->multiLine=multiLine;
   res->theChar=(char*)(multiLine->mline.buf[0]);
   res->fd=-1;
   res->fline=NULL;
#ifdef MIXED_ARITHMETIC
   res->fNativeLine=NULL;
   res->nativeX=NULL;
   res->maxX=NULL;
   res->allMPvariables=NULL;
#endif
   res->x=NULL;
   res->f=NULL;

   res->ep=NULL;
   res->ep_i=NULL;
   res->wasCut=0;
   res->pstack=NULL;
   res->flags=0;

#ifdef WITH_OPTIMIZATION
   res->flags|=FL_WITHOPTIMIZATION;
#endif
   res->tline=NULL;
   res->condChain=NULL;
   res->allCondChain=NULL;
   res->allocatedCondChains=NULL;
#ifdef MIXED_ARITHMETIC
   res->constStrings=NULL;
   res->newConstantsPool.full=0;
#endif
   res->rtTriad.length=0;
   res->ctTriad.max=0;
   res->ctTriad.free=0;
   res->ctTriad.theScan=res;

   return res;
}/*newScannerMultiStr*/

int initMultiLine(multiLine_t *multiLine)
{
char *c;
   initCollect(&multiLine->mline);
   
   c=(char*)mmap (0, MULTILINESIZE+1,
        PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(c == MAP_FAILED){
      halt(-10, "mmap failed\n");
      return 1;
   }
   c[0]='\0';
   pushCell(&multiLine->mline,(void*)c);

   multiLine->nbuf=multiLine->free=0;

   return 0;
}/*initMultiLine*/

int addToMultiLine(multiLine_t *multiLine, char *str)
{
char *m=(char*)(multiLine->mline.buf[multiLine->mline.fill-1])+multiLine->free;
   while(*str!='\0'){   
      for(;multiLine->free<MULTILINESIZE;multiLine->free++)
         if( (*m++=*str++)=='\0' ){
            str--;
            break;
         }
      if(multiLine->free==MULTILINESIZE){
         ((char*)(multiLine->mline.buf[multiLine->mline.fill-1]))[MULTILINESIZE]='\0';
         m=(char*)mmap (0, MULTILINESIZE+1, PROT_READ|PROT_WRITE,
             MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
         if(m== MAP_FAILED){
            halt(-10, "mmap failed\n");
            return 1;
         }/*if(m == MAP_FAILED)*/
         pushCell(&multiLine->mline,m);
         multiLine->free=0;
      }/*if(multiLine->free==MULTILINESIZE)*/
   }/*while(*str!='\0')*/
   return 0;
}/*addToMultiLine*/

void destroyMultiLine(multiLine_t *multiLine)
{
void *tmp;
   while(multiLine->mline.fill != 0 ){
      tmp=popCell(&multiLine->mline);
      if(tmp!=NULL)
         munmap(tmp,MULTILINESIZE+1);
   }/*while(multiLine->mline.fill != 0 )*/
   free(multiLine->mline.buf);
   multiLine->mline.buf=NULL;
   multiLine->free=0;
   multiLine->nbuf=0;
}/*destroyMultiLine*/


static SC_INLINE void ctTriadRealloc( ct_triad_t *t)
{
   SC_INT bytesForInt, newBytesForInt;
   SC_INT newMax=t->max * 2;
   char *newOperation;
   SC_INT *newOperand;

   newOperation = (char*)mmap (0, newMax,
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(newOperation == MAP_FAILED)
      halt(-10, "mmap failed\n");
   memcpy(newOperation,t->operation,t->max);
   munmap(t->operation,t->max);
   t->operation=newOperation;
   
   
   bytesForInt=t->max*sizeof(SC_INT);
   newBytesForInt=newMax*sizeof(SC_INT);

   newOperand=(SC_INT*)mmap (0, newBytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(newOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");
   memcpy(newOperand,t->triadType,bytesForInt);
   munmap(t->triadType,bytesForInt);
   t->triadType=newOperand;


   newOperand=(SC_INT*)mmap (0, newBytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(newOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");
   memcpy(newOperand,t->lastUsing,bytesForInt);
   munmap(t->lastUsing,bytesForInt);
   t->lastUsing=newOperand;

   newOperand=(SC_INT*)mmap (0, newBytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(newOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");
   memcpy(newOperand,t->refCounter,bytesForInt);
   munmap(t->refCounter,bytesForInt);
   t->refCounter=newOperand;

   newOperand=(SC_INT*)mmap (0, newBytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(newOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");
   memcpy(newOperand,t->firstOperand,bytesForInt);
   munmap(t->firstOperand,bytesForInt);
   t->firstOperand=newOperand;

   newOperand=(SC_INT*)mmap (0, newBytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(newOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");
   memcpy(newOperand,t->secondOperand,bytesForInt);
   munmap(t->secondOperand,bytesForInt);
   t->secondOperand=newOperand;
   t->max=newMax;
}/*ctTriadRealloc*/

static SC_INLINE void ctTriadInit(SC_INT hashSize, ct_triad_t *t)
{
   SC_INT bytesForInt;
   t->max=getpagesize();
   t->free=1;

   t->operation = (char*)mmap (0, t->max,
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(t->operation == MAP_FAILED)
      halt(-10, "mmap failed\n");

   bytesForInt=t->max*sizeof(SC_INT);

   t->triadType=(SC_INT*)mmap (0, bytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(t->triadType == MAP_FAILED)
      halt(-10, "mmap failed\n");

   t->lastUsing=(SC_INT*)mmap (0, bytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(t->lastUsing == MAP_FAILED)
      halt(-10, "mmap failed\n");

   t->refCounter=(SC_INT*)mmap (0, bytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(t->refCounter == MAP_FAILED)
      halt(-10, "mmap failed\n");

   t->firstOperand=(SC_INT*)mmap (0, bytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(t->firstOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");
   t->secondOperand=(SC_INT*)mmap (0, bytesForInt, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(t->secondOperand == MAP_FAILED)
      halt(-10, "mmap failed\n");

   hash_init(hashSize,&t->hashTable, t);
}/*ctTriadInit*/

static SC_INLINE void ctTriadFree( ct_triad_t *t)
{
   SC_INT bytesForInt=t->max * sizeof(SC_INT);
      munmap(t->operation,t->max);
      munmap(t->triadType,bytesForInt);
      munmap(t->lastUsing,bytesForInt);
      munmap(t->refCounter,bytesForInt);
      munmap(t->firstOperand,bytesForInt);
      munmap(t->secondOperand,bytesForInt);
      t->free=t->max=0;
      if(t->hashTable.theSizeInBytes > 0 )
         hash_destroy(&t->hashTable);
      t->theScan = NULL;
}/*ctTriadFree*/

static SC_INLINE void rtTriadInit(rt_triad_t *t, SC_INT theLength)
{
   t->operations=(char*)malloc(theLength);
   if(t->operations==NULL)
      halt(-10,"malloc failed\n");
   t->operands=(rt_triadaddr_t*)malloc(theLength*sizeof(rt_triadaddr_t));
   if(t->operands==NULL)
      halt(-10,"malloc failed\n");
#ifdef MIXED_ARITHMETIC
   t->nativeOperands=(rt_triadaddr_t*)malloc(theLength*sizeof(rt_triadaddr_t));
   if(t->nativeOperands==NULL)
      halt(-10,"malloc failed\n");
#endif
   t->length=theLength;
}/*rtTriadInit*/


#if ARITHMETIC != NATIVE
static SC_INLINE void rtTriadClear(rt_triad_t *t)
{
   char *rto=t->operations;
   char *rtoStop=t->operations+t->length;/*Note, in runline the 
       last iteration is out of the loop so it is shorter by 1.*/
   rt_triadaddr_t *rta=t->operands;
   for(;rto<rtoStop;rto++,rta++){
      switch(*rto){
         case ROP_JMP:
         case ROP_LEQ:
         case ROP_L:
         case ROP_EQ:
         case ROP_GEQ:
         case ROP_G:
         case ROP_NEQ:
            continue;
         case ROP_IPOW_F:
            clearIFvar(&(rta->aIF.result));
#if 0
            clearIFvar(&(rta->aIF.result),"aIF-");
#endif
            break;
         default:
         clearIFvar(&(rta->aF.result));
#if 0
         clearIFvar(&(rta->aF.result),"aF-");
#endif
      }/*switch(*rto)*/
   }/*for(;rto<rtoStop;rto++,rta++)*/
}/*rtTriadClear*/
#endif

static SC_INLINE void rtTriadFree(rt_triad_t *t)
{
  free(t->operands);
  t->operands=NULL;
#ifdef MIXED_ARITHMETIC
  free(t->nativeOperands);
  t->nativeOperands=NULL;
#endif
  free(t->operations);
  t->operations=NULL;
  t->length=0;
}/*rtTriadFree*/

/*if l<0, does not deallocate array itself:*/
static SC_INLINE void freeIFArray(INTERNAL_FLOAT *f,SC_INT l)
{
   if(l<0){
#if ARITHMETIC != NATIVE
      SC_INT i;
      l=-l;
      for(i=0;i<l;i++)
         clearIFvar(f+i);
#endif
   }
   else{
#if ARITHMETIC != NATIVE
      SC_INT i;
      for(i=0;i<l;i++)
         clearIFvar(f+i);
#endif
      free(f);
   }
}/*freeIFArray*/

void destroyScanner( scan_t *theScan)
{
   if(theScan==NULL)
      return;
#ifndef NO_FILES
   if(theScan->fd)/*Not a stdin*/
      close(theScan->fd);
#endif
   if(theScan->multiLine!=NULL)
      destroyMultiLine(theScan->multiLine);

   if(theScan->ctTriad.max!=0)
     ctTriadFree(&theScan->ctTriad);

   if(theScan->rtTriad.length!=0){
#if ARITHMETIC != NATIVE
      if( theScan->flags & FL_NEVERTLINE )
         rtTriadClear(&theScan->rtTriad);
#endif
      rtTriadFree(&theScan->rtTriad);
   }/*if(theScan->rtTriad.length!=0)*/
#if ARITHMETIC == MPFR
   mpfr_free_cache();
#endif
   if(theScan->f!=NULL)
      free(theScan->f);

   if(theScan->x!=NULL)
      freeIFArray(theScan->x,theScan->nx);
#ifdef MIXED_ARITHMETIC
   if(theScan->nativeX!=NULL)
      free(theScan->nativeX);
   if(theScan->maxX!=NULL)
      free(theScan->maxX);
#endif

   if(theScan->ep!=NULL)
      free(theScan->ep);
   if(theScan->ep_i!=NULL)
      free(theScan->ep_i);

   if(theScan->fline!=NULL){
      freeIFArray(theScan->fline->buf,theScan->fline->fill);
      free(theScan->fline);
   }/*if(theScan->fline!=NULL)*/

#ifdef MIXED_ARITHMETIC
   if(theScan->fNativeLine!=NULL){
      free(theScan->fNativeLine->buf);
      free(theScan->fNativeLine);
   }
   if(theScan->allMPvariables!=NULL){
      free(theScan->allMPvariables->buf);
      free(theScan->allMPvariables);
   }
#endif
   if(theScan->pstack!=NULL){
      free(theScan->pstack->buf);
      free(theScan->pstack);
   }

   if(theScan->tline!=NULL){
      void *tmp;
     /*
       structure of tline:
       [0] the 128-bytes aligned array of NREGS+2 INTERNAL_FLOAT's
       [1] the address of the allocated block with this array
       [2] the array of OFTUSED+1 INTERNAL_FLOAT's
       [3] and further -- arrays of TLINE_POOL INTERNAL_FLOAT's.
      */

      /*We can't free buf[0] since it is shifted, so
        only clear it, put the length negative:*/
      freeIFArray(theScan->tline->buf[0],-(NREGS+2));
      theScan->tline->buf[0]=NULL;
      /*And now free:*/
      free(theScan->tline->buf[1]);
      theScan->tline->buf[1]=NULL;
      /*Free often used pool:*/
      freeIFArray(theScan->tline->buf[2],OFTUSED+1);
      theScan->tline->buf[2]=NULL;
      /*Now only TLINE_POOL:*/
      while(  (tmp=popCell(theScan->tline))!=NULL )
         freeIFArray((INTERNAL_FLOAT*)tmp,TLINE_POOL);
      free(theScan->tline->buf);
      free(theScan->tline);
   }

   if(theScan->condChain!=NULL){
      free(theScan->condChain->buf);
      free(theScan->condChain);
   }
   if(theScan->allCondChain!=NULL){
      free(theScan->allCondChain->buf);
      free(theScan->allCondChain);
   }
   if(theScan->allocatedCondChains!=NULL){
      void *tmp;
      while(  (tmp=popCell(theScan->allocatedCondChains))!=NULL )
         free(tmp);
      free(theScan->allocatedCondChains->buf);
      free(theScan->allocatedCondChains);
   }
#ifdef MIXED_ARITHMETIC
   if(theScan->constStrings!=NULL){
      free(theScan->constStrings->buf);
      free(theScan->constStrings);
   }
   if(theScan->newConstantsPool.full > 0 ){
      destroyMPool(&theScan->newConstantsPool);
      theScan->newConstantsPool.full=0;
   }
#endif
   free(theScan);
   totalInputLength=-1; /*let it be processed normally by the next masterScan()*/
}/*destroyScanner*/

static SC_INLINE int ReadNextBlock(scan_t *theScan)
{
#ifndef NO_FILES
ssize_t r;
   if(theScan->fd<0){
      halt(14,"ReadNextBlock: unexpected end of input\n");
      return 1;
   }
   while( ( (r=read(theScan->fd,theScan->buf,INBUFSIZE))<0 ) && (errno == EINTR)  );
   if(r<1){
      halt(14,"ReadNextBlock: unexpected end of input\n");
      return 1;
   }
   theScan->buf[r]='\0';
   theScan->theChar=theScan->buf;
   return 0;
#endif
   if(theScan->multiLine!=NULL){
      munmap(theScan->multiLine->mline.buf[theScan->multiLine->nbuf],MULTILINESIZE+1);
      theScan->multiLine->mline.buf[theScan->multiLine->nbuf]=NULL;
      theScan->theChar=theScan->multiLine->mline.buf[++theScan->multiLine->nbuf];
      if(theScan->multiLine->nbuf>=theScan->multiLine->mline.fill){
         halt(14,"ReadNextBlock: unexpected end of input\n");
         return 1;
      }
   }/*if(theScan->multiLine!=NULL)*/
   return 0;
}/*ReadNextBlock*/

static SC_INLINE int skipUntilEndOfComment(scan_t *theScan)
{
   while(*(theScan->theChar)!='}'){
      (theScan->theChar)++;
      if( *(theScan->theChar) == '\0' ){
#ifdef FOR_MATHEMATICA
         if(theScan->multiLine == NULL)
            return 0;
#endif
         if(ReadNextBlock(theScan))
            return 1;
      }
   }
   /*Here *(theScan->theChar)=='}'*/
   (theScan->theChar)++;
   /*Note* here *(theScan->theChar) may be '\0', no problem,
     nextCharComplex() will invoke ReadNextBlock(). */
   return 0;
}/*skipUntilEndOfComment*/

/* The function returns the current character corresponding to *theChar,
   ignoring white spaces (<=' ' ), and sets the pointer theChar to the next
   position. The double \ (combination \\) is treated as one \ character.
   The single \ at the end of the line is ignored:
 */

static char nextCharComplex(scan_t *theScan)
{
   for(;;){
      switch(*(theScan->theChar)){
         case '\0':/*Each time reading a new buffer we add '\0' at the end*/
#ifdef FOR_MATHEMATICA
            if(theScan->multiLine == NULL)
               return '\0';
#endif

            if(ReadNextBlock(theScan))
               return '\0';
            continue;
         case '{':
            if(skipUntilEndOfComment(theScan))
               return '\0';
            continue;
         case '\\':/*Look at the next element: exceptions are '\n' and '\\'*/
            switch((theScan->theChar)[1]){
               case '\0':/*The block is expired*/
#ifdef FOR_MATHEMATICA
                  if(theScan->multiLine == NULL)
                     return '\0';
#endif
                  (theScan->theChar)++;
                  if(ReadNextBlock(theScan))
                     return '\0';
                  if(*(theScan->theChar)=='\n'){
                     (theScan->theChar)++;
                     continue;
                  }else{
                     if(*(theScan->theChar)=='\\')
                        (theScan->theChar)++;
                     return '\\';
                  }
               case '\n':
                  (theScan->theChar)+=2;
                  break;
              case '\\':
                  (theScan->theChar)+=2;
                  return '\\';
            }/*switch((theScan->theChar)[1])*/
            break;
      }/*switch(*(theScan->theChar))*/
      if(*(theScan->theChar)>' ')
         break;
      (theScan->theChar)++;
   }/*for(;;)*/
   return *(theScan->theChar)++;
}/*nextCharComplex*/

static SC_INLINE char nextChar(scan_t *theScan)
{

   if(*(theScan->theChar) > ' ' ){
      switch(*(theScan->theChar)){
         case '\\':
         case '{':
             return nextCharComplex(theScan);
         default:
             return *(theScan->theChar)++;
      }
   }
   while(*(theScan->theChar) == ' ')(theScan->theChar)++;
   if(*(theScan->theChar) > ' ' ){
      switch(*(theScan->theChar)){
         case '\\':
         case '{':
             return nextCharComplex(theScan);
         default:
             return *(theScan->theChar)++;
      }
   }
   return nextCharComplex(theScan);
}/*nextChar*/

static void parseError(char *msg)
{
   halt(16,"Parse error: %s\n",msg);
}

static SC_INLINE char parseInt(int *n,scan_t *theScan)
{
char c;
   for(;;)switch(c=nextChar(theScan)){
      case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
      case '8':case '9':
        *n=*n*10+c-'0';
        break;
      default:
        return c;
   }/*for(;;)switch(c=nextChar(theScan))*/
}/*parseInt*/

static SC_INLINE char parseLong(SC_INT *n,scan_t *theScan)
{
char c;
   for(;;)switch(c=nextChar(theScan)){
      case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
      case '8':case '9':
        *n=*n*10+c-'0';
        break;
      default:
        return c;
   }/*for(;;)switch(c=nextChar(theScan))*/
}/*parseLong*/


/* *maxL is the maximum length of the string:*/
static SC_INLINE char parseIntN(int *n,int *maxL,scan_t *theScan)
{
char c;
   for(;(*maxL)>0;)switch(c=nextChar(theScan)){
      case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
      case '8':case '9':
        *n=*n*10+c-'0';
        (*maxL)--;
        break;
      default:
        return c;
   }/*for(;;)switch(c=nextChar(theScan))*/
   return c;
}/*parseIntN*/

static SC_INLINE char parseInt2(int *e,int *n,scan_t *theScan)
{
char c;
   *e=1;
   for(;;)switch(c=nextChar(theScan)){
      case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
      case '8':case '9':
        (*e)*=10;
        *n=*n*10+c-'0';
        break;
      default:
        return c;
   }/*for(;;)switch(c=nextChar(theScan))*/
}/*parseInt2*/

static SC_INLINE char *parseConst(char *buf,int l,scan_t *theScan)
{
   for(;l>0;l--,buf++)switch(*buf=nextChar(theScan)){
      case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
      case '8':case '9':
        break;
      default:
        return buf;
   }/*for(;l>0;l--,buf++)switch(*buf=nextChar(theScan))*/
   parseError("Too long number");
   return NULL;
}/*parseConst*/

static SC_INLINE int simpleScan(int iniOp, scan_t *theScan);


static SC_INLINE char parseFloat(FLOAT *r,int n, scan_t *theScan)
{
int m=MAX_FLOAT_LENGTH;
char buf[MAX_FLOAT_LENGTH];
char *c, ret;
   sprintf(buf,"%d",n);
   for(c=buf;(*c)!='\0';c++)m--;
   
   for(;m>0;m--)switch((*c)=nextChar(theScan)){
      case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
      case '8':case '9':
         c++;
         break;
      default:
         m=-1;
   }/*for(;m>0;m--)switch((*c)=nextChar(theScan))*/
   if(m==0)
      return '\0';
   if((*c) == '.'){
      c++;
      (*c)=nextChar(theScan);
      if(  ((*c)=='e')||(*(c)=='E') ){
         c++;
         (*c)=nextChar(theScan);
      }
      if(  ((*c)=='-')||(*(c)=='+') ){
         c++;
         (*c)=nextChar(theScan);
      }
      if( ( (*c) < '0' ) || ( (*c) > '9' )  )
         return '\0';
      c++;
      for(m=MAX_FLOAT_LENGTH-(c-buf);m>0;m--)switch((*c)=nextChar(theScan)){
         case '0':case '1':case '2':case '3':case '4':case '5':case '6':case '7':
         case '8':case '9':
            c++;
            break;
         default:
            m=-1;         
      }
      if(m==0)
         return '\0';
   }/*if((*c) == '.')*/
   ret = *c;
   *c='\0';
   *r=atof(buf);
   return ret;
}/*parseFloat*/

#ifdef STATISTICS_OUT
long int allT=0;
long int newT=0;
#endif


/*returns  the triad number*/
static SC_INLINE SC_INT sCtTriad(int op,SC_INT triadType,SC_INT operand1,
                                  SC_INT operand2,ct_triad_t *t)
{
SC_INT n=t->free;
   if (n >= t->max)
      ctTriadRealloc(t);
   t->operation[n]=op;
   t->triadType[n]=triadType;
   t->firstOperand[n]=operand1;
   t->secondOperand[n]=operand2;
#ifdef STATISTICS_OUT
      allT++; newT++;
#endif
   t->free++;
   return n;
}/*sCtTriad*/

/*returns  the triad number*/
static SC_INLINE SC_INT tandsCtTriad(int op,SC_INT triadType,SC_INT operand1,SC_INT operand2,ct_triad_t *t)
{
SC_INT n=t->free;
SC_INT res;
   if (n >= t->max)
      ctTriadRealloc(t);
   t->operation[n]=op;
   t->firstOperand[n]=operand1;
   t->triadType[n]=triadType;
   t->secondOperand[n]=operand2;

#ifdef WITH_OPTIMIZATION
   if(t->theScan->flags & FL_WITHOPTIMIZATION) {
      res=hash_tands(n, &t->hashTable);
#ifdef STATISTICS_OUT
      allT++;
#endif
      if(n == res){/*A new triad was installed*/
         t->free++;
#ifdef STATISTICS_OUT
         newT++;
#endif
      }/*if(n == res)*/
   }/*if(t->theScan->flags & FL_WITHOPTIMIZATION)*/
   else{
      res=n;
      t->free++;
   }
#else
   res=n;
   t->free++;
#endif
   return res;
}/*tandsCtTriad*/


/*x[a]:*/
#define MK_X(a) (-(a))
/*f[a]:*/
#define MK_F(a) (-(theScan->nx +(a)))
/*fline->buf[a]:*/
#define MK_FLOAT(a) (-(theScan->nx +theScan->nf+1+(a)))

static SC_INLINE INTERNAL_FLOAT *getAddress(SC_INT op,scan_t *theScan)
{
   if(op > 0){/*Triad*/
      SC_INT res;
      if( theScan->flags & FL_NEVERTLINE)
         return &(theScan->rtTriad.operands[op-1].aF.result);
      res=theScan->ctTriad.lastUsing[op];
      if(res>0)/*tline*/
         return (INTERNAL_FLOAT*)(theScan->tline->buf[ res / TLINE_POOL]) + (res % TLINE_POOL);

      if(res<0)/*often used*/
         return (INTERNAL_FLOAT*)(theScan->tline->buf[2])-res-1;

      /*The result is contained by the triad:*/
      return &(theScan->rtTriad.operands[op-1].aF.result);
   }/*if(op > 0)*/
   if(op >= - theScan->nx) /*x*/
      return theScan->x -op;
   if(op >= -theScan->nx-theScan->nf)/*f*/
      return getAddress(theScan->f[-op-theScan->nx],theScan);
   /*Note, getAddress might be called recursively, its ok!
     E.g., the function body consists only from one 
     another function*/
   /*fline*/
   return theScan->fline->buf -op-(theScan->nx +theScan->nf+1);
}/*getAddress*/

#ifdef MIXED_ARITHMETIC
static SC_INLINE FLOAT *getNativeAddress(SC_INT op,scan_t *theScan)
{
   if(op > 0)/*Triad*/
      return &(theScan->rtTriad.nativeOperands[op-1].aFN.result);

   if(op >= - theScan->nx) /*x*/
      return theScan->nativeX -op;
   if(op >= -theScan->nx-theScan->nf)/*f*/
      return getNativeAddress(theScan->f[-op-theScan->nx],theScan);
   /*Note, getNativeAddress might be called recursively, its ok!
     E.g., the function body consists only from one 
     another function*/
   /*fline*/
   return theScan->fNativeLine->buf -op-(theScan->nx +theScan->nf+1);
}/*getNativeAddress*/
#endif

#define MKFP(OP,c) if(c==NULL)return R##OP##_F;return R##OP##_P

static SC_INLINE char mkROPlabels(char OPlabel,INTERNAL_FLOAT *c)
{
   switch(OPlabel){
      case OP_NEG:
         MKFP(OP_NEG,c);
      case OP_INV:
         MKFP(OP_INV,c);
      case OP_LOG:
         MKFP(OP_LOG,c);
      case OP_CPY:
         MKFP(OP_CPY,c);
      case OP_JMP:
         return ROP_JMP;
      case OP_LEQ:
         return ROP_LEQ;
      case OP_L:
         return ROP_L;
      case OP_EQ:
         return ROP_EQ;
      case OP_GEQ:
         return ROP_GEQ;
      case OP_G:
         return ROP_G;
      case OP_NEQ:
         return ROP_NEQ;
      case OP_CPY2:
         MKFP(OP_CPY2,c);
      case OP_IPOW:
         MKFP(OP_IPOW,c);
      case OP_IPOW2:
         MKFP(OP_IPOW2,c);
      case OP_IPOW3:
         MKFP(OP_IPOW3,c);
      case OP_IPOW4:
         MKFP(OP_IPOW4,c);
      case OP_IPOW5:
         MKFP(OP_IPOW5,c);
      case OP_IPOW6:
         MKFP(OP_IPOW6,c);
      case OP_POW:
         MKFP(OP_POW,c);
      case OP_MINUS:
         MKFP(OP_MINUS,c);
      case OP_DIV:
         MKFP(OP_DIV,c);
      case OP_PLUS:
         MKFP(OP_PLUS,c);
      case OP_MUL:
         MKFP(OP_MUL,c);
   }/*switch(OPlabel)*/
   halt(-10,"mkROPlabels %d: %p", OPlabel, c);
   return 0;
}/*mkROPlabels*/

/*Puts indices of triads with top OFTUSED theScan->ctTriad.refCounter
  to the array oftUsed:*/
static SC_INLINE void fillOftUsed(SC_INT *oftUsed, scan_t *theScan)
{
   SC_INT k,i;
   SC_INT *refCounter=theScan->ctTriad.refCounter;
   SC_INT *oftUsedStop=oftUsed+OFTUSED;

   memcpy(theScan->ctTriad.triadType,
          theScan->ctTriad.refCounter,
          theScan->ctTriad.free * sizeof(SC_INT));
   /*Get the largest k-th element:*/
   k=selectK(theScan->ctTriad.triadType+1, theScan->ctTriad.free-1, OFTUSED);
   for(i=1;(oftUsed<oftUsedStop)&&(i<theScan->ctTriad.free);i++){
      if(refCounter[i] >= k)
         *oftUsed++=i;
   }
   /*To be on the safe side:*/
  while(oftUsed<oftUsedStop)*oftUsed++=0;

}/*fillOftUsed*/

#ifdef MIXED_ARITHMETIC
static SC_INLINE void storeAddress(scan_t *theScan, INTERNAL_FLOAT *cell)
{
pushCell(theScan->allMPvariables,cell);
}/*storeAddress*/
#endif

static SC_INLINE void *allocatTLinePool(int *tlineInd)
{
   void *tmp=malloc(TLINE_POOL*sizeof(INTERNAL_FLOAT));
   if(tmp==NULL)
      halt(-10,"allocatTLinePool: malloc fails\n");
   /*will be deallocated in destroyScanner()*/
#if ARITHMETIC != NATIVE
   initIFArray(tmp,TLINE_POOL);
   /*will be cleared in destroyScanner()*/
#endif
   *tlineInd=0;
   return tmp;
}/*allocatTLinePool*/


#define BUILDTR_INCR_ALL_POINTERS myIndex++,\
    cto++,cta1++,cta2++,ctNotUsed1++,ctNotUsed2++,ctLastUsing++

/*Translator from triads to tetrads:*/
static SC_INLINE int buildRtTriad(scan_t *theScan)
{
char *cto=theScan->ctTriad.operation+1;
char *ctoStop=theScan->ctTriad.operation+theScan->ctTriad.free;
SC_INT *ctLastUsing=theScan->ctTriad.lastUsing+1;
SC_INT *ctNotUsed1=theScan->ctTriad.triadType+1;
SC_INT *ctNotUsed2=theScan->ctTriad.refCounter+1;
SC_INT *cta1=theScan->ctTriad.firstOperand+1;
SC_INT *cta2=theScan->ctTriad.secondOperand+1;
SC_INT oftUsed[OFTUSED];
SC_INT oftUsedInd=0;
INTERNAL_FLOAT *oftUsedPtr=NULL;/*allocated space will be deallocated 
              in the destructor of the tline collection*/
SC_INT myIndex=1;
int tlineInd=0;
void *regsQBuf[NREGS+3];
qFL_t regsQ=QFL0,/*Very often*/
      tlineQ=QFL0;

INTERNAL_FLOAT *ptr1,*ptr2,*ptr3;
#ifdef MIXED_ARITHMETIC
FLOAT *ptr1N,*ptr2N;
#endif
char *rto;
rt_triadaddr_t *rta;
#ifdef MIXED_ARITHMETIC
rt_triadaddr_t *rtaN;
#endif

   /*if the number of triads is < INDIRECT_ADDRESSING_THRESHOLD,
     the triad result is stored by value, if it is >=
     INDIRECT_ADDRESSING_THRESHOLD, the triad result is stored by 
     address:
   */

#if ARITHMETIC != NATIVE
   theScan->flags|=FL_NEVERTLINE;
#else
   if(theScan->ctTriad.free < INDIRECT_ADDRESSING_THRESHOLD)
      theScan->flags|=FL_NEVERTLINE;
#endif
   theScan->x=malloc(theScan->nx * sizeof(INTERNAL_FLOAT));

   if(theScan->x == NULL)
      halt(-10,"malloc failed\n");

#ifdef MIXED_ARITHMETIC
   theScan->nativeX=malloc(theScan->nx * sizeof(FLOAT));
   if(theScan->nativeX == NULL)
      halt(-10,"malloc failed\n");
#endif

#if ARITHMETIC != NATIVE
    initIFArray(theScan->x,theScan->nx);
#endif
#ifdef MIXED_ARITHMETIC
   {/*Block*/
      int i;
      for(i=0; i<theScan->nx; i++)
         pushCell(theScan->allMPvariables,theScan->x+i);
   }/*Block*/
#endif
   if( !(theScan->flags & FL_NEVERTLINE) ){

      /*Put OFTUSED most popular triads to the array oftUsed*/
      if(OFTUSED>0)
         fillOftUsed(oftUsed, theScan);

      /*We need not triadType anymore, reuse it as a "notUsed1" pointer:*/
      memset(ctNotUsed1,0,(theScan->ctTriad.free-1)*sizeof(SC_INT));
      /*We need not refCounter anymore, reuse it as a "notUsed2" pointer:*/
      memset(ctNotUsed2,0,(theScan->ctTriad.free-1)*sizeof(SC_INT));

      qFLInit(NULL,getpagesize()/sizeof(void*),&tlineQ,QFL_REALLOC_IF_FULL);

      theScan->tline=initCollect(NULL);

      /*Allocate special storage for frequently rewrittables:*/
      /*Attention! in the tline collection will be placed the part 
        aligned to 128! We can't free this space in this 
        routine (and put NULL to tline->buf[0]) since
        the array is needed during evaluation. Instead,
        we add the second head to the collection in the cell 1.
        In the tline destructor (see destroyScanner()) we just put
        tline->buf[0]=NULL.
      */
      /*use prt1 and ptr2 which are still unused:*/
      ptr2=ptr1=malloc((NREGS+2)*sizeof(INTERNAL_FLOAT)+128);
      if(ptr1==NULL)
         halt(-10,"NREGS: malloc fails\n");
      /*will be deallocated in destroyScanner()*/
      /*Align to most probable cacheline size (128):*/
      if( (long int)ptr1 % 128 )
         ptr1=(INTERNAL_FLOAT*)( ((~127l) & (long int)ptr1) + 128);
      /*Very dangerous! 
        First,this should be in the cell number 0!
        Second, we put here not the head of the allocated block,
        the head will be placed into the second cell, buf[1]:
      */
      pushCell( theScan->tline, ptr1);/*put the storage*/
      pushCell( theScan->tline, ptr2);/*put the head*/
#if ARITHMETIC != NATIVE
      initIFArray(ptr1,NREGS+2);
      /*Will be cleared in destroyScanner()*/
#endif
      /*This queue can't overfull:*/
      qFLInit(regsQBuf,NREGS+3,&regsQ,QFL_FAIL_IF_FULL);
      {/*Block*/
         long int i;
         for(i=NREGS-1;i>0;i--)
            qFLPushLifo(&regsQ, (void*)i);
      }/*Block*/
      /*Allocate special storage for often used:*/
      oftUsedPtr=malloc((OFTUSED+1)*sizeof(INTERNAL_FLOAT));
      if(oftUsedPtr==NULL)
         halt(-10,"OFTUSED: malloc fails\n");
      /*will be deallocated in destroyScanner()*/
#if ARITHMETIC != NATIVE
      initIFArray(oftUsedPtr,OFTUSED+1);
      /*Will be cleared in destroyScanner()*/
#endif
      /*Attention! This MUST be placed into cell number 2! If you change it,
        you have to change the index 2 of theScan->tline->buf[2] in
        the routine getAddress()!:*/
      pushCell( theScan->tline, oftUsedPtr);
      
      pushCell( theScan->tline, allocatTLinePool(&tlineInd));
   }/*if( !(theScan->flags & FL_NEVERTLINE) )*/

   rtTriadInit(&theScan->rtTriad, theScan->ctTriad.free-1);

   rto=theScan->rtTriad.operations;
   rta=theScan->rtTriad.operands;
#ifdef MIXED_ARITHMETIC
   rtaN=theScan->rtTriad.nativeOperands;
#endif
   for(;cto<ctoStop;BUILDTR_INCR_ALL_POINTERS){
      if(*cto == OP_JMP){
         *(rto++)=ROP_JMP;
         rta->aJ = *cta1;
         rta++;
#ifdef MIXED_ARITHMETIC
         rtaN->aJ = *cta1;
         rtaN++;
#endif
         continue;
      }/*if(*cto == OP_JMP)*/
      /*if the bit FL_NEVERTLINE is set in theScan->flags, use only aF addresses, i.e. the 
         result will be stored directly in triads.*/
      ptr3=NULL;/*ptr3 will indicate whether the result will be 
                  stored in a triad (if ptr3 == NULL), or the 
                  triad contains the address of the result, ptr3.*/
      if( !(theScan->flags & FL_NEVERTLINE) ){
         /* Free triads which are not used anymore.
            There are at most two of them, with indices ctNotUsed1
            and ctNotUsed2. The field lastUsing of these triads contains
            the information of where the reslt was stored:
            lastUsing > 0 corresponds to the tline element, tlineQ
            <=0 -- the result was stored directly or somewhere else.
         */
         {/*block*/
            SC_INT **nuInd,*notUsedIndBuf[3]={ctNotUsed1,ctNotUsed2,NULL};
            for(nuInd=notUsedIndBuf;*nuInd!=NULL;nuInd++){
               if( **nuInd != 0){
                  SC_INT ind=(theScan->ctTriad.lastUsing)[**nuInd];
                  if(ind>0){
                     /*sometimes size of SC_INT is less than the 
                       size of a pointer, so we introduce a new variable
                       tmp suitable for casting:*/
                     long int tmp=ind;
                     if(tmp<NREGS)/*the result was in regsQ:*/
                        qFLPushLifo(&regsQ, (void*)(tmp));
                     else/*the result was in tline:*/
                        qFLPushLifo(&tlineQ, (void*)(tmp));
                  }/*if(ind>0)*/
                  /*else -- the either 0, i.e., direct result,
                    or <0 -- often used*/
               }/*if( **nuInd != 0)*/
            }/*for(nuInd=notUsedIndBuf;nuInd!=NULL;nuInd++)*/
         }/*block*/
         /*Now set up the ctNotUsed field for the corresponding triad:*/
         if(*ctLastUsing != 0){
            /*Last using may be theScan->ctTriad.free - 1, then just ignore it.*/
            if(*ctLastUsing < theScan->ctTriad.free - 1){
               SC_INT *tmp=theScan->ctTriad.triadType+(*ctLastUsing + 1);
               if(*tmp!=0)/*first one is opccupied, use the second one:*/
                  tmp=theScan->ctTriad.refCounter+(*ctLastUsing + 1);
               *tmp=myIndex;
            }/*if(*ctLastUsing < theScan->ctTriad.free - 1)*/
         }/*if(*ctLastUsing != 0)*/
         /*Now determine where the result will be stored.
           We have to decide, whether it will be stored directly, of
           placed to tline.
           */
         if(*ctLastUsing >0){
            if((OFTUSED>0)&&(myIndex == oftUsed[oftUsedInd])){/*often occures*/
               ptr3=oftUsedPtr+oftUsedInd;
               oftUsedInd++;
               *ctLastUsing=-oftUsedInd;
            }/*if(myIndex == oftUsed[oftUsedInd])*/
            if(   (ptr3 == NULL) && ((NREGS>0)&&(*ctLastUsing - myIndex < NREGS-1) )  ){
                 /*short distance, use regs block*/
               long int tmp=(long int)qFLPop(&regsQ);
               /*this queue can't be empty! But sometimes this happen...*/
               if(tmp > 0){
                  ptr3=(INTERNAL_FLOAT*)(theScan->tline->buf[0]) + tmp;
                  *ctLastUsing=tmp;
               }
            }/*else if (*ctLastUsing - myIndex < NREGS )*/
            if(ptr3==NULL){
               void *tmp=qFLPop(&tlineQ);
                  if(tmp!=NULL){
                     ptr3=(INTERNAL_FLOAT*)(theScan->tline->buf[((SC_INT)tmp) / TLINE_POOL]);
                     ptr3+=(((long int)tmp) % TLINE_POOL);
                     *ctLastUsing= (long int)tmp;
                  }/*if(tmp!=NULL)*/
                  else{
                     if(tlineInd>=TLINE_POOL)
                        pushCell( theScan->tline, allocatTLinePool(&tlineInd));
                     ptr3=(INTERNAL_FLOAT*)(theScan->tline->buf[theScan->tline->fill-1])+tlineInd;
                     *ctLastUsing=((theScan->tline->fill-1)*TLINE_POOL+tlineInd);
                     tlineInd++;
                  }/*else*/
            }/*if(ptr3==NULL)*/
         }/*if(*ctLastUsing >0)*/
      }/*if( !(theScan->flags & FL_NEVERTLINE) )*/

      /*Now we know, where the result will be stored: in an address,
        pointed by ptr3, or (if ptr3=NULL) directly in triad. Some 
        operations (like comparison, jump) have no results, also 
        the last triad. ctLastUsing for them is 0.
       */
#if ARITHMETIC != NATIVE
      if( theScan->flags & FL_NEVERTLINE )switch(*cto){
         case OP_JMP:
         case OP_LEQ:
         case OP_L:
         case OP_EQ:
         case OP_GEQ:
         case OP_G:
         case OP_NEQ:
            break;
         case OP_IPOW:
            initIFvar(&(rta->aIF.result));
#ifdef MIXED_ARITHMETIC
            pushCell(theScan->allMPvariables,&(rta->aIF.result));
#endif
            break;
         default:
            initIFvar(&(rta->aF.result));
#ifdef MIXED_ARITHMETIC
            pushCell(theScan->allMPvariables,&(rta->aF.result));
#endif
      }/*if( theScan->flags & FL_NEVERTLINE )switch(*cto)*/
#endif/*if ARITHMETIC != NATIVE*/
      ptr1=getAddress(*cta1,theScan);
#ifdef MIXED_ARITHMETIC
      ptr1N=getNativeAddress(*cta1,theScan);
#endif
      /*OP_IPOW requires integer for a second argument:*/
      if(*cto == OP_IPOW){
        if(ptr3==NULL){
           rta->aIF.firstOperand=ptr1;
           rta->aIF.secondOperand=*cta2;
#ifdef MIXED_ARITHMETIC
           rtaN->aIFN.firstOperand=ptr1N;
           rtaN->aIFN.secondOperand=*cta2;
#endif
        }/*if(ptr3==NULL)*/
        else{
           rta->aIP.firstOperand=ptr1;
           rta->aIP.secondOperand=*cta2;
           rta->aIP.result=ptr3;
        }
      }/*if(*cto == OP_IPOW)*/
      else{
         ptr2=getAddress(*cta2,theScan);
         if(ptr3==NULL){
#ifdef MIXED_ARITHMETIC
            ptr2N=getNativeAddress(*cta2,theScan);
            rtaN->aFN.firstOperand=ptr1N;
            rtaN->aFN.secondOperand=ptr2N;
#endif            
            rta->aF.firstOperand=ptr1;
            rta->aF.secondOperand=ptr2;
         }/*if(ptr3==NULL)*/
         else{
            rta->aP.firstOperand=ptr1;
            rta->aP.secondOperand=ptr2;
            rta->aP.result=ptr3;
         }
      }/*if(*cto == OP_IPOW)...else*/
      rta++;
#ifdef MIXED_ARITHMETIC
      rtaN++;
#endif
      switch(*cto){
         case OP_CPY2:
            *(rto++)=mkROPlabels(OP_CPY2,ptr3);
            break;
         case OP_CPY:
            *(rto++)=mkROPlabels(OP_CPY,ptr3);
            break;
         case OP_LOG:
            *(rto++)=mkROPlabels(OP_LOG,ptr3);
            break;
         case OP_NEG:
            *(rto++)=mkROPlabels(OP_NEG,ptr3);
            break;
         case OP_INV:
            *(rto++)=mkROPlabels(OP_INV,ptr3);
            break;
         case OP_POW:
            *(rto++)=mkROPlabels(OP_POW,ptr3);
            break;
         case OP_IPOW:
            *(rto++)=mkROPlabels(OP_IPOW,ptr3);
            break;
         case OP_IPOW2:
            *(rto++)=mkROPlabels(OP_IPOW2,ptr3);
            break;
         case OP_IPOW3:
            *(rto++)=mkROPlabels(OP_IPOW3,ptr3);
            break;
         case OP_IPOW4:
            *(rto++)=mkROPlabels(OP_IPOW4,ptr3);
            break;
         case OP_IPOW5:
            *(rto++)=mkROPlabels(OP_IPOW5,ptr3);
            break;
         case OP_IPOW6:
            *(rto++)=mkROPlabels(OP_IPOW6,ptr3);
            break;
         case OP_MINUS:
            *(rto++)=mkROPlabels(OP_MINUS,ptr3);
            break;
         case OP_DIV:
            *(rto++)=mkROPlabels(OP_DIV,ptr3);
            break;
         case OP_PLUS:
            *(rto++)=mkROPlabels(OP_PLUS,ptr3);
            break;
         case OP_MUL:
            *(rto++)=mkROPlabels(OP_MUL,ptr3);
            break;
         case OP_LEQ:            
            *(rto++)=mkROPlabels(OP_LEQ,ptr3);
            break;
         case OP_L:
            *(rto++)=mkROPlabels(OP_L,ptr3);
            break;
         case OP_EQ:
            *(rto++)=mkROPlabels(OP_EQ,ptr3);
            break;
         case OP_GEQ:
            *(rto++)=mkROPlabels(OP_GEQ,ptr3);
            break;
         case OP_G:
            *(rto++)=mkROPlabels(OP_G,ptr3);
            break;
         case OP_NEQ:
            *(rto++)=mkROPlabels(OP_NEQ,ptr3);
      }/*switch(*cto)*/
   }/*for(;cto<ctoStop;BUILDTR_INCR_ALL_POINTERS)*/
   if( !(theScan->flags & FL_NEVERTLINE) ){
      qFLDestroy(&tlineQ);
   }/*if( !(theScan->flags & FL_NEVERTLINE) )*/

   return 0;
}/*buildRtTriad*/

/*Sets triad fields lastUsing and refCounter*/
static SC_INLINE void markTriad(SC_INT triad,SC_INT op1,SC_INT op2,scan_t *theScan)
{
   /*What about functions?:*/
   while( (op1< - theScan->nx)&&(op1 >= -theScan->nx-theScan->nf) )
      op1=theScan->f[-op1-theScan->nx];
   while( (op2< - theScan->nx)&&(op2 >= -theScan->nx-theScan->nf) )
      op2=theScan->f[-op2-theScan->nx];

   if(op1>0){
      if(theScan->ctTriad.lastUsing[op1]<triad){
         theScan->ctTriad.lastUsing[op1]=triad;
         (theScan->ctTriad.refCounter[op1])++;
      }
   }/*if(op1>0)*/
   if((op1!=op2)&&(op2>0)){
      if(theScan->ctTriad.lastUsing[op2]<triad){
         theScan->ctTriad.lastUsing[op2]=triad;
         (theScan->ctTriad.refCounter[op2])++;
      }
   }/*if((op1!=op2)&&(op2>0))*/
}/*markTriad*/

/*Builds a triad without hash and stack pollution.
  It returns a triad index, or 0 if fails:*/
static SC_INLINE SC_INT addOpNoHashNoStack(int op,SC_INT triadType,scan_t *theScan)
{
   SC_INT op1,op2;
   if( (op > 0 )&&(op<=OP_CPY) ){ /*unary*/
      op1=op2=popInt(theScan->pstack);
   }else{
      op2=popInt(theScan->pstack);
      op1=popInt(theScan->pstack);
   }
   if(op1==0||op2==0){
      parseError("Stack underflow");
      return 0;
   }
   /*reuse triadType:*/
   triadType=sCtTriad(op,triadType,op1,op2,&theScan->ctTriad);
   if(op!=OP_JMP) 
      markTriad(triadType,op1,op2,theScan);
   return triadType;
}/*addOpNoHashNoStack*/

static SC_INLINE int addCpyOp(scan_t *theScan)
{
SC_INT op=popInt(theScan->pstack),res;
   /*-1 is a non-optimisible triad type:*/
   res=sCtTriad(OP_CPY,-1,op,op,&theScan->ctTriad);
   addInt(theScan->pstack,res);
   markTriad(res,op,op,theScan);
   return 0;
}/*addCpyOp*/
      
static SC_INLINE SC_INT checkOp(int op, SC_INT op1, SC_INT op2, scan_t *theScan)
{
SC_INT triad;
SC_INT *s,*c;
   /*Form a new triad:*/
   triad=sCtTriad(op,0,op1,op2,&theScan->ctTriad);
   /*check if there is such a triad in the global scope:*/
   triad=hash_check(triad, &theScan->ctTriad.hashTable);
   theScan->ctTriad.free--;
#ifdef STATISTICS_OUT
   allT--; newT--;
#endif
   if(triad!=0)
      return triad;

   /*No is such a triad in the global scope, look it up in the chain:*/
   s=theScan->allCondChain->buf+theScan->allCondChain->fill;
   c=theScan->allCondChain->buf;
   for(;c<s;c++){
      /*Form a new triad:*/
      triad=sCtTriad(op,*c,op1,op2,&theScan->ctTriad);
      /*check if there is such a triad in the chain:*/
      triad=hash_check(triad, &theScan->ctTriad.hashTable);
      theScan->ctTriad.free--;
#ifdef STATISTICS_OUT
      allT--; newT--;
#endif
      if(triad!=0)/*found*/
         return triad;
   }/*for(;c<s;c++)*/
   /*not found*/
   return 0;
}/*checkOp*/

static SC_INLINE int addOp( int op, scan_t *theScan)
{
   SC_INT op1,op2,triad;
   if( (op > 0 )&&(op<=OP_CPY) ){ /*unary*/
      op1=op2=popInt(theScan->pstack);
   }else{
      op2=popInt(theScan->pstack);
      op1=popInt(theScan->pstack);
   }
   if(op1==0||op2==0){
      parseError("Stack underflow");
      return 1;
   }
#ifdef WITH_OPTIMIZATION
   if(theScan->flags & FL_WITHOPTIMIZATION) {
      if(theScan->allCondChain->fill > 0){/*"if" is active*/
            triad=checkOp(op,op1,op2,theScan);
            if(triad==0)/*not found -- install with the latest chain*/
              triad=tandsCtTriad(op,
                theScan->allCondChain->buf[theScan->allCondChain->fill-1],
                op1,op2,&theScan->ctTriad);            
            /*now triad is not 0*/
      }/*if(theScan->allCondChain->fill > 0)*/
      else
         triad=tandsCtTriad(op,0,op1,op2,&theScan->ctTriad);
   }/*if(theScan->flags & FL_WITHOPTIMIZATION)*/
   else
#endif /*#ifdef WITH_OPTIMIZATION*/
      triad=tandsCtTriad(op,0,op1,op2,&theScan->ctTriad);
   addInt(theScan->pstack,triad);
   markTriad(triad,op1,op2,theScan);
   return 0;
}/*addOp*/

/*Sets theScan->condChain and theScan->allCondChain:*/
static SC_INLINE int createCondChain(SC_INT cCond,scan_t *theScan)
{
   SC_INT theChain;
   addInt(theScan->condChain,cCond);
   /*The first element is used as a buffer length:*/
   theScan->condChain->buf[0]=theScan->condChain->fill;

   pushCell( theScan->allocatedCondChains, (void *)(theScan->condChain->buf) );

   theChain=hash_check(1-theScan->allocatedCondChains->fill, &theScan->ctTriad.hashTable);

   if(theChain==0){/*not found.*/ 
      /*Allocate, storage, copy from theScan->condChain and set into the hash:*/
      int l=(theScan->condChain->buf[0])*sizeof(SC_INT);
         theScan->allocatedCondChains->buf[theScan->allocatedCondChains->fill - 1]=
         (void*)memcpy(malloc(l),theScan->condChain->buf,l);
         theChain=hash_tands(1-theScan->allocatedCondChains->fill,&theScan->ctTriad.hashTable);
         if(theChain!=1-theScan->allocatedCondChains->fill){
            parseError("createCondChain: internal error");
            return 1;
         }/*if(theChain!=1-theScan->allocatedCondChains->fill*/
   }/*if(theChain==0)*/
   else/*Forget it:*/
      theScan->allocatedCondChains->fill--;
   addInt(theScan->allCondChain,-theChain);
   return 0;
}/*createCondChain*/

static SC_INLINE int scanIf(scan_t *theScan)
{
char c;
int op=0;
SC_INT jumpPosition, resultFirstBranch,cCond;
int memScanFlags=theScan->flags & FL_WITHOPTIMIZATION;

   if(nextChar(theScan)!='f'){
      parseError("'f' expected");
      return 1;
   }
   if(nextChar(theScan)!='('){
      parseError("'(' expected");
      return 1;
   }
   /*Condition is scanned withoput optimizations:*/
   theScan->flags &= ~FL_WITHOPTIMIZATION;
   if(simpleScan(0,theScan)<0)
      return 1;
   c=nextChar(theScan);
   switch (c){
      case '<':
         if( (c=nextChar(theScan)) == '=' ){
            op=OP_LEQ;
            c=nextChar(theScan);
         }else
            op=OP_L;
         break;
      case '=':
         op=OP_EQ;
         c=nextChar(theScan);
         break;
      case '>':
         if( (c=nextChar(theScan)) == '=' ){
            op=OP_GEQ;
            c=nextChar(theScan);
         }else
            op=OP_G;
         break;
      case '!':
         op=OP_NEQ;
         c=nextChar(theScan);
         break;
      default:
         parseError("Condition expected");
         return 1;
   }/*switch (c)*/
   if(c != '('){
      parseError("'(' expected");
      return 1;
   }/*if(c != '(')*/
   if(simpleScan(0,theScan)<0)
       return 1;

   /*Restore the optimization bit:*/
   if(memScanFlags)
      theScan->flags |= FL_WITHOPTIMIZATION;
   else
      theScan->flags &= ~FL_WITHOPTIMIZATION;


   /*Condition is scanned. Build the triad:*/
   cCond=addOpNoHashNoStack(op,-1,theScan);
   if( cCond == 0 ){
      parseError("scanIf: internal error");
      return 1;
   }
   /*Use jumpPosition as a tmp var:*/
   jumpPosition=hash_check(cCond, &theScan->ctTriad.hashTable);
   if(jumpPosition!=0)/*Already was*/
      cCond=jumpPosition;
   else /*The first "if", nstall it:*/
      if(cCond!=hash_tands(cCond,&theScan->ctTriad.hashTable)){
         parseError("scanIf: internal error");
         return 1;
      }

   /*Build the condition chain for the first branch:*/

   /*Now look the chain in the hash:*/
   
   if(createCondChain(cCond,theScan))
      return 1;
  /*now theScan->condChain and theScan->allCondChain are set properly*/

   /*Now we add fictitious jump instruction*/
   addInt(theScan->pstack,1);
   jumpPosition=addOpNoHashNoStack(OP_JMP,-1,theScan);
   if(jumpPosition == 0 )
      return 1;
   /*Scan the first branch:*/
   if(nextChar(theScan)!= '(' ){
      parseError("'(' expected");
      return 1;
   }

   if(simpleScan(0,theScan)<0)
       return 1;
   /*If no triads are generated, we have to a put copy:
     indeed, the second branch might need something changeable.
     The starting position we have in jumpPosition.*/
   if(theScan->ctTriad.free - 1 == jumpPosition)/*No triads*/
      addCpyOp(theScan);

   /*Now in the stack we have the address ot the first branch,
     and this is the triad address. Store it:*/
   resultFirstBranch=theScan->pstack->buf[theScan->pstack->fill-1];

   if(nextChar(theScan)!= '(' ){
      parseError("'(' expected");
      return 1;
   }

   /*Now we know the addres for the first jump. Note, there will be 
     one triad more in this branch (namely, "jump"), so ctTriad.free is ok.
     Install it:*/

   theScan->ctTriad.firstOperand[jumpPosition]=theScan->ctTriad.free-jumpPosition;
   theScan->ctTriad.secondOperand[jumpPosition]=theScan->ctTriad.free-jumpPosition;


   /*Add  the jump instruction*/
   addInt(theScan->pstack,1);/*fictitiuos address*/
   jumpPosition=addOpNoHashNoStack(OP_JMP,-1,theScan);
   if(jumpPosition == 0 )
      return 1;

   /*Process the second branch, the same condition, but negative */
   /* Remove the first branch:*/
   theScan->allCondChain->fill--;
   theScan->condChain->fill--;
   if(createCondChain(-cCond,theScan))
      return 1;
   /*now theScan->condChain and theScan->allCondChain are set properly*/

   /*Scan the second branch:*/   
   if(simpleScan(0,theScan)<0)
       return 1;

   /*Now we know the addres for the second jump. Note, there will be 
     one triad more in this branch (namely, "cpy2"), so ctTriad.free is ok.
     Install it:*/
   theScan->ctTriad.firstOperand[jumpPosition]=theScan->ctTriad.free-jumpPosition;
   theScan->ctTriad.secondOperand[jumpPosition]=theScan->ctTriad.free-jumpPosition;

   /*Now bild copy triad which copies the result of the scond branch to the 
     result of the first branch:*/
   jumpPosition= popInt(theScan->pstack);
   if(jumpPosition==0 ){
      parseError("scanIf: stack underflow");
      return 1;
   }

   sCtTriad(OP_CPY2,
         -1,
         jumpPosition,
         resultFirstBranch,
         &theScan->ctTriad);

   /*Note, resultFirstBranch is in the stack!*/

   /* Remove the second branch chain:*/
   theScan->allCondChain->fill--;
   theScan->condChain->fill--;

   return 0;
}/*scanIf*/


/*
   simpleScan() collects terms into an array of the following structure.
   This is the local variable of simpleScan(), and is reset by memset() in 
   the beginning of simpleScan()
*/
typedef struct termStruct_struct{
SC_INT addrr;

SC_INT nTriads;
SC_INT nMul;
SC_INT nDiv;
SC_INT nPow;
SC_INT nLog;
SC_INT nIf;
} termStruct_t;

typedef struct cDiad{
char operation;/* OP_MUL or OP_DIV */
/* 'x' -- just x;
   'X' -- "extra" x obtained after reduction of the power to 1;(not used now)
   'f' -- just f;
   'p' -- power;
   'l' -- log;
   'i' -- if;
   '\0' -- from nested scan;
   'c' -- coefficient: */
char theType;
SC_INT addr;
}cDiad_t;
/*The idea is to sort diads as mentioned, similar diads are sorted
according addr.*/

#define D_X 0
#define D_XX 1
#define D_F 2
#define D_P 3
#define D_L 4
#define D_I 5
#define D_0 6
#define D_C 7


#define MAX_DIADS_IN_CHAIN 32
#define MK_DLINK(D) if(chainRoots[D]==-1)chainRoots[D]=dCounter;\
                    else chainLinks[chainCurrentLink[D]]=dCounter;\
                    chainCurrentLink[D]=dCounter;\
                    chainLinks[dCounter]=-1


static SC_INLINE int processDiadChain(
                                       int *nDiadChains, 
                                       char *from,
                                       char *links,
                                       cDiad_t *diad,
                                       scan_t *theScan)
{
   if(*from < 0 )
      return 0;
   if(*nDiadChains==0){
      /*if nDiadChains==0, then the first operation may be '/'. This is very rare situation, 
        so we just add 1 to the stack, if so.*/
      if(diad[*from].operation == OP_DIV){
         addInt(theScan->pstack,MK_FLOAT(1));/*First operation is '/'*/
      }/*if(diad[*from].operation == OP_DIV)*/
      else{
         addInt(theScan->pstack,diad[*from].addr);
         *from=links[*from];
      }
      *nDiadChains=1;
   }/*if(*nDiadChains==0)*/
   while(*from >= 0){
      addInt(theScan->pstack,diad[*from].addr);
      if( addOp(diad[*from].operation,theScan) )
         return 1;
      *from=links[*from];
   }/*while(*from >= 0)*/
   return 0;
}/*processDiadChain*/

static SC_INLINE int processDiads(
                         int nDiadChains,
                         cDiad_t *diad,
                         char *chainRoots,
                         char *chainLinks,
                         scan_t *theScan)
{
   if(processDiadChain(&nDiadChains,chainRoots+D_X,chainLinks,diad,theScan))
       return 1;
   if(processDiadChain(&nDiadChains,chainRoots+D_C,chainLinks,diad,theScan))
       return 1;
   if(processDiadChain(&nDiadChains,chainRoots+D_F,chainLinks,diad,theScan))
       return 1;
   if(processDiadChain(&nDiadChains,chainRoots+D_P,chainLinks,diad,theScan))
       return 1;
   if(processDiadChain(&nDiadChains,chainRoots+D_L,chainLinks,diad,theScan))
       return 1;
   if(processDiadChain(&nDiadChains,chainRoots+D_I,chainLinks,diad,theScan))
       return 1;
   if(processDiadChain(&nDiadChains,chainRoots+D_0,chainLinks,diad,theScan))
       return 1;
   return 0;
}/*processDiads*/

#ifndef NO_MATH
static SC_INLINE char *evaluateString(char *str)
{
const char *theString;
   if( MLEvaluate(stdlink,str)==0 )
      return NULL;
   if( MLGetString(stdlink, &theString)==0 )
      return NULL;
   return (char *)theString;
}/*evaluateString*/
#endif

#ifdef MIXED_ARITHMETIC
static SC_INLINE int pushConst(char *str,int invAlso,scan_t *theScan)
{
SC_INT tmp=mpoolAdd(str,&(theScan->newConstantsPool),strlen(str)+1);
   if(tmp<0)
      return -1;
   addInt(theScan->constStrings, tmp);
   if(invAlso)
      addInt(theScan->constStrings,-1);
   return 0;
}/*pushConst*/
#endif

static SC_INLINE SC_INT installConstant(char *str,char *val,scan_t *theScan)
{
SC_INT cpos;
char *p;
FLOAT r=strtod(val,&p);

   if(p==NULL)
      return -1;
   if(*p!='\0')
      if(*p!='\n')
         return -1;

#ifdef MIXED_ARITHMETIC
   cpos=theScan->fNativeLine->fill;
   addNativeFloat(theScan->fNativeLine,r);
   addNativeFloat(theScan->fNativeLine,1.0l/r);
   if( pushConst(val,1,theScan)<0 )
      return -1;
#endif

#if ARITHMETIC == NATIVE
   cpos=theScan->fline->fill;
   addFloat(theScan->fline,r);
   addFloat(theScan->fline,1.0l/r);
#endif
   installTrie(str,&cpos,&(theScan->newConstantsTrie));
   return cpos;
}/*installConstant*/

static SC_INLINE SC_INT evaluateConstant(char *str,scan_t *theScan)
{
#ifndef NO_MATH
char *theAnswer;
#endif
SC_INT cpos=checkTrie(str,&(theScan->newConstantsTrie));
#ifndef NO_MATH
   if(cpos<0){
      char buf[128];
      sprintf(buf,"N[%s,4096]",str);
      theAnswer=evaluateString(buf);
      if(theAnswer==NULL){
         halt(16," Can't evaluate '%s'\n",buf);
         return -1;
      }
      cpos=installConstant(str,theAnswer,theScan);
      if(cpos<0)
         return -1;
      MLReleaseString(stdlink, theAnswer);
   }
#endif
   return cpos;
}/*evaluateConstant*/

/*Returns the position in fline or -1 on error. If val==NULL, assumes that val = str. If invAlso!=0,
  stores also the inverse value to the next after cpos cell:*/
static SC_INLINE SC_INT tAndSConst(char *str,char *val,int invAlso,scan_t *theScan)
{
SC_INT cpos;
char *p;
FLOAT r;

   cpos=checkTrie(str,&(theScan->newConstantsTrie));
   if(cpos>=0)
      return cpos; /*found*/

   if(val==NULL)
      val=str;

   r=strtod(val,&p);

   if(p==NULL)
      return -1;
   if(*p!='\0')
      if(*p!='\n')
         return -1;
#ifdef MIXED_ARITHMETIC
   cpos=theScan->fNativeLine->fill;
   addNativeFloat(theScan->fNativeLine,r);
   if(invAlso)
      addNativeFloat(theScan->fNativeLine,1.0l/r);
   if( pushConst(val,invAlso,theScan)<0 )
      return -1;
#endif
#if ARITHMETIC == NATIVE
   cpos=theScan->fline->fill;
   addFloat(theScan->fline,r);
   if(invAlso)
      addFloat(theScan->fline,1.0l/r);
#endif
   installTrie(str,&cpos,&(theScan->newConstantsTrie));
   return cpos;
}/*tAndSConst*/

static SC_INLINE char scanPolyGamma(scan_t *theScan)
{
   char c,*id=ID_CONST_POLYGAMMA;
   char polyGamma[128];
   int i,l=strlen(id);
   SC_INT cpos=-1;

   /*If ID_CONST_POLYGAMMA == "PolyGamma" then "Pol" is scanned, now scan "yGamma":*/
   for(i=3;i<l;i++){
      c=nextChar(theScan);
      if(c != id[i]){
         halt(16,"Parse error: %c expected instead of %c\n",id[i],c);
         return'\0';
      }
   }/*for(i=1;i<l;i++)*/
   i=l=0;
   if(nextChar(theScan)!=O_B){parseError("'"O_B_S"' expected");return '\0';}
   if(parseInt(&i,theScan)!=','){parseError("',' expected");return '\0';}
   if(parseInt(&l,theScan)!=C_B){parseError("'"C_B_S"' expected");return '\0';}
   /*now in i and l we have a weight*/

   sprintf(polyGamma,"PolyGamma[%d,%d]",i,l);

   if( (i==1)&&(l==1) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA1_1,1,theScan);
   else if( (i==2)&&(l==1) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA2_1,1,theScan);
   else if( (i==2)&&(l==2) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA2_2,1,theScan);
   else if( (i==2)&&(l==3) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA2_3,1,theScan);
   else if( (i==2)&&(l==4) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA2_4,1,theScan);
   else if( (i==3)&&(l==1) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA3_1,1,theScan);
   else if( (i==3)&&(l==2) )
      cpos=tAndSConst(polyGamma,STR_CONST_POLYGAMMA3_2,1,theScan);
   else
      cpos=evaluateConstant(polyGamma,theScan);

   if(cpos<0){
      halt(16,"Parse error: can't evaluate %s\n",polyGamma);
      return '\0';
   }

   addInt(theScan->pstack,MK_FLOAT(cpos));
   return nextChar(theScan);
}/*scanPolyGamma*/

/*
   simpleScan() collects terms into the 'termStruct' with the corresponding 
   prefix operations in 'operation':
*/
static SC_INLINE char scanTerm(int iniOp,
                               termStruct_t *termStruct, 
                               char *operation,
                               scan_t *theScan)
{
   cDiad_t diad[MAX_DIADS_IN_CHAIN];
   char chainRoots[8];/*indexes of diad for the first 'x',..., or -1 if absent*/
   char chainCurrentLink[8];/*the last link, or -1*/
   char chainLinks[MAX_DIADS_IN_CHAIN];/*index of the next element, or -1 if no.*/
   int dCounter=0;
   int nDiadChains=0;
   char c;
   int n;
      diad[0].operation=OP_MUL;
      memset(chainRoots,-1,8);
      memset(chainCurrentLink,-1,8);
      /*No reasons to initialize chainLinks*/
      for(;;){
         c=nextChar(theScan);
         if(c=='\\'){/*Mathematica60 bug: it sends "\\012" instead of '\n'*/
            if( (nextChar(theScan)!='0')||
                (nextChar(theScan)!='1')||
                (nextChar(theScan)!='2')
              ){ 
                  parseError("Unexpected '\\'");
                  return '\0';
               }
              c=nextChar(theScan);
         }/*if(c=='\\')*/

         /*The operand:*/
         diad[dCounter].theType='\0';
         switch(c){
            case '-':
              *operation=OP_MINUS;
              /*no break*/
            case '+':/*ignore leading '+':*/
              continue;
            case 'i':
              if(scanIf(theScan))
                 return '\0';
              termStruct->nIf++;
              diad[dCounter].theType='i';
              MK_DLINK(D_I);
              c=nextChar(theScan);
              break;
            case '(':
              if(simpleScan(0,theScan)<0)
                 return '\0';
              MK_DLINK(D_0);
              c=nextChar(theScan);
              break;
            case 'L':
              if(nextChar(theScan)!='o'){parseError("'o' expected");return '\0';}
              if(nextChar(theScan)!='g'){parseError("'g' expected");return '\0';}
            case 'l':
              if(nextChar(theScan)!=O_B){
                 parseError("'"O_B_S"' expected");
                 return '\0';
              }
              if(simpleScan(OP_LOG,theScan)<0)
                 return '\0';
              diad[dCounter].theType='l';
              MK_DLINK(D_L);
              termStruct->nLog++;
              c=nextChar(theScan);
              break;
            case 'f':
            case 'x':
              if(nextChar(theScan)!=O_B){
                 parseError("'"O_B_S"' expected");
                 return '\0';
              }
              n=0;
              if(parseInt(&n,theScan)!=C_B){
                 parseError("'"C_B_S"' expected");
                 return '\0';
              }
              /*Attention! The case x() will be treated as x[0]!*/
              if(c=='x'){
                 addInt(theScan->pstack,MK_X(n));
                 diad[dCounter].theType='x';
                 MK_DLINK(D_X);
              }else{
                 addInt(theScan->pstack,MK_F(n));
                 diad[dCounter].theType='f';
                 MK_DLINK(D_F);
              }
              c=nextChar(theScan);
              break;
            case 'P':
              /*Either P, or PolyGamma, or Power*/
              c=nextChar(theScan);
              if(c=='o'){
                 c=nextChar(theScan);
                 if(c=='w'){
                    if(nextChar(theScan)!='e'){parseError("'er[' expected");return '\0';}
                    if(nextChar(theScan)!='r'){parseError("'r[' expected");return '\0';}

                    case 'p':
                       if(nextChar(theScan)!=O_B){
                          parseError("'"O_B_S"' expected");
                          return '\0';
                       }

                       if(simpleScan(OP_POW,theScan)<0)
                          return '\0';

                       diad[dCounter].theType='p';
                       MK_DLINK(D_P);

                        termStruct->nPow++;
                        c=nextChar(theScan);
                        break;
                 }/*if(c=='w')*/
                 else if(c == ID_CONST_POLYGAMMA[2]){
                    /*"PolyGamma", see a macros *_CONST_POLYGAMMA* 
                       in the file constants.h*/
                    diad[dCounter].theType='c';
                    MK_DLINK(D_C);
                    c=scanPolyGamma(theScan);
                    if(c=='\0')
                       return '\0';
                    break;
                 }/*else if(c == ID_CONST_POLYGAMMA[3])*/
              }/*if(c=='o')*/
              else{/*"pi", see a macros *_CONST_PI in the file constants.h*/
                 SC_INT cpos=tAndSConst("Pi",STR_CONST_PI,1,theScan);
                 if(cpos<0){
                    parseError("Can't evaluate Pi");
                    return '\0';
                 }
                 diad[dCounter].theType='c';
                 MK_DLINK(D_C);
                 addInt(theScan->pstack,MK_FLOAT(cpos));
                 break;
              }/*if(c=='o')...else*/
              parseError("'P','Power' or 'PolyGamma' expected");
              return '\0';
            case '0':case '1':case '2':case '3':case '4':
            case '5':case '6':case '7':case '8':case '9':
              diad[dCounter].theType='c';
              MK_DLINK(D_C);
              {/*Block*/
                 char *ptrC,constBuf[CONST_BUF];
                 SC_INT cpos;
                 constBuf[0]=c;
                 ptrC=parseConst(constBuf+1,CONST_BUF-2,theScan);
                 if(ptrC == NULL)
                    return '\0';
                 c=*ptrC;
                 if(c=='.'){
                    ptrC=parseConst(ptrC+1,CONST_BUF-(ptrC-constBuf)-2,theScan);
                    if(ptrC == NULL)
                       return '\0';
                    c=*ptrC;
                 }
                 *ptrC='\0';/*Now in constBuf[] we have the string*/
                 /*Now if isFloat!=0 then it is a floating point const.
                   Also if l>8 we couldn't use integers for native arithmetics.*/
                 cpos=tAndSConst(constBuf,NULL,1,theScan);
                 if(cpos<0){
                    halt(16,"Parse error: can't convert %s to number\n",constBuf);
                    return '\0';
                 }/*if(cpos<0)*/
/*#ifndef MIXED_ARITHMETIC*/
/*Problems with increasing precision for inverse constants.
  Switch off for non-native arithmetics.*/

                 /*Replace division by multiplication:*/
                 if(diad[dCounter].operation == OP_DIV){
                    diad[dCounter].operation = OP_MUL;
                    addInt(theScan->pstack,MK_FLOAT(cpos+1));
                 }
                 else
/*#endif*/
                    addInt(theScan->pstack,MK_FLOAT(cpos));
              }/*Block*/
              break;
            default:
                 halt(16,"Parse error: unexpected %c\n",c);
                 return '\0';
         }/*switch(c)*/
         diad[dCounter].addr=popInt(theScan->pstack);
         dCounter++;
         if(dCounter >= MAX_DIADS_IN_CHAIN){
            if(processDiads(
                    nDiadChains,
                    diad,
                    chainRoots,
                    chainLinks,
                    theScan)
              )
               return '\0';
            nDiadChains++;
            dCounter=0;
         }

         /*The operation:*/
         switch(c){
            case '*':
               diad[dCounter].operation=OP_MUL;
               termStruct->nMul++;
               break;
            case '/':
               diad[dCounter].operation=OP_DIV;
               termStruct->nDiv++;
               break;
            case ',':
               if( iniOp!= OP_POW ){
                  parseError("Unexpected ','");
                  return '\0';
               }
               /*No break!*/
            case '+':
            case '-':
            case ')':
            case ']':
            case ';':

               if(processDiads(
                    nDiadChains,
                    diad,
                    chainRoots,
                    chainLinks,
                    theScan)
                 )
                  return '\0';
               nDiadChains++;
               dCounter=0;
               return c;
            default:
              halt(16,"Parse error: unexpected %c\n",c);
              return '\0';
         }/*switch(c)*/
      }/*for(;;)*/
}/*scanTerm*/

/*if xVar>0, this is the power of bare x[xVar]:*/
static SC_INLINE char processPower(SC_INT xVar,scan_t *theScan)
{
char c;
int n,neg=0;
   do{
      c=nextChar(theScan);
      switch(c){
            case '-':
               neg=-1;
               break;
            case '0':case '1':case '2':case '3':case '4':
            case '5':case '6':case '7':case '8':case '9':              
              n=c-'0';/*first char is scanned already!*/
              c=parseInt(&n,theScan);
              if(c == '.'){
                 int n2=0;
                 int e;
                 c=parseInt2(&e,&n2,theScan);
                 if(n2!=0){/*Fractional part present*/
                    FLOAT r=(FLOAT)n+((FLOAT)n2)/e;
                    if(neg<0)
                       r*=-1.0;

#ifdef MIXED_ARITHMETIC
                    addInt(theScan->pstack,MK_FLOAT(theScan->fNativeLine->fill));
                    addNativeFloat(theScan->fNativeLine,r);
                    addNativeFloat(theScan->fNativeLine,1.0l/r);
                    {/*Block*/
                       char buf[MAX_FLOAT_LENGTH];
                       sprintf(buf,"%d.%d",n,n2);
                       if(pushConst(buf,1,theScan))
                          return '\0';
                    }/*Block*/
#else
                    addInt(theScan->pstack,MK_FLOAT(theScan->fline->fill));
                    addFloat(theScan->fline,r);
#endif
                    if(addOp(OP_POW ,theScan))
                       return 0;
                    return c;
                 }/*if(n2!=0)*/
              }/*if(c == '.')*/
              /*Integer power!*/
              if(n==0){
                 /*a^0=1, so pop one element and put 1:*/
                 popInt(theScan->pstack);
                 addInt(theScan->pstack,MK_FLOAT(1));                 
              }else if (n>1){
                switch(n){
                   case OP_IPOW2:
                      if(addOp(OP_IPOW2 ,theScan))
                         return 0;
                      break;
                   case OP_IPOW3:
                      if(addOp(OP_IPOW3 ,theScan))
                         return 0;
                      break;
                   case OP_IPOW4:
                      if(addOp(OP_IPOW4 ,theScan))
                         return 0;
                      break;
                   case OP_IPOW5:
                      if(addOp(OP_IPOW5 ,theScan))
                         return 0;
                      break;
                   case OP_IPOW6:
                      if(addOp(OP_IPOW6 ,theScan))
                         return 0;
                      break;
                   default: 
                      addInt(theScan->pstack,n);
                      if(addOp(OP_IPOW ,theScan))
                         return 0;
                }/*switch(n)*/
              }/*if (n>1)*/
              if(neg<0){
                 if(addOp(OP_INV ,theScan))
                     return 0;
#ifdef MIXED_ARITHMETIC
                 if(xVar>0){
                    if (theScan->maxX[xVar] == 0)
                        theScan->maxX[0]++;
                    if (theScan->maxX[xVar] < n)
                        theScan->maxX[xVar]= n;
                 }/*if(xVar>0)*/
#endif
              }/*if(neg<0)*/
              return c;
      }/*switch(c)*/
   }while(neg<0);
   return c;
}/*processPower*/


static SC_INLINE int cmpTerms(termStruct_t *term1,termStruct_t *term2)
{
   return (
      (term1->nTriads == term2->nTriads)&&
      (term1->nMul == term2->nMul)&&
      (term1->nDiv == term2->nDiv)&&
      (term1->nPow == term2->nPow)&&
      (term1->nLog == term2->nLog)&&
      (term1->nIf == term2->nIf)
   ); 
}/*cmpTerms*/

static SC_INLINE void couplePairs(size_t memStack,scan_t *theScan)
{
   if(theScan->pstack->fill == memStack+1){
      addOp(OP_PLUS,theScan);
      return;
   }
   if(theScan->pstack->fill == memStack+2){
      addOp(OP_PLUS,theScan);
      addOp(OP_PLUS,theScan);
      return;
   }
   /**/
   {
      SC_INT mem=popInt(theScan->pstack);
         addOp(OP_PLUS,theScan);
         couplePairs(memStack,theScan);  
         addInt(theScan->pstack,mem);
         addOp(OP_PLUS,theScan);
   }
}/*couplePairs*/

/*
  Initial state: one term in the stack.Final state: one term in the stack:
 */
static SC_INLINE int processSameTerms(int i,
                                      termStruct_t *term, 
                                      char *operations, 
                                      int termCounter,
                                      scan_t *theScan)
{
   size_t memStack=theScan->pstack->fill;
   char lastOp='\0';
   SC_INT lastAddrr=0;
   termStruct_t *cterm=term+i;
   char *coperation=operations+i;
   for(;;){
      for(i++;i<termCounter;i++)if( (term[i].addrr!=0) &&(cmpTerms(cterm,term+i)) )
            break;
      if(i==termCounter){/*non-paired term*/
         lastAddrr=cterm->addrr;
         cterm->addrr=0;
         lastOp=*coperation;
         break;
      }/*if(i==termCounter)*/
      /*Now cterm and term+i is a pair.*/
      if(*coperation==OP_MINUS){
         if(operations[i]==OP_MINUS)
            operations[i]=OP_PLUS;
         else
            operations[i]=OP_MINUS;
         addInt(theScan->pstack,cterm->addrr);
         addInt(theScan->pstack,term[i].addrr);
         addOp(operations[i],theScan);
         addOp(*coperation,theScan);
      }else{/* *coperation==OP_PLUS */
         addInt(theScan->pstack,cterm->addrr);
         addInt(theScan->pstack,term[i].addrr);
         addOp(operations[i],theScan);/*now the result is in the stack*/
         /* *coperation is OP_PLUS, we know!*/
      }/*if(coperation==OP_MINUS)*/
      cterm->addrr=0;
      term[i].addrr=0;

      /*Try to find next element:*/
      for(i++;i<termCounter;i++)if( (term[i].addrr!=0) &&(cmpTerms(cterm,term+i)) )
            break;
      if(i==termCounter)/*not found*/
         break;
      /*Process the rest of the chain:*/
      cterm=term+i;
      coperation=operations+i;
   }/*for(;;)*/
   /*Now in the stack may be several operand all with positve coefficients.*/
   if(memStack != theScan->pstack->fill)
      couplePairs(memStack,theScan);
   if(lastOp!='\0'){
      addInt(theScan->pstack,lastAddrr);
      addOp(lastOp,theScan);
   }
   return 0;
}/*processSameTerms*/

/*
  Initial state: no terms in the stack.The first operation must be OP_PLUS.
  Final state: one term in the stack:
 */
static SC_INLINE int addSameTerms(int i,
                                      termStruct_t *term,
                                      char *operations,
                                      int termCounter,
                                      scan_t *theScan)
{
termStruct_t *cterm=term+i;
   for(i++;i<termCounter;i++)
      if(cmpTerms(cterm,term+i))
         break;
   if(i==termCounter){/*Only one such a term, just add it to the stack:*/
      addInt(theScan->pstack,cterm->addrr);
      cterm->addrr=0;
      return 0;
   }/*if(i==termCounter)*/
   /*At least one pair*/
   addInt(theScan->pstack,cterm->addrr);
   addInt(theScan->pstack,term[i].addrr);
   addOp(operations[i],theScan);/*No reason for the error checking*/
   /*Need not term->addrr anymore, mark as used:*/
   cterm->addrr=0;
   term[i].addrr=0;
   for(i++;i<termCounter;i++)
      if(cmpTerms(cterm,term+i))
         break;
   if(i==termCounter)/*There are no more such terms*/
      return 0;
   /*There are other terms in a chain -- normal process:*/
   return processSameTerms(i,term,operations,termCounter,theScan);
}/*addSameTerms*/

static SC_INLINE int processTermsChain(termStruct_t *term, 
                                       char *operations,
                                       int termCounter,
                                       size_t memStack,
                                       scan_t *theScan)
{
   SC_INT i;
   switch(termCounter){
      case 0:
         return 0;
      case 1:
         if(memStack!=theScan->pstack->fill){/*continuation!*/
            /*The stack contains somethnig from the previous buffer,
              add the current term:*/
            addInt(theScan->pstack,term->addrr);
            addOp(*operations,theScan);
         }else{/*just put to the stack*/
            addInt(theScan->pstack,term->addrr);
            if(*operations == OP_MINUS)
               addOp(OP_NEG,theScan);
         }
         return 0;
      case 2:
         if(memStack==theScan->pstack->fill){
            /*The only result from the expression*/
            addInt(theScan->pstack,term->addrr);
            if(*operations == OP_MINUS)
               addOp(OP_NEG,theScan);         
            addInt(theScan->pstack,term[1].addrr);
            addOp(operations[1],theScan);/*No reason for the error checking*/
         }else{
            /*The stack contains somethnig from the previous buffer,
              perform binary operation and add the result to the stack:*/
            addInt(theScan->pstack,term->addrr);
            if(*operations==OP_MINUS){
               if(operations[1]==OP_MINUS)
                  operations[1]=OP_PLUS;
               else
                  operations[1]=OP_MINUS;
            }/*if(coperation==OP_MINUS)*/
            addInt(theScan->pstack,term[1].addrr);
            addOp(operations[1],theScan);/*No reason for the error checking*/
            addOp(*operations,theScan);/*No reason for the error checking*/
         }
         return 0;
      default:
         break;
   }/*switch(termCounter)*/

   if(memStack==theScan->pstack->fill){
      /*Stack is empty, there may be a problem with negative terms,
        so we try to find a group with leading positive term:*/
      for(i=0; i<termCounter;i++)
         if( (operations[i]==OP_PLUS)&&(term[i].nIf==0) )
            break;
      if(i<termCounter)/*there is at least one term with a positive coefficient*/
         addSameTerms(i,term, operations,termCounter,theScan);
      else /*all terms are negative, put 0 into the stack*/
         addInt(theScan->pstack,MK_FLOAT(0));
   }/*if(memStack==theScan->pstack->fill)*/

   /*Now we have something in the stack.*/

   for(i=0; i<termCounter;i++)if( (term[i].addrr!=0)&&(term[i].nIf==0) )
      processSameTerms(i,term, operations,termCounter,theScan);

   for(i=0; i<termCounter;i++)if(term[i].addrr!=0)
      processSameTerms(i,term, operations,termCounter,theScan);

   return 0;
}/*processTermsChain*/

/*
   simpleScan collects terms  into the array term, with the corresponding prefix operations 
   in the array operations. Each time the term chain become longer than MAX_TERMS_IN_CHAIN,
   it builds triads.
   Returns number of processed terms but not bigger than MAX_TERMS_IN_CHAIN.
   The point is that we need this value only to distinguish single term and many terms
   and try to avoid integer overflow.
*/
#define MAX_TERMS_IN_CHAIN 1024
static SC_INLINE int simpleScan(int iniOp, scan_t *theScan)
{
char c='\0';
termStruct_t term[MAX_TERMS_IN_CHAIN];
char operations[MAX_TERMS_IN_CHAIN];
int termCounter=0,wasChainProcessed=0;
size_t memStack=theScan->pstack->fill;
      operations[0]=OP_PLUS;
      while(c!='q'){
         memset(term+termCounter,0,sizeof(termStruct_t));
         c=scanTerm(iniOp,term+termCounter,operations+termCounter,theScan);
         if(c=='\0')
            return -1;
         term[termCounter].addrr=popInt(theScan->pstack);
         termCounter++;

         if(termCounter>=MAX_TERMS_IN_CHAIN){
            processTermsChain(term,operations,termCounter,memStack,theScan);
            termCounter=0;wasChainProcessed=1;
         }
         switch(c){
            case ',':/*p(... , */
               processTermsChain(term,operations,termCounter,memStack,theScan);
#ifdef MIXED_ARITHMETIC
               if(wasChainProcessed || termCounter!=1 )
                  c=processPower(0,theScan);
               else{
                  SC_INT xVar=-term[0].addrr;
                  if( (xVar < 1)|| (xVar > theScan->nx) )
                     xVar=0;
                  c=processPower(xVar,theScan);
               }
#else
               c=processPower(0,theScan);
#endif
               if(c!=C_B){
                  if(c!='\0')
                     parseError("'"C_B_S"' expected");
                  return -1;
               }
               return (wasChainProcessed)?MAX_TERMS_IN_CHAIN:termCounter;
            case '-':
              operations[termCounter]=OP_MINUS;
              break;
            case '+':
              operations[termCounter]=OP_PLUS;
              break;
            case ')':
            case ']':
            case ';':      
               c='q';/*Quit from the while()*/
              break;
            default:
              halt(16,"Parse error: unexpected %c\n",c);
              return -1;
         }/*switch(c)*/
      }/*while(c!='q')*/

      processTermsChain(term,operations,termCounter,memStack,theScan);
      if(iniOp == OP_LOG )
         if(addOp(OP_LOG,theScan))
            return -1;
      return (wasChainProcessed)?MAX_TERMS_IN_CHAIN:termCounter;
}/*simpleScan*/

#ifdef MIXED_ARITHMETIC
static SC_INLINE int intlog2(FLOAT x){return ceil(log(x)*1.44269504088896l);}

static SC_INLINE int getPrecision(scan_t *theScan)
{
int i;
FLOAT prod=1.0;
   if(theScan->maxX[0]>0){
      for(i=1;i<theScan->nx; i++){
         if(theScan->maxX[i]>0) prod*=ipow(g_mpsmallX,theScan->maxX[i]);
      }
   }
   if(prod>g_mpthreshold)
      prod=g_mpthreshold/2.0;

   if(prod<ABS_MONOM_MIN)
      prod=ABS_MONOM_MIN;/*To avoid overflow*/

   if(g_mpminChanged==0)
      g_mpmin=prod;/*else -- user defined*/
   return intlog2(1.0/prod)+g_mpPrecisionShift;
}/*getPrecision*/
#endif

int masterScan(scan_t *theScan)
{
int nx=0,nf=0,i;
char c=nextChar(theScan);
/*The scanner is initialised, the pline is not yet*/
   switch(c){
            case '0':case '1':case '2':case '3':case '4':
            case '5':case '6':case '7':case '8':case '9':
              /*Parse the total length of the input:*/
              if(theScan->multiLine ==NULL){/*the actual value of the length*/
                 SC_INT newLen=c-'0';/*first char is scanned already!*/
                 if(parseLong(&newLen,theScan)==';'){
                    if(totalInputLength<0)
                       totalInputLength=newLen;
                    break;
                 }
              }
              else{/*get the total length of the input from scanned multiLine,
                     ingore the value from the stream:*/
                 SC_INT l=c-'0';/*first char is scanned already!*/
                 if(parseLong(&l,theScan)==';'){/*Just drop out the value of l!*/
                       totalInputLength=theScan->multiLine->free;
                       if(theScan->multiLine->mline.fill>0)
                          totalInputLength+=
                          (theScan->multiLine->mline.fill - 1)*MULTILINESIZE;
                    break;
                 }
              }
            default:
                 halt(16,"Parse error: the input length with ; expected instead of %c\n",c);
                 return 1;
   }/*switch(c)*/
   
   switch(c=nextChar(theScan)){
            case '\\':/*Mathematica60 bug: it sends "\\012" instead of '\n'*/
               if( (nextChar(theScan)!='0')||
                   (nextChar(theScan)!='1')||
                   (nextChar(theScan)!='2')
                 ){
                    parseError("Unexpected '\\'");
                    return 1;
               }
                 c=nextChar(theScan);
                 if( (c<'0')||(c>'9') ){
                    parseError("Number expected");
                    return 1;
                 }
            case '0':case '1':case '2':case '3':case '4':
            case '5':case '6':case '7':case '8':case '9':
              nx=c-'0';/*first char is scanned already!*/
              if(parseInt(&nx,theScan)==';')
                 break;
            default:                 
                 halt(16,"Parse error: number of x with ; expected instead of %c\n",c);
                 return 1;
   }/*switch(c)*/
   nx++;/*f counted from 1, not 0!*/
   switch(c=nextChar(theScan)){
            case '\\':/*Mathematica60 bug: it sends "\\012" instead of '\n'*/
               if( (nextChar(theScan)!='0')||
                   (nextChar(theScan)!='1')||
                   (nextChar(theScan)!='2')
                 ){
                    parseError("Unexpected '\\'");
                    return 1;
               }
                 c=nextChar(theScan);
                 if( (c<'0')||(c>'9') ){
                    parseError("Number expected");
                    return 1;
                 }
            case '0':case '1':case '2':case '3':case '4':
            case '5':case '6':case '7':case '8':case '9':
              nf=c-'0';/*first char is scanned already!*/
              if(parseInt(&nf,theScan)==';')
                 break;
            default:                 
                 halt(16,"Parse error: number of f with ; expected instead of %c\n",c);
                 return 1;
   }/*switch(c)*/
   nf++;/*f counted from 1, not 0!*/

#ifdef MIXED_ARITHMETIC
   theScan->maxX=malloc(sizeof(SC_INT)*(nx+1));
   memset(theScan->maxX,0,sizeof(SC_INT)*(nx+1));
#endif

   theScan->nx=nx;
   theScan->nf=nf;

   ctTriadInit(estimateHashSize(theScan),&theScan->ctTriad);

#ifdef MIXED_ARITHMETIC
   theScan->fNativeLine=initCollectNativeFloat(INIT_FLOAT_SIZE, NULL);
   theScan->constStrings=initCollectInt(NULL);
   initMPool(&(theScan->newConstantsPool),1024);
   theScan->allMPvariables=initCollect(NULL);
#else
   theScan->fline=initCollectFloat(INIT_FLOAT_SIZE, NULL);
#endif

   if(initTrie(&(theScan->newConstantsTrie))<0)
      halt(-20,"initTrie failed -- not enough memory?\n");

   /*0 and 1 are pre-set:*/
   tAndSConst("0",NULL,0,theScan);
   /*also 1/1 since we replace / by *: for case ##/1*/
   tAndSConst("1",NULL,1,theScan);

   theScan->pstack=initCollectInt(NULL);

   theScan->condChain=initCollectInt(NULL);
   addInt(theScan->condChain,0);/*First element will be used as a buffer length*/

   theScan->allCondChain=initCollectInt(NULL);

   theScan->allocatedCondChains=initCollect(NULL);
   pushCell(theScan->allocatedCondChains,NULL);
   /*indexes must be counted from 1! 0 is a NULL pointer*/

   theScan->f=malloc(theScan->nf * sizeof(SC_INT));

   if(theScan->f==NULL)
      halt(-10,"malloc failed\n");

   /*Scan all f first:*/
   for(i=1;i<nf;i++){
      if(simpleScan(0,theScan)<0)
         return 1;/*Failure*/
      /*Store result and remove it from the stack:*/
      theScan->f[i]=popInt(theScan->pstack);
   }/*for(i=1;i<nf;i++)*/
   if (simpleScan(0,theScan)<0)
      return 1;/*Failure*/
   if(theScan->ctTriad.free - 1 != theScan->pstack->buf[theScan->pstack->fill-1] )
   /*Last triad is not the resulting one (e.g., first "if", or pure number)*/
      addCpyOp(theScan);
   /*Now convert ctTrian into rtTriad:*/
   hash_destroy(&(theScan->ctTriad.hashTable));

#ifdef MIXED_ARITHMETIC
   i=getPrecision(theScan);
   if (g_default_precisionChanged == 0)
      g_default_precision=i;

   mpfr_set_default_prec(g_default_precision);

   theScan->fline=initCollectFloat(theScan->fNativeLine->fill, NULL);

   addInt(theScan->constStrings,0);/*To be on the safe side*/

   for(i=0; i <theScan->fNativeLine->fill; i++){
      char *str=(char*) (theScan->newConstantsPool.pool + theScan->constStrings->buf[i]);
      pushCell(theScan->allMPvariables,theScan->fline->buf+theScan->fline->fill);
      addInternalFloatStr(theScan->fline,str);
      if(theScan->constStrings->buf[i+1] == -1){
         /*Store also */
         /*Initialise the last cell -- note, no reallocation may occur!:*/
         reservFloat(theScan->fline);
         /*fline->buf[1] is 1, fline->buf[fline->fill-1] is 'val', we
           store 1/val into fline->buf[fline->fill]:*/
#if ARITHMETIC == MPFR      
         mpfr_div(theScan->fline->buf[theScan->fline->fill],
               theScan->fline->buf[1],
               theScan->fline->buf[theScan->fline->fill-1],GMP_RNDN);
#endif
         pushCell(theScan->allMPvariables,theScan->fline->buf+theScan->fline->fill);
         theScan->fline->fill++;/*Store the cell*/
         i++;/*Skip next "-1" in theScan->constStrings*/
      }/*if(theScan->constStrings->buf[i+1] == -1)*/
   }/*for(i=0; i <theScan->fNativeline->fill; i++)*/
/*@@@*/
#if 0
   if(theScan->constStrings != NULL){
     free(theScan->constStrings->buf);
     free(theScan->constStrings);
     theScan->constStrings=NULL;
   }
   if(theScan->newConstantsPool.full>0){
      destroyMPool(&(theScan->newConstantsPool));
      theScan->newConstantsPool.full=0;
   }
#endif
/*@@@*/
#endif
   destroyTrie(&(theScan->newConstantsTrie));
   if(buildRtTriad(theScan))
      return 1;/*Failure*/
   ctTriadFree(&theScan->ctTriad);
   /*Now in theScan->rtTriad we have translated triads.*/
   /*The only p`roblem -- with the x array, it must be copied to 
     theScan->x before evaluation.*/
   return 0;
}/*masterScan*/
