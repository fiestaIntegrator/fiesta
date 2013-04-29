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

/*The special CIntegrate-specific hash table, see comments in the file hash.h*/

#include "comdef.h"
/*mmap/munmap:*/
#ifdef WIN
#include "mman.h"
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

#include <string.h>

#include "scanner.h"

#include "hash.h"

/*#define DEBUG*/

HASH_INT primes[]={1021,2039,4093,8191,16381,32749,65521,131071,
                    262139,524287,1048573,2097143,4194301,8388593,33554393};
/*33554393 is the Prime[2063689], corresponds to 2^25 which is (*8) 256 MB*/

hash_table_t *hash_init(HASH_INT theSize,hash_table_t *theTable, ct_triad_t *ctTriad)
{
int ps=getpagesize();
int i;
   theTable->theSize=0;
   for(i=0; i<sizeof(primes)/sizeof(HASH_INT); i++)
      if(primes[i]>=theSize){
         theTable->theSize=primes[i];
         break;
      }/*if(primes[i]>=theSize)*/

   if(theTable->theSize==0)
      theTable->theSize=primes[i-1];
   
   theTable->theSizeInBytes=theTable->theSize*sizeof(hash_cell_t);

   if(theTable->theSizeInBytes % ps !=0 ){/*Not aligned to the page*/
      /*Page size alingment:*/
      theTable->theSizeInBytes=(theTable->theSizeInBytes/ps+1)*ps;
   }/*if(theTable->theSizeInBytes % ps !=0 )*/
   theTable->table=(hash_cell_t*)mmap (0, theTable->theSizeInBytes, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(theTable->table == MAP_FAILED)
      halt(-10, "mmap failed\n");

   memset(theTable->table, 0, theTable->theSizeInBytes);
   theTable->scratch=NULL;
   theTable->free=0;
   theTable->max=0;
   theTable->allocatedInBytes=0;
   theTable->ctTriad=ctTriad;
   return theTable;
}/*hash_init*/

void hash_destroy(hash_table_t *theTable)
{
   munmap(theTable->table,theTable->theSizeInBytes);
   munmap(theTable->scratch,theTable->allocatedInBytes);
   theTable->theSizeInBytes=0;
}/*hash_destroy*/

/*special comparison function to compare two "condition chains":*/
static HASH_INLINE int cmpAllocatedCondChains(SC_INT *v1, SC_INT *v2)
{
 int i=*v1-1;
   if(*v1!=*v2) return 0;
   for(;i>0;i--)if(*v1++ != *v2++)return 0;
   return 1;
}/*cmpAllocatedCondChains*/

/*The comparison function, compares two triads c1 and c2 (if c1>0 and c2>0)
of two condition chains (is c1<0 and c2<0) Returns 1 if triads co-incide, otherwise, 0:*/
static HASH_INLINE int hash_cmp( SC_INT c1, SC_INT c2,hash_table_t *theTable)
{ 
SC_INT *ct_firstOperand, *ct_secondOperand,
       *ct_triadType;
char   *ct_operation;

   if((c1<0)&&(c2<0)){/*The condition chain*/
      return cmpAllocatedCondChains(
                (SC_INT*)theTable->ctTriad->theScan->allocatedCondChains->buf[-c1],
                (SC_INT*)theTable->ctTriad->theScan->allocatedCondChains->buf[-c2]
             );
   }
   if((c1<0)||(c2<0))
      return 0;
   ct_firstOperand = theTable->ctTriad->firstOperand;
   ct_secondOperand =theTable->ctTriad->secondOperand;
   ct_triadType=theTable->ctTriad->triadType;
   ct_operation = theTable->ctTriad->operation;

   if(ct_triadType[c1] != ct_triadType[c2])
      return 0;

   if(ct_firstOperand[c1] == ct_firstOperand[c2]){
      if(ct_secondOperand[c1] == ct_secondOperand[c2])
         return (ct_operation[c1] == ct_operation[c2]);
      return 0;
   }else if(ct_operation[c1]!=ct_operation[c2])
      return 0;
   if(ct_operation[c1]>0)/*Non-commuting operation*/
      return 0;
   if(ct_firstOperand[c1] == ct_secondOperand[c2])
      return (ct_firstOperand[c2] == ct_secondOperand[c1]);
   return 0;
}/*hash_cmp*/

/*Calculates hash function for a triad c, if c>0 or the "condition chain" 
number -c, if c<0:*/
static HASH_INLINE HASH_INT hash_hash(HASH_INT c, hash_table_t *theTable)
{
   if(c<0){/*The condition chain*/
      SC_INT ret=5381;
      SC_INT *v=(SC_INT*)theTable->ctTriad->theScan->allocatedCondChains->buf[-c];
      SC_INT i=*v;
      /*Algotintm "Times 33 with XOR" DJBX33X by Daniel J. Bernstein:*/
      for(;i>0;i--){
         ret = ((ret << 5) + ret) ^ *v++;
      }/*for(;i>0;i--)*/
      return (unsigned HASH_INT)ret % theTable->theSize;
   }/*if(c<0)*/
   /*The triad.*/
   return (unsigned HASH_INT) (
           ( (theTable->ctTriad->triadType[c]+3) << 1 )+
           theTable->ctTriad->firstOperand[c]+
           theTable->ctTriad->secondOperand[c]+
           theTable->ctTriad->operation[c]
           ) % theTable->theSize;
}/*hash_hash*/

/*Returns address of a triad, or 0*/
HASH_INT hash_check(HASH_INT triad, hash_table_t *theTable)
{
   HASH_INT theIndex=hash_hash(triad,theTable);
   hash_cell_t *theElement= theTable->table+theIndex;

   if(theElement->triad != 0)for(;;){
      if( hash_cmp(triad, theElement->triad, theTable) )
         return theElement->triad;
      if(theElement->next == 0)
         break;
      theElement= theTable->scratch+theElement->next;
   }/*if(theElement->triad != 0)for(;;)*/
   return 0;
}/*hash_check*/

/*allocs the scratch array:*/
static HASH_INLINE void hash_alloc(hash_table_t *theTable)
{
   /*Allocate multiple both page size and the size of a cell,
     this will allow us each time on re-allocation just to double
     all the sizes:*/
   theTable->max = getpagesize();
   theTable->allocatedInBytes=theTable->max * sizeof(hash_cell_t);
   theTable->free = 1;/*Reserv 0 for a zero pointer*/
   theTable->scratch=(hash_cell_t*)mmap (0, theTable->allocatedInBytes, 
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(theTable->scratch == MAP_FAILED)
      halt(-10, "mmap failed\n");
}/*hash_alloc*/

/*Doubls the size of the scratch array:*/
static HASH_INLINE void hash_realloc(hash_table_t *theTable)
{
      HASH_INT newAllocatedInBytes=theTable->allocatedInBytes * 2;
      hash_cell_t *newScratch;

      newScratch=(hash_cell_t*)mmap (0, newAllocatedInBytes,
         PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
      if(newScratch== MAP_FAILED)
         halt(-10, "mmap failed\n");

      /*Note, re-allocation may happen only when all the scratch is consumed,
        so use allocatedInBytes instead of free*sizeof(hash_cell_t): */
      memcpy(newScratch,theTable->scratch,theTable->allocatedInBytes);
      if(munmap(theTable->scratch,theTable->allocatedInBytes))
         halt(-10, "munmap failed\n");

      theTable->max *=2;
      theTable->allocatedInBytes=newAllocatedInBytes;
      theTable->scratch=newScratch;
}/*hash_realloc*/

#ifdef DEBUG
HASH_INT maxC=0;
HASH_INT avC=0;
HASH_INT nC=0;
#endif

/*
"Test And Set", always returns address of
a triad. If it is absent, the function installs it first.
 Note, the
function only installs the triad to the table, it does nothing with
the array ctTriad! Usual usage see the file scanner.c:tandsCtTriad() :
first, add a triad to the array ctTriad, then call the function
hash_tand.  If this function returns another number, this mens, this
triad already installed, and the triad counter will not be incremented
so the next triad will reuse the same element of the ctTriad array.
*/

HASH_INT hash_tands(HASH_INT triad, hash_table_t *theTable)
{
   hash_cell_t *theElement= theTable->table+hash_hash(triad,theTable);
   int isScratch;
   if(theElement->triad == 0)/*No such a triad -- install it and return:*/
      return theElement->triad=triad;
   if(theElement->next == 0){
      if( hash_cmp(triad, theElement->triad, theTable) )
         return theElement->triad;
      isScratch=0;/*theElement points somewhere into the table*/
   }else{
#ifdef DEBUG
      HASH_INT maxCNew=0;
#endif
      isScratch=1;/*theElement points somewhere into the scratch*/
      for(;;){
         if( hash_cmp(triad, theElement->triad, theTable) )
            return theElement->triad;
         if(theElement->next == 0)
            break;
         theElement= theTable->scratch+theElement->next;
#ifdef DEBUG
         maxCNew++;
#endif
      }/*for(;;)*/
#ifdef DEBUG
      nC++;
      avC+=maxCNew;
      if(maxCNew>maxC)
         maxC=maxCNew;
#endif
   }/*else*/
   /*Install a new element in scratch:*/
   if(theTable->free == theTable->max){/*(re)allocate the scratch array:*/
      if(theTable->max == 0)
         hash_alloc(theTable);
      else{
         if(isScratch){/*theElement points somewhere into the scratch*/
            /*Store the position:*/
            HASH_INT n=theElement-theTable->scratch;
               hash_realloc(theTable);
               /*The scratch was changed, restore the position:*/
               theElement=theTable->scratch+n;
         }else/*theElement has no relation to the scratch, just realloc the scratch:*/
            hash_realloc(theTable);
      }/*else*/
   }/*if(theTable->free == theTable->max)*/
   /*Install the triad and return:*/
   theTable->scratch[theTable->free].triad=triad;
   theElement->next=theTable->free;
   theTable->scratch[theTable->free].next=0;
   theTable->free++;
   return triad;
}/*hash_tands*/

