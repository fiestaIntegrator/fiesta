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
  Inline implementation of the generic queues, both LIFO and FIFO.
  A queue is a cyclic buffer with two pointers, "head" and "tail".
  if head==tail then the queue is empty. If head+1==tail, the queue is full.
  Elements are always extracted from the tail.
  If new elements are plased to the head, then the queue is a FIFO queue,
  if new elements are plased to the tail, then the queue is a LIFO queue.

  The queue is created by the constructor qFLInit. The last argument of
  the constructor "isFixed", defines the behaviour of the queue when the
  buffer is full:

  If isFixed ==0, the queue is autmatically reallocates the full buffer.
  If isFixed>0, the queue overwrites the buffer when it is full.
  If isFixed <0, the push operation fails with -2 when the buffer is full.
  Corresponding macros may be used:
  QFL_REALLOC_IF_FULL, QFL_OVERWRITE_IF_FULL and QFL_FAIL_IF_FULL.

  Public functions:
  qFLInit -- constructor;
  qFLDestroy -- destructor;
  qFLReset -- empty the queue without destroying it;
  qFLPop -- extract one element from the queue;
  qFLPushFifo -- put a new element to the FIFO queue;
  qFLPushLifo -- put a new element to the LIFO queue.

  The macro QFL0 is a static initializer for the structure qFL_struct,
  which may be used to prevent comiler warnings:
  qFL_t q=QFL0;
*/

#ifndef QUEUE_H
#define QUEUE_H 1

#include <stdlib.h>
#include <string.h>

#include "comdef.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef COM_INLINE
#define QFL_INLINE COM_INLINE
#else
#define QFL_INLINE inline
#endif

#ifdef COM_INT
#define FL_INT COM_INT
#else
#define FL_INT int
#endif

/*
  If isFixed>0, the queue overwrites the buffer when it is full.
  If isFixed <0 -- when the buffer is full, push operation fails with -2:
*/
#define QFL_FAIL_IF_FULL (-1)
#define QFL_OVERWRITE_IF_FULL 1
#define QFL_REALLOC_IF_FULL 0

#define QFL_MAX_LENGTH 65536

#define QFL_INI_LENGTH 1024


#ifndef qFLAlloc
#define qFLAlloc(theSize) malloc(theSize)
#define qFLFree(theScratch,theLength) free(theScratch)
#define qFLAllocFail NULL
#endif

typedef struct qFL_struct{
   void **scratch;
   FL_INT len;
   /*if head==tail then the queue is empty. If head+1==tail, the queue is full.*/
   FL_INT head;/*Put new cells for Fifo*/
   FL_INT tail;/*Pop cells and put Lifo*/
   int isFixed;/*0 if reallocatable, otherwise, non-rellocatable:
                 >0 -- overwrite when full, <0 -- fail with -2*/
}qFL_t;
/*Static initializer for the structure qFL_struct:*/
#define QFL0 {NULL,0,0,0,0}

/*Expects theFL->head == theFL->tail:*/
static QFL_INLINE int qFLRealloc(qFL_t *theFL)
{
int oldLen=theFL->len;
void** tmp;

   if(theFL->len >= QFL_MAX_LENGTH )
      return -1;

   theFL->len *= 2;

   tmp=(void**) qFLAlloc(theFL->len*sizeof(void*));
   if(tmp == qFLAllocFail)
       return -1;
   memcpy(tmp,theFL->scratch,theFL->head*sizeof(void*));
   memcpy(tmp+theFL->head+oldLen,theFL->scratch+theFL->tail,
           (oldLen-theFL->tail)*sizeof(void*));
   theFL->tail+=oldLen;
   qFLFree(theFL->scratch,oldLen*sizeof(void*));
   theFL->scratch=tmp;
   return 0;
}/*qFLRealloc*/

/*Puts the cell to the head:*/
static QFL_INLINE FL_INT qFLPushFifo(qFL_t *theFL, void *cell)
{
FL_INT ret=theFL->head;
   theFL->scratch[ret]=cell;
   theFL->head++;
   if(theFL->head == theFL->len)
      theFL->head=0;
    if(theFL->head == theFL->tail){/*The buffer is full*/
      if(theFL->isFixed == 0){
         if(qFLRealloc(theFL))
            return -1;
         return ret;
      }/*if(theFL->isFixed == 0)*/
      if(theFL->isFixed>0){/*next time overwrite the tail*/
         if(++theFL->tail== theFL->len)
           theFL->tail=0;
         return ret;
      }
      /*theFL->isFixed<0 -- rollback head and fail with -2:*/
      if(--theFL->head<0)
         theFL->head=theFL->len-1;
      return -2;
   }/*if(theFL->head == theFL->tail)*/
   return ret;
}/*qFLPushFifo*/

/*Do LIFO placement (puts the cell to the tail):*/
static QFL_INLINE FL_INT qFLPushLifo(qFL_t *theFL, void *cell)
{
   if(--theFL->tail <0)
      theFL->tail=theFL->len-1;
   if(theFL->head == theFL->tail){/*The buffer is full*/
      if(theFL->isFixed == 0){/*realloc*/
         if(qFLRealloc(theFL))
            return -1;
      }/*if(theFL->isFixed == 0)*/
      else if(theFL->isFixed>0){/*overwrite the head*/
         if(--theFL->head<0)
            theFL->head = theFL->len-1;
      }/*else if(theFL->isFixed>0)*/
      else{/*theFL->isFixed<0, rollback the tail and fail*/
         if(++theFL->tail == theFL->len)
            theFL->tail = 0;
         return -2;
      }
   }/*if(theFL->head == theFL->tail)*/
   theFL->scratch[theFL->tail]=cell;
   return theFL->tail;
}/*qFLPushLifo*/

static QFL_INLINE void *qFLPop(qFL_t *theFL)
{
   void *ret;

   if(theFL->head == theFL->tail)
      return NULL;/*empty queue*/
   ret=theFL->scratch[theFL->tail];
   if(++theFL->tail == theFL->len)
      theFL->tail=0;
   return ret;
}/*qFLPop*/

/*Initializes the queue theFL. If it iz NULL, allocates it.
  If scratch==NULL, allocates the scratch, otherwise, use the given array.
  DO NOT CALL destroy, DO NOT USE isFixed == 0 for the external array!
  If iniLength> 0, uses this value for initialization, otherwise, uses default value.
  If isFixed ==0, the queue is autmatically reallocates the full buffer.
  If isFixed>0, the queue overwrites the buffer when it is full.
  If isFixed <0 -- when the buffer is full, push operation fails with -2:
*/
static QFL_INLINE qFL_t *qFLInit(void ** scratch,FL_INT iniLength, qFL_t *theFL,int isFixed)
{
  if(theFL == NULL){
     theFL=malloc(sizeof(qFL_t));
     if(theFL == NULL)
        return NULL;
  }
  theFL->isFixed=isFixed;
  if(iniLength<=0)
     theFL->len=QFL_INI_LENGTH;
  else
     theFL->len=iniLength;
  theFL->head=theFL->tail=0;
  if(scratch!=NULL)
     theFL->scratch=scratch;
  else{
     theFL->scratch=(void**) qFLAlloc( (theFL->len)*sizeof(void*));
     if(theFL->scratch == qFLAllocFail)
        return NULL;
  }
  return theFL;
}/*qFLInit*/

static QFL_INLINE void qFLReset(qFL_t *theFL)
{
   theFL->head=theFL->tail=0;
}/*qFLReset*/

static QFL_INLINE void qFLDestroy(qFL_t *theFL)
{
   qFLFree(theFL->scratch,theFL->len*sizeof(void*));
   theFL->scratch=NULL;
   theFL->len=0;
   theFL->head=theFL->tail=0;
}/*qFLDestroy*/

#endif


/*#define EXAMPLE*/

#ifdef EXAMPLE

#include <stdio.h>

int main(void)
{
int i,j,k;
char c[128];
qFL_t q;
int buf[1024];
   for(i=0;i<1024;i++)
      buf[i]=i;

   qFLInit(3,&q,QFL_REALLOC_IF_FULL);

/*   
#define QFL_FAIL_IF_FULL (-1)
#define QFL_OVERWRITE_IF_FULL 1
#define QFL_REALLOC_IF_FULL 0
*/ 
   for(;;){
      printf("+number to push, -number to stack, p to pop, r to reset, q to quit:");

      fgets(c,128,stdin);
      if(c[0]=='q')
         break;
      if(c[0]=='p'){
         int *a=(int*)qFLPop(&q);
         if(a!=NULL)
            printf(" Pop: %d\n",*a);
         else
            printf(" Pop: NULL\n");
         continue;
      }
      if(c[0]=='r'){
         printf(" Reset\n");
         qFLReset(&q);
         continue;
      }
      i=0;
      sscanf(c,"%d",&i);
      if(i<0){
         j=qFLPushLifo(&q, buf-i);
         printf(" Stack %d to %d\n",-i,j);
      }else{
         j=qFLPushFifo(&q, buf+i);
         printf(" Fifo %d to %d\n",i,j);
      }
   }/*for(;;)*/
   qFLDestroy(&q);
   return 0;
   
}
#ifdef __cplusplus
}
#endif

#endif
