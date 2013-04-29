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

#ifndef TRIES_H
#define TRIES_H 1


#ifdef __cplusplus
extern "C" {
#endif

/*

This file contains online implementation of prefix tree for ASCII-Z
keys and non-negative integer 'values'. The storage is mmap-based
pool, there are also functions for placing arbitrary object to the
pool so the integer 'value' stored in the a trie may be a reference to
the informaton stored in the pool.

No node in the tree stores the key associated with that node; instead,
its position in the tree shows what key it is associated with. All the
descendants of a node have a common prefix of the string associated
with that node, and the root is associated with the empty
string. Values are not associated with every node, only with
leaves and some inner nodes that carry 'value'>=0.

All objects are allocated in the mmap-based pool, all "pointers" are
integer index of this pool, NULL pointer is -1.

The node structure is the following:

{value,key,next,right}

'key' is the currently variable positon of a key, 'right' (if any)
points to the next node with the same prefix, 'next' (if any) points
to the next node with the next prefix. Example:

NULL,'a',next NULL
           |
            -> NULL 'a' next,right -> NULL 'b' next NULL
                         |                     |
                         |                      ->3,'1',NULL,NULL
                          -> 1,'1',next,right -> 2,'b',NULL,NULL
                                    |
                                     -> 4,'x',NULL,NULL

This trie contains the following strings with associated numbers:
aa1  1
aa1x 4
aab  2
ab1  3

Chains through right links form ordered lists.

*/
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <errno.h>
#include <string.h>

#include <unistd.h>
#include <sys/mman.h>

#include "triesi.h"

static TR_INLINE TR_INT initMPool(mpool_t *mpool,int initSize)
{
   mpool->full=getpagesize()*initSize;
   mpool->fill=0;
   mpool->pool = (char*)mmap (0, mpool->full,
      PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
   if(mpool->pool == MAP_FAILED){
      halt(-10, "initMPool: mmap failed\n");
      return -1;
   }
   return mpool->full;
}/*initMPool*/

static TR_INLINE void destroyMPool(mpool_t *mpool)
{
   munmap(mpool->pool,mpool->full);
}/*destroyMPool*/

static TR_INT mpoolAdd(void *ptr, mpool_t *mpool,TR_INT theSize)
{
TR_INT tmp,origSize=theSize;
   if(theSize % 8)
      theSize=(theSize+8)/8*8;
   if( (mpool->fill+theSize)>=mpool->full ){
      char *newPool;
      tmp=mpool->full*2;
      newPool=(char*)mmap (0,tmp,
         PROT_READ|PROT_WRITE, MAP_PRIVATE|MAP_ANONYMOUS, 0, 0);
      if(newPool== MAP_FAILED){
         halt(-10, "mpoolAdd: mmap failed\n");
         return -1;
      }
      memcpy(newPool,mpool->pool,mpool->fill);
      munmap(mpool->pool,mpool->full);
      mpool->pool=newPool;
      mpool->full=tmp;
   }
   memcpy(mpool->pool+mpool->fill,ptr,origSize);
   tmp=mpool->fill;
   mpool->fill+=theSize;
   return tmp;
}/*mpoolAdd*/

static TR_INLINE int initTrie(trieRoot_t *trie)
{
   if(initMPool(&(trie->mpool),16)<0)
      return -1;
   trie->root.value=-1;
   trie->root.key='\0';
   trie->root.right=-1;
   trie->root.next=-1;
   return 1;
}/*initTrie*/

static TR_INLINE void destroyTrie(trieRoot_t *trie)
{
   destroyMPool(&(trie->mpool));
}/*destroyTrie*/

static TR_INT buildTrie(char *str, TR_INT value,trieRoot_t *trie,TR_INT right)
{
trie_t tbuf;
trie_t *t;
TR_INT root;

   tbuf.value=-1;
   tbuf.key=*str;
   tbuf.right=right;

   right=root=mpoolAdd(&tbuf, &(trie->mpool),sizeof(trie_t));
   while( *++str !='\0'){
      TR_INT tmp;
      tbuf.value=-1;
      tbuf.key=*str;
      tbuf.right=-1;
      tmp=mpoolAdd(&tbuf, &(trie->mpool),sizeof(trie_t));
      t=(trie_t *)(trie->mpool.pool+right);
      t->next=right=tmp;
   }
   t=(trie_t *)(trie->mpool.pool+right);
   t->next=-1;
   t->value=value;
   return root;
}

/*returns -1 if fails, 1 is install a new or 0 if found already existing.
In the latter case the value of '*value' is set to the existing value:*/
static int installTrie(char *str, TR_INT *value,trieRoot_t *trie)
{
trie_t *t=&(trie->root);
TR_INT mem,tmp;
   if(t->right<0){
      t->right=buildTrie(str,*value, trie, -1);
      return 1;
   }
   mem= ((char*)&(t->right)) - trie->mpool.pool;
   t=(trie_t *)(trie->mpool.pool+t->right);
   while(*str != '\0'){
      if (*str < t->key){
         tmp=buildTrie(str,*value, trie, *((TR_INT*)(trie->mpool.pool+mem)));
         *((TR_INT*)(trie->mpool.pool+mem))=tmp;
         return 1;
      }
      if(*str == t->key){
         str++;
         if(*str == '\0'){
            if(t->value >= 0){
               *value = t->value;
               return 0;
            }/*if(t->value >= 0)*/
            t->value = *value;
            return 1;
         }/*if(*str == '\0')*/
         if(t->next < 0){
            /*we can't do just 
              t->next = buildTrie(str,*value,trie,-1);
              since the actual address of t might be 
              changed by buildTrie
             */
            mem= ((char*)&(t->next))-trie->mpool.pool;
            tmp = buildTrie(str,*value,trie,-1);
            *((TR_INT*)(trie->mpool.pool+mem))=tmp;
            return 1;
         }/*if(t->next < 0)*/
         mem= ((char*)&(t->next))-trie->mpool.pool;
         t=(trie_t *)(trie->mpool.pool+t->next);
         continue;
      }/*if(*str == t->key)*/
      /* *str > t->key*/
      if(t->right < 0){
         /*we can't do just 
           t->right = buildTrie(str,*value,trie,-1);
           since the actual address of t might be 
           changed by buildTrie
          */
         mem= ((char*)&(t->right))-trie->mpool.pool;
         tmp=buildTrie(str,*value, trie, -1);
         *((TR_INT*)(trie->mpool.pool+mem))=tmp;
         return 1;
      }/*if(t->right < 0)*/
      mem= ((char*)&(t->right))-trie->mpool.pool;
      t=(trie_t *)(trie->mpool.pool+t->right);
   }/*while(*str != '\0')*/
   return -1;
}/*installTrie*/

/*returns -1 if fails, othervise, found value:*/
static TR_INT checkTrie(char *str, trieRoot_t *trie)
{
trie_t *t=&(trie->root);

   if(t->right<0)
      return -1;
   t=(trie_t *)(trie->mpool.pool+t->right);
   while(*str != '\0'){
      if (*str < t->key)
         return -1;
      if(*str == t->key){
         str++;
         if(*str == '\0'){
            if(t->value >= 0)
               return t->value;
            return -1;
         }/*if(*str == '\0')*/
         if(t->next < 0)
            return -1;
         t=(trie_t *)(trie->mpool.pool+t->next);
         continue;
      }/*if(*str == t->key)*/
      /* *str > t->key*/
      if(t->right < 0)
         return -1;
      t=(trie_t *)(trie->mpool.pool+t->right);
   }/*while(*str != '\0')*/
   return -1;
}/*checkTrie*/

#ifdef __cplusplus
}
#endif

#endif
