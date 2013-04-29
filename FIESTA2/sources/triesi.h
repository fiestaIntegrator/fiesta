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

#ifndef TRIESI_H
#define TRIESI_H 1

#ifdef __cplusplus
extern "C" {
#endif

#include "comdef.h"

#ifdef COM_INLINE
#define TR_INLINE COM_INLINE
#else
#define TR_INLINE inline
#endif

#ifdef COM_INT
#define TR_INT COM_INT
#else
#define TR_INT int
#endif

/*
triesi.h is interface to prefix tree implementation in tries.h
*/


typedef struct trie_struct{
   TR_INT value; /*pointer to a pool position*/
   char key;
   TR_INT right;/*horizontal link*/
   TR_INT next;/*vertical link*/
}trie_t;

typedef struct mpool_struct{
   char *pool;
   TR_INT fill;/*actual filled*/
   TR_INT full;/*allocated*/
}mpool_t;

typedef struct trieRoot_struct{
   trie_t root;
   mpool_t mpool;
}trieRoot_t;

#ifdef __cplusplus
}
#endif

#endif
