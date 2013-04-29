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

The special CIntegrate-specific hash table. It stores triads.

The main structure, hash_table_t, contains two arrays of hash_cell_t:
"table" and "scratch"; table is the main table, it stores triad
indexes (of the array ctTriad) in the fields table[h].triad, where "h"
is a hash address of the triad.  If table[h].triad == 0, the
corresponding triad is absent and to be installed to the array ctTriad
(this operation is performed not here but in the file scanner.c). If
the field table[h].next is not 0 then there are more triads with the
same hash address h. These additional cells are placed to the array
"scratch". The chain starts with table[h], and continues in
scratch[table[h].next]. The chain is terminated with
scratch[].next==0.

Sometimes the same hash table is used to store "condition chain",
which is the index of another buffer, available in this structure as
ctTriad->theScan->allocatedCondChains->buf. Indices of this buffer are
representetd as negative numbers, see e.g. the function
hash.c:hash_cmp().

In order to avoid possible memory segmentation, both the table and
scratch arrays allocated using mmap(), not malloc(). The array scratch
is reallocated when it becomes too small.

Public functions:
hash_init -- constructor, reurns the pointer to the initialized hash
table.  Parameters: theSize -- the size of the table, it will be
adjusted to the nearest prime number from the array hash.c:primes;
theTable -- the allocated hash table; ctTriad -- the array of triads,
see the file scanner.c.

hash_destroy(theTable) -- destructor; destroys the hash table
theTable.

hash_check(triad,theTable) -- returns the triad index if the triad is
installed in the hash table theTable, otherwise, 0.

hash_tand(triad,theTable) -- "Test And Set", always returns address of
a triad. If it is absent, the function installs it first. Note, the
function only installs the triad to the table, it does nothing with
the array ctTriad! Usual usage see the file scanner.c:tandsCtTriad() :
first, add a triad to the array ctTriad, then call the function
hash_tand.  If this function returns another number, this mens, this
triad already installed, and the triad counter will not be incremented
so the next triad will reuse the same element of the ctTriad array.
*/


#ifndef HASH_H
#define HASH_H

#include "comdef.h"
#include "scanner.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef COM_INLINE
#define HASH_INLINE COM_INLINE
#else
#define HASH_INLINE inline
#endif

#ifdef COM_INT
#define HASH_INT COM_INT
#else
#define HASH_INT int
#endif

/*From scanner.h:*/
struct ct_triad_struct;

typedef struct hash_cell_struct {
   HASH_INT triad;
   HASH_INT next;
}hash_cell_t;

typedef struct hash_table_struct {
   hash_cell_t *table;
   HASH_INT theSize;
   HASH_INT theSizeInBytes;
   hash_cell_t *scratch;
   HASH_INT free; /*first unused index of scratch*/
   HASH_INT max; /* number of allocated cells in scratch */
   HASH_INT allocatedInBytes; /*size of allocated scratch in bytes*/
   struct ct_triad_struct *ctTriad;
}hash_table_t;


hash_table_t *hash_init(HASH_INT theSize,hash_table_t *theTable, 
                         struct ct_triad_struct *ctTriad);


void hash_destroy(hash_table_t *theTable);

/*Returns address of a triad, or 0*/
HASH_INT hash_check(HASH_INT triad, hash_table_t *theTable);

/*TestAndSet, always returns address of a triad. If it is 
  absent, the function installs it first:*/
HASH_INT hash_tands(HASH_INT triad, hash_table_t *theTable);

#ifdef __cplusplus
}
#endif
#endif
