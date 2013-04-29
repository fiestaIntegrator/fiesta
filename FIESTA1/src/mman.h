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
   Simple Windows stubs for mmap, munmap, getpagesize and off_t.
*/

#ifndef MMAN_H
#define MMAN_H
#include "comdef.h"
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PROT_READ 0
#define PROT_WRITE 0
#define MAP_PRIVATE 0
#define MAP_ANONYMOUS 0
#define MAP_FAILED NULL

#ifndef COM_INLINE
#define COM_INLINE inline
#endif

typedef int off_t;

static COM_INLINE void *mmap(void *start, 
                             size_t length, 
                             int prot , 
                             int flags, 
                             int fd, 
                             off_t offset)
{
   return malloc(length);
}/*mmap*/

static COM_INLINE int munmap(void *start, size_t length)
{
   free(start);
   return 0;
}/*munmap*/

static COM_INLINE int getpagesize(void)
{
return 4096;
}

#ifdef __cplusplus
}
#endif

#endif
