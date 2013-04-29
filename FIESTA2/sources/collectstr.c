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
  Tools for manipulation by several types of 
  re-allocatable arrays, see comments in the file 
  collectstr.h
*/

#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "collectstr.h"

collect_t *initCollect( collect_t *theCollection)
{
   if(theCollection == NULL)
       theCollection=malloc(sizeof(collect_t));
   theCollection->buf=malloc( (theCollection->top=INIT_STR_SIZE)*sizeof(void*) );
   theCollection->fill=0;
   return theCollection;
}/*initCollect*/

void reallocCollect(collect_t *theCollection )
{
   theCollection->top*=2;
   theCollection->buf=realloc(theCollection->buf,(theCollection->top)*sizeof(void*));
}/*reallocCollec*/


collectStr_t *initCollectStr( collectStr_t *theSTring )
{
   if(theSTring == NULL)
       theSTring=malloc(sizeof(collectStr_t));
   *(theSTring->buf=malloc(theSTring->top=INIT_STR_SIZE))='\0';
   theSTring->fill=0;
   return theSTring;
}/*initCollectStr*/

void reallocCollectStr(collectStr_t *theSTring )
{
   theSTring->top*=2;
   theSTring->buf=realloc(theSTring->buf,theSTring->top);
}/*reallocCollectStr*/

collectInt_t *initCollectInt( collectInt_t *theInt )
{
   if(theInt == NULL)
       theInt=malloc(sizeof(collectInt_t));
   theInt->buf=malloc((theInt->top=INIT_INT_SIZE)*sizeof(CS_INT));
   theInt->fill=0;
   return theInt;
}/*initCollectInt*/

void reallocCollectInt(collectInt_t *theInt )
{
   theInt->top*=2;
   theInt->buf=realloc(theInt->buf,(theInt->top)*(sizeof(CS_INT)));
}/*reallocCollectInt*/

collectFloat_t *initCollectFloat( size_t iniSize, collectFloat_t *theFloat )
{
   if(theFloat == NULL)
       theFloat=malloc(sizeof(collectFloat_t));
   theFloat->buf=malloc((theFloat->top=iniSize)*sizeof(INTERNAL_FLOAT));
   theFloat->fill=0;
   return theFloat;
}/*initCollectFloat*/

void reallocCollectFloat(collectFloat_t *theFloat )
{
   theFloat->top*=2;
   theFloat->buf=realloc(theFloat->buf,(theFloat->top)*(sizeof(INTERNAL_FLOAT)));
}/*reallocCollectFloat*/

#ifdef MIXED_ARITHMETIC

void reallocCollectNativeFloat(collectNativeFloat_t *theFloat )
{
   theFloat->top*=2;
   theFloat->buf=realloc(theFloat->buf,(theFloat->top)*(sizeof(FLOAT)));
}/*reallocCollectNativeFloat*/

collectNativeFloat_t *initCollectNativeFloat( size_t iniSize, collectNativeFloat_t *theFloat )
{
   if(theFloat == NULL)
       theFloat=malloc(sizeof(collectNativeFloat_t));
   theFloat->buf=malloc((theFloat->top=iniSize)*sizeof(FLOAT));
   theFloat->fill=0;
   return theFloat;
}/*initCollectNativeFloat*/

#endif
