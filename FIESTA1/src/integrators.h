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

#ifndef INTEGRATE_H
#define INTEGRATE_H 1

#include "scanner.h"
#include "runline.h"


#define MAX_INTEGRATOR 15

#include "vegasf.h"
#include <cuba.h>


#ifdef __cplusplus
extern "C" {
#endif


/*Use a global variable since there are no closures in C:*/
extern scan_t *g_fscan;/*defined in CIntegrate.c*/

/*The structure containing the answer and the error:*/
typedef struct result_struct{
    double s1; /* integral from all iterations */
    double s2; /* err. from all iterations */
}result_t;

typedef int (*integrator_t)(result_t *result);

extern integrator_t g_integratorArray[MAX_INTEGRATOR];
extern int g_currentIntegrator;
extern char *g_integratorNamesArray[MAX_INTEGRATOR];
extern int g_topIntegrator;

int initIntegrators(void);

int setPar(int, char *, char *);
int getAllPars(int, char *);

#ifdef __cplusplus
}
#endif
#endif
