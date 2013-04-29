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

#ifndef COMDEF_H
#define COMDEF_H

//#define WIN 1

#ifdef WIN
#define VEGAS VEGAS
#define VEGAS1 VEGAS1
#define VEGAS2 VEGAS2
#define VEGAS3 VEGAS3
#define RESULT RESULT
#define COM_INLINE 
#else
#define VEGAS vegas_
#define VEGAS1 vegas1_
#define VEGAS2 vegas2_
#define VEGAS3 vegas3_
#define RESULT result_
#define COM_INLINE inline
#endif


#define  NO_FILES 1
#define  FOR_MATHEMATICA 1

#define WITH_OPTIMIZATION 1
#ifdef WITH_OPTIMIZATION
/*#define STATISTICS_OUT 1*/
#endif
/*#define COM_INT int*/
#define COM_INT long int
#define FLOAT double

#endif
