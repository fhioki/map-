/**
 * Copyright (c) 2014 mdm <marco.masciola@gmail.com>
 *
 * This file is part of MAP++.
 *
 * MAP++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MAP++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MAP++. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _MAPSYS_H
#define _MAPSYS_H

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#if defined(_WIN32) || defined(_WIN64)
#  include <Windows.h>
#  include <tchar.h>
#else
#  include <unistd.h>
#endif


#ifdef _WIN64 
#  include "stdbool.h"
#  define map_snprintf _snprintf
#  define map_strcat(a,b,c) strcat_s(a,b,c)
#  define map_end_color(a)
#  define MAP_EXTERNCALL __declspec( dllexport )
#elif _WIN32 
#  include "stdbool.h"
#  define map_snprintf _snprintf
#  define map_strcat(a,b,c) strcat_s(a,b,c)
#  define map_end_color(a)
#  define MAP_EXTERNCALL __declspec( dllexport )
#else
#  include <stdbool.h>
#  define map_snprintf snprintf
#  define map_strcat(a,b,c) strncat(a,c,b)
#  define map_end_color(a) map_end_unix_color(a)
#  define MAP_EXTERNCALL 
#endif


#ifndef BUILD_DEFS_H
#  define BUILD_DEFS_H
#  define BUILD_YEAR_CH0 (__DATE__[ 7])
#  define BUILD_YEAR_CH1 (__DATE__[ 8])
#  define BUILD_YEAR_CH2 (__DATE__[ 9])
#  define BUILD_YEAR_CH3 (__DATE__[10])
#  define BUILD_MONTH_IS_JAN (__DATE__[0] == 'J' && __DATE__[1] == 'a' && __DATE__[2] == 'n')
#  define BUILD_MONTH_IS_FEB (__DATE__[0] == 'F')
#  define BUILD_MONTH_IS_MAR (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'r')
#  define BUILD_MONTH_IS_APR (__DATE__[0] == 'A' && __DATE__[1] == 'p')
#  define BUILD_MONTH_IS_MAY (__DATE__[0] == 'M' && __DATE__[1] == 'a' && __DATE__[2] == 'y')
#  define BUILD_MONTH_IS_JUN (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'n')
#  define BUILD_MONTH_IS_JUL (__DATE__[0] == 'J' && __DATE__[1] == 'u' && __DATE__[2] == 'l')
#  define BUILD_MONTH_IS_AUG (__DATE__[0] == 'A' && __DATE__[1] == 'u')
#  define BUILD_MONTH_IS_SEP (__DATE__[0] == 'S')
#  define BUILD_MONTH_IS_OCT (__DATE__[0] == 'O')
#  define BUILD_MONTH_IS_NOV (__DATE__[0] == 'N')
#  define BUILD_MONTH_IS_DEC (__DATE__[0] == 'D')
#  define BUILD_MONTH_CH0 (__DATE__[ 0])
#  define BUILD_MONTH_CH1 (__DATE__[ 1])
#  define BUILD_MONTH_CH2 (__DATE__[ 2])
#  define BUILD_DAY_CH0 ((__DATE__[4] >= '0') ? (__DATE__[4]) : '0')
#  define BUILD_DAY_CH1 (__DATE__[ 5])
#endif // BUILD_DEFS_H


#ifdef DEBUG
#  define checkpoint() printf("Checkpoint: Line %d in file %s\n",__LINE__,__FILE__);
#else
#  define checkpoint() 
#endif // DEBUG


#define PROGNAME "MAP++ (Mooring Analysis Program++)"
#define PROGVERSION "1.00.00a"
#define CHECKERRQ(code) if(success!=MAP_SAFE) { *ierr=map_set_universal_error("", map_msg, ierr, code ); break; } 
#define CHECKERRK(code) if(success!=MAP_SAFE) { *ierr=map_set_universal_error("", map_msg, ierr, code ); } 
#define MAPFREE(obj) if(obj!=NULL) { free( obj ); obj=NULL; } 
#define MAP_STR_LEN 1024
#define MapReal double /* this can be set to double, long double, ect as necessary */
#define DEG2RAD 0.01745329251 /*  pi/180  */
#define RAD2DEG 57.2957795131 /*  180/pi  */
#define ARCSINH(x) log(x+sqrt(1+x*x))
#define SPACE_LENGTH 12
#define MACHINE_EPSILON 1e-16
#define MAP_RETURN if(*ierr==MAP_FATAL) { return MAP_FATAL; }; return MAP_SAFE;


/**
 * Text Coloring (OS dependant)
 *
 * Used for text coloring in the terminal. If we are on a non-Unix OS, then:
 *   -- set the strings to "" (empty) so that garbage is not printed
 */
#ifdef __posix 
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#elif __linux
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#elif __unix
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#elif __APPLE__
#  define MAP_COLOR_RED "\033[1;31m"
#  define MAP_COLOR_YELLOW "\033[1;33m"
#  define MAP_COLOR_BLUE "\033[1;34m"
#  define MAP_COLOR_END "\033[0m"
#else
#  define MAP_COLOR_RED ""
#  define MAP_COLOR_YELLOW ""
#  define MAP_COLOR_BLUE ""
#  define MAP_COLOR_END ""
#endif

void __get_machine_name( char* machineName );

#endif /* _MAPSYS_H */
