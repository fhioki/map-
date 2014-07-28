/***************************************************************************
 *   Copyright (C) 2014 mdm                                                *
 *   marco[dot]masciola at gmail                                           *
 *                                                                         *
 *   MAP++ is free software; you can redistribute it and/or modify it      *
 *   under the terms of the GNU General Public License as published by     *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.           *
 ***************************************************************************/


#include "map.h"
#include "maperror.h"
#include "protos.h"
#include "pyprotos.h"
#include "MAP_Types.h"


/**
 * Fortran binding routine
 * SUBROUTINE MAP_SetCableLibraryData( interf ) bind(C,name='set_cable_library_data')   
 */
MAP_EXTERNCALL void set_cable_library_data(MAP_InitInputType_t* initFortType)
{
  InitializationData* init=(InitializationData*)initFortType->object; 
  int N=init->librarySize;
  int length=strlen(initFortType->libraryInputLine)+1;

  init->libraryInputString = (char**)realloc(init->libraryInputString, sizeof(char*)*(N+1));
  init->libraryInputString[N]=malloc(sizeof(char)*length);    
  strcpy(init->libraryInputString[N], initFortType->libraryInputLine);
  init->librarySize++;
};


/**
 * Fortran binding routine
 * SUBROUTINE MAP_SetNodeData( interf ) bind(C,name='set_node_data')
 */
MAP_EXTERNCALL void set_node_data(MAP_InitInputType_t* initFortType)
{
  InitializationData* init=initFortType->object; 
  int N=init->nodeSize;
  int length=strlen(initFortType->nodeInputLine)+1;
  
  init->nodeInputString = (char**)realloc(init->nodeInputString, sizeof(char*)*(N+1));
  init->nodeInputString[N]=malloc(sizeof(char)*length);    
  strcpy(init->nodeInputString[N], initFortType->nodeInputLine);
  init->nodeSize++;
};


/**
 * Fortran binding routine
 * SUBROUTINE MAP_SetElementData( interf ) bind(C,name='set_element_data')
 */
MAP_EXTERNCALL void set_element_data(MAP_InitInputType_t* initFortType)
{
  InitializationData* init=initFortType->object; 
  int N=init->elementSize;
  int length=strlen(initFortType->elementInputLine)+1;

  init->elementInputString = (char**)realloc(init->elementInputString, sizeof(char*)*(N+1));
  init->elementInputString[N]=malloc(sizeof(char)*length);    
  strcpy(init->elementInputString[N], initFortType->elementInputLine);
  init->elementSize++;
};


/**
 * Fortran binding routine
 * SUBROUTINE MAP_SetSolverOptions( interf ) bind(C,name='set_solver_options')
 */
MAP_EXTERNCALL void set_solver_options(MAP_InitInputType_t* initFortType)
{
  InitializationData* init=initFortType->object; 
  int N=init->solverOptionsSize;
  int length=strlen(initFortType->optionInputLine)+1;

  init->solverOptionsString = (char**)realloc(init->solverOptionsString, sizeof(char*)*(N+1));
  init->solverOptionsString[N]=malloc(sizeof(char)*length);    
  strcpy(init->solverOptionsString[N], initFortType->optionInputLine);
  init->solverOptionsSize++;
};


/**
 * Fortran binding routine
 * SUBROUTINE MAP_set_initinput_to_null(interf,msg,err) bind(C,name='set_init_to_null') 
 */
MAP_EXTERNCALL void set_init_to_null(MAP_InitInputType_t* initFortType, char *map_msg, MAP_ERROR_CODE *ierr )
{
  InitializationData* init = initFortType->object; 
  init->libraryInputString=NULL;
  init->nodeInputString=NULL;
  init->elementInputString=NULL;
  init->solverOptionsString=NULL;
  init->expandedNodeInputString=NULL; 
  init->expandedElementInputString=NULL;
  init->summaryFileName=NULL;  
  init->sizeOfFullNodeString=0; 
  init->sizeOfFullElementString=0; 
  init->librarySize=0;
  init->nodeSize=0;
  init->elementSize=0;
  init->solverOptionsSize=0;
};


/**
 * @todo no longer needs to be called in fortran since we are no longer setting up a 
 *       mirror variable in the C OtherStates struct. However, this need to be called 
 *       in Python since it's how we set environment variables interactively. 
 *
 * Fortran binding routine
 * SUBROUTINE MAP_set_depth( interf, val ) bind(C,name='fcall_set_sea_depth')  
 */
MAP_EXTERNCALL void fcall_set_sea_depth( MAP_ParameterType_t* param, MapReal d )
{
  param->depth=d;
};


/**
 * @todo no longer needs to be called in fortran since we are no longer setting up a 
 *       mirror variable in the C OtherStates struct. However, this need to be called 
 *       in Python since it's how we set environment variables interactively. 
 *
 * Fortran binding routine
 * SUBROUTINE MAP_set_sea_density( interf, val ) bind(C,name='fcall_set_sea_density')  
 */
MAP_EXTERNCALL void fcall_set_sea_density( MAP_ParameterType_t* param, MapReal rho )
{
  param->rhoSea=rho;
};


/**
 * @todo no longer needs to be called in fortran since we are no longer setting up a 
 *       mirror variable in the C OtherStates struct. However, this need to be called 
 *       in Python since it's how we set environment variables interactively. 
 *
 * Fortran binding routine
 * SUBROUTINE MAP_set_gravity( interf, val ) bind(C,name='fcall_set_gravity')  
 */
MAP_EXTERNCALL void fcall_set_gravity(MAP_ParameterType_t* param, MapReal gravity)
{
  param->g = gravity;
};
