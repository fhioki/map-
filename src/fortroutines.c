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
 * SUBROUTINE MAP_set_summary_file_name(interf,msg,err) bind(C,name='set_summary_file_name')
 *
 * @todo: need to free summaryFileName. This is done in delete_all_init_data(...), 
 *        should be called in Fortran routines
 */
MAP_EXTERNCALL void set_summary_file_name ( MAP_InitInputType_t* initFortType, char *map_msg, MAP_ERROR_CODE *ierr ) {  
  int length=0;
  char *temp=NULL;
  InitializationData* init=initFortType->object; 

  if(init->summaryFileName!=NULL) MAPFREE(init->summaryFileName);

  temp=strtok( initFortType->summaryFileName, " \n&");
  length=strlen(temp)+1;
  init->summaryFileName=(char*)malloc(sizeof(char)*length);
  strcpy(init->summaryFileName, temp);
  initFortType->summaryFileName[0]=0;  
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
  init->fullNodeInputString=NULL; 
  init->fullElementInputString=NULL;
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


/**
 *
 */
MAP_ERROR_CODE get_iteration_output_stream(MAP_OutputType_t *yType, MAP_OtherStateType_t* otherType, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int count=0;
  ModelData* data=otherType->object;
  VarTypePtr* varPtr=NULL;
  VarType* var=NULL;

  if (!yType->wrtOutput) { /* if NULL, then wrtOutput is not initialized */
    /* first set size of out list */
    yType->wrtOutput_Len = (int)list_size(&data->yList->outListPtr);
    yType->wrtOutput_Len += (int)list_size(&data->yList->outList);
    
    /* allocation array */
    yType->wrtOutput = malloc(sizeof(double)*yType->wrtOutput_Len);
  };

  list_iterator_start(&data->yList->outListPtr);
  while (list_iterator_hasnext(&data->yList->outListPtr)) { 
    varPtr = (VarTypePtr*)list_iterator_next(&data->yList->outListPtr);
    yType->wrtOutput[count] = *varPtr->value;
    count++;
  };
  list_iterator_stop(&data->yList->outListPtr);     

  list_iterator_start(&data->yList->outList);
  while (list_iterator_hasnext(&data->yList->outList)) { 
    var = (VarType*)list_iterator_next(&data->yList->outList);
    yType->wrtOutput[count] = var->value;
    count++;
  };
  list_iterator_stop(&data->yList->outList);     
  return MAP_SAFE;
};


/**
 *
 */
MAP_EXTERNCALL void fcall_get_header_string(int* N, char** strArr, MAP_OtherStateType_t* otherType)
{ 
  int count=0;    
  ModelData* data=otherType->object;
  VarTypePtr* varPtr=NULL;
  VarType* var=NULL;

  list_iterator_start(&data->yList->outListPtr);
  while (list_iterator_hasnext(&data->yList->outListPtr)) { 
    varPtr = (VarTypePtr*)list_iterator_next(&data->yList->outListPtr);
    strcpy(strArr[count],varPtr->name);
    count++;
  };
  list_iterator_stop(&data->yList->outListPtr);     

  list_iterator_start(&data->yList->outList);
  while (list_iterator_hasnext(&data->yList->outList)) { 
    var = (VarType*)list_iterator_next(&data->yList->outList);
    strcpy(strArr[count],var->name);
    count++;
  };
  list_iterator_stop(&data->yList->outList);     
};


/**
 *
 */
MAP_EXTERNCALL void fcall_get_unit_string(int* N, char** strArr ,MAP_OtherStateType_t* otherType )
{ 
  int count=0;    
  ModelData* data=otherType->object;
  VarTypePtr* varPtr=NULL;
  VarType* var=NULL;

  list_iterator_start(&data->yList->outListPtr);
  while (list_iterator_hasnext(&data->yList->outListPtr)) { 
    varPtr = (VarTypePtr*)list_iterator_next( &data->yList->outListPtr );
    strcpy(strArr[count],varPtr->units);
    count++;
  };
  list_iterator_stop(&data->yList->outListPtr);     

  list_iterator_start(&data->yList->outList);
  while (list_iterator_hasnext(&data->yList->outList)) { 
    var = (VarType*)list_iterator_next( &data->yList->outList );
    strcpy(strArr[count],var->units);
    count++;
  };
  list_iterator_stop(&data->yList->outList);     
};
