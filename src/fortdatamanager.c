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


extern const char MAP_ERROR_STRING[][1024];


/**
 *
 */
MAP_EXTERNCALL InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* new_data = NULL;
  map_reset_universal_error(map_msg, ierr);
  new_data = (InitializationData*)malloc( sizeof(InitializationData) );
  if (new_data == NULL) 
  {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_4);    
    return new_data;
  } else {
    *ierr = MAP_SAFE;
    return new_data;
  };
};


/**
 *
 */
MAP_EXTERNCALL void MAP_InitInput_Delete(InitializationData* init)
{
  MAPFREE( init ); 
};


// /**
//  * @todo: should not be called. Delete this in the fortran code
//  */
// MAP_EXTERNCALL void* MAP_Input_Create( char* map_msg, MAP_ERROR_CODE *ierr) 
// {
//   // InputData* new_data = NULL;  
//   // map_reset_universal_error( map_msg, ierr );
//   // new_data = (InputData*)malloc( sizeof(InputData) );
//   // if ( new_data == NULL ) 
//   // {
//   //   *ierr = map_set_universal_error( "", map_msg, ierr, MAP_FATAL_5 );    
//   //   return new_data;
//   // } else {
//   //   *ierr = MAP_SAFE;
//   //   return new_data;
//   // };
// };
//  
//  
// /**
//  * @todo: slated for removal; deprecated when above funtion is removed
//  */
// MAP_EXTERNCALL void MAP_Input_Delete(/*InputData* u*/ void* none)
// {
//   // MAPFREE( u );
// };


// /**
//  * @todo: should not be called. Delete this in the fortran code
//  */
// MAP_EXTERNCALL void* MAP_Param_Create( char* map_msg, MAP_ERROR_CODE *ierr )
// {
//   // ParameterData* new_data = NULL;
//   // map_reset_universal_error( map_msg, ierr );
//   // new_data = (ParameterData*)malloc( sizeof(ParameterData) );
//   // if ( new_data == NULL ) 
//   // {
//   //   *ierr = map_set_universal_error( "", map_msg, ierr, MAP_FATAL_6 );    
//   //   return new_data;
//   // } else {
//   //   *ierr = MAP_SAFE;
//   //   return new_data;
//   // };
// };
// 
// 
// /**
//  * @todo: slated for removal; deprecated when above funtion is removed
//  */
// MAP_EXTERNCALL void MAP_Param_Delete(/*ParameterData* p*/ void* none)
// {
//   // MAPFREE( p );
// };


// /**
//  * @todo: should not be called. Delete this in the fortran code
//  */
// MAP_EXTERNCALL void* MAP_ContState_Create( char* map_msg, MAP_ERROR_CODE *ierr )
// {
//   // ContinuousData* new_data = NULL;
//   // map_reset_universal_error( map_msg, ierr );
//   // new_data = (ContinuousData*)malloc( sizeof(ContinuousData) );
//   // if ( new_data == NULL ) 
//   // {
//   //   *ierr = map_set_universal_error( "", map_msg, ierr, MAP_FATAL_7 );    
//   //   return new_data;
//   // } else {
//   //   *ierr = MAP_SAFE;
//   //   return new_data;
//   // };
// };
// 
// 
// /**
//  * @todo: slated for removal; deprecated when above funtion is removed
//  */
// MAP_EXTERNCALL void MAP_ContState_Delete(/*ContinuousData* x*/ void* none)
// {
//   // MAPFREE( x );
// };


// /**
//  * @todo: should not be called. Delete this in the fortran code
// */
// MAP_EXTERNCALL void* MAP_ConstrState_Create( char* map_msg, MAP_ERROR_CODE *ierr )
// {
//   // ConstraintData* new_data = NULL;
//   // map_reset_universal_error( map_msg, ierr );
//   // new_data = (ConstraintData*)malloc( sizeof(ConstraintData) );
//   // if ( new_data == NULL ) 
//   // {
//   //   *ierr = map_set_universal_error( "", map_msg, ierr, MAP_FATAL_8 );    
//   //   return new_data;
//   // } else {
//   //   *ierr = MAP_SAFE;
//   //   return new_data;
//   // };
// };
// 
// 
// /**
//  * @todo: slated for removal; deprecated when above funtion is removed
//  */
// MAP_EXTERNCALL void MAP_ConstrState_Delete(/*ConstraintData* z*/ void* none)
// {
//   // MAPFREE( z );
// };


MAP_EXTERNCALL ModelData *MAP_OtherState_Create(char *map_msg, MAP_ERROR_CODE *ierr)
{
  ModelData* new_data=NULL;
  map_reset_universal_error(map_msg, ierr);
  new_data = (ModelData*)malloc( sizeof(ModelData) );
  if (new_data==NULL) 
  {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_43);    
    return new_data;
  } else {
    new_data->sizeOfCableLibrary=-9999;
    new_data->sizeOfElements=-9999;
    new_data->sizeOfNodes=-9999;    
    return new_data;
  };
};


/**
 * @todo: slated for removal; deprecated when above funtion is removed
 */
MAP_EXTERNCALL void MAP_OtherState_Delete(ModelData* data)
{
  MAPFREE(data);
};


// /**
//  * @todo: should not be called. Delete this in the fortran code
//  */
// MAP_EXTERNCALL void* MAP_Output_Create(char *map_msg, MAP_ERROR_CODE *ierr)
// {
//   // OutputData* new_data = NULL;
//   // map_reset_universal_error( map_msg, ierr );
//   // new_data = (OutputData*)malloc( sizeof(OutputData) );
//   // if ( new_data == NULL ) 
//   // {
//   //   *ierr = map_set_universal_error( "", map_msg, ierr, MAP_FATAL_10 );    
//   //   return new_data;
//   // } else {
//   //   *ierr = MAP_SAFE;
//   //   return new_data;
//   // };
// };
// 
// 
// /**
//  * @todo: slated for removal; deprecated when above funtion is removed
//  */ 
// MAP_EXTERNCALL void MAP_Output_Delete(/*OutputData* y*/ void* none)
// {
//   //MAPFREE( y );
// };


// /**
//  *
//  */
// MAP_EXTERNCALL void* MAP_InitOutput_Create(char *map_msg, MAP_ERROR_CODE *ierr)
// {
//   // InitializationOutputData* new_data = NULL;
//   // map_reset_universal_error(map_msg, ierr);
//   // new_data = (InitializationOutputData*)malloc( sizeof(InitializationOutputData) );
//   // if ( new_data == NULL ) 
//   // {
//   //   *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_11);    
//   //   return new_data;
//   // } else {
//   //   *ierr = MAP_SAFE;
//   //   return new_data;
//   // };
// };
// 
// 
// /**
//  *
//  */
// MAP_EXTERNCALL void MAP_InitOutput_Delete(/*InitializationOutputData* io*/void* none)
// {
//   // MAPFREE(io);
// };
