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
#include "initialization.h"
#include "protos.h"

/**
 * @file 
 */


/** @addtogroup FortranCall */
/* @{ */
MAP_EXTERNCALL void set_sea_depth(MAP_ParameterType_t* p_type, const MapReal depth)
{
  p_type->depth = depth;
};


MAP_EXTERNCALL void set_sea_density(MAP_ParameterType_t* p_type, const MapReal rho)
{
  p_type->rhoSea = rho;
};


MAP_EXTERNCALL void set_gravity(MAP_ParameterType_t* p_type, double gravity)
{
  p_type->g = gravity;
};
/* @} */


void initialize_init_type_to_null(MAP_InitInputType_t* init_type)
{
  /* initialize the native Fortran/C types */
  init_type->gravity = -999.9;
  init_type->seaDensity = -999.9;
  init_type->depth = -999.9;
  init_type->fileName[0] = 0;
  init_type->summaryFileName[0] = 0;
  init_type->libraryInputLine[0] = 0;
  init_type->nodeInputLine[0] = 0;
  init_type->elementInputLine[0] = 0;
  init_type->optionInputLine[0] = 0;
};


void initialize_init_data_to_null(InitializationData* init_data)
{
  /* initialize the MAP initialization internal data strcture */
  init_data->libraryInputString = NULL;
  init_data->nodeInputString = NULL;
  init_data->elementInputString = NULL;
  init_data->solverOptionsString = NULL;
  init_data->expandedNodeInputString = NULL; 
  init_data->expandedElementInputString = NULL;
  init_data->summaryFileName = NULL;  
  init_data->sizeOfFullNodeString = 0; 
  init_data->sizeOfFullElementString = 0; 
  init_data->librarySize = 0;
  init_data->nodeSize = 0;
  init_data->elementSize = 0;
  init_data->solverOptionsSize = 0;
};


void initialize_model_data_to_null(ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  model_data->MAP_SOLVE_TYPE = -999;
  model_data->yList = NULL; 
  model_data->sizeOfCableLibrary = 0;
  model_data->sizeOfElements = 0;
  model_data->sizeOfNodes = 0;    
  initialize_solver_data_to_null(&model_data->outerSolveData, &model_data->solverData);    
  initialize_vessel_to_null(&model_data->vessel, map_msg, ierr);    
  initialize_model_options_to_defaults(model_data);
};


void initialize_solver_data_to_null(MinPackDataOuter* mpOuter, MinPackDataInner* mpInner)
{
  mpInner->m = 2;
  mpInner->n = 2;
  mpInner->factor = 1.0E2;             
  mpInner->factor = 1.0E2;             
  mpInner->ldfjac = 2; /* number of columns in fjac */
  mpInner->mode = 1;             
  mpInner->nprint = 2;           
  mpInner->info = 0;             

  mpOuter->jac = NULL;
  mpOuter->x = NULL;
  mpOuter->b = NULL;
  mpOuter->l = NULL;
  mpOuter->u = NULL;
  mpOuter->y = NULL;
};


void initialize_model_options_to_defaults(ModelData* data) 
{
  data->modelOptions.innerFTol = 1e-6;
  data->modelOptions.innerGTol = 1e-6;
  data->modelOptions.innerXTol = 1e-6;
  data->modelOptions.innerMaxIts = 500;
  data->modelOptions.sizeOfRepeatAngles = 0;
  data->modelOptions.repeatAngles = NULL;
  data->modelOptions.integrationDt = 0.01;
  data->modelOptions.kbLm = 3.0E6;
  data->modelOptions.cbLm = 3.0E5;
  data->modelOptions.waveKinematics = false;

  data->outerSolveData.fd = BACKWARD_DIFFERENCE;
  data->outerSolveData.pg = false;
  data->outerSolveData.tol = 1e-6;
  data->outerSolveData.epsilon = 1e-3;
  data->outerSolveData.maxIts = 500;
  data->outerSolveData.jac = NULL;
  data->outerSolveData.x = NULL;
  data->outerSolveData.b = NULL;
  data->outerSolveData.l = NULL;
  data->outerSolveData.u = NULL;
}; 


void initialize_vessel_to_null(Vessel* floater, char* map_msg, MAP_ERROR_CODE* ierr)
{

  floater->xi = NULL;
  floater->yi = NULL;
  floater->zi = NULL;
    
  floater->displacement.x.name = NULL;
  floater->displacement.y.name = NULL;
  floater->displacement.z.name = NULL;
  floater->refOrigin.x.name = NULL;
  floater->refOrigin.y.name = NULL;
  floater->refOrigin.z.name = NULL;
  floater->lineSumForce.fx.name = NULL;
  floater->lineSumForce.fy.name = NULL;
  floater->lineSumForce.fz.name = NULL;
  floater->orientation.phi.name = NULL;
  floater->orientation.the.name = NULL;
  floater->orientation.psi.name = NULL;
    
  floater->displacement.x.units = NULL;
  floater->displacement.y.units = NULL;
  floater->displacement.z.units = NULL;
  floater->refOrigin.x.units = NULL;
  floater->refOrigin.y.units = NULL;
  floater->refOrigin.z.units = NULL;
  floater->lineSumForce.fx.units = NULL;
  floater->lineSumForce.fy.units = NULL;
  floater->lineSumForce.fz.units = NULL;
  floater->orientation.phi.units = NULL;
  floater->orientation.the.units = NULL;
  floater->orientation.psi.units = NULL;
};


MAP_EXTERNCALL InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* new_data = NULL;
  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(InitializationData));
  if (new_data == NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_4);    
    return new_data;
  } else {
    initialize_init_data_to_null(new_data);
    *ierr = MAP_SAFE;
    return new_data;
  };
};


MAP_EXTERNCALL MAP_InitInputType_t* map_create_init_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitInputType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = malloc(sizeof(MAP_InitInputType_t));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_4);    
  } else {    
    initialize_init_type_to_null(new_data); /* set F2C types to null */
    new_data->object = NULL;
    new_data->object = (InitializationData*)(uintptr_t)MAP_InitInput_Create(map_msg, ierr);
  };
  return new_data;
};


MAP_EXTERNCALL ModelData* MAP_OtherState_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* new_data = NULL;
  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(ModelData));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_43);    
    return new_data;
  } else {
    initialize_model_data_to_null(new_data, map_msg, ierr);
    *ierr = MAP_SAFE;
    return new_data;
  };
};

MAP_EXTERNCALL MAP_OtherStateType_t* map_create_other_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OtherStateType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = malloc(sizeof(MAP_OtherStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_43);    
  } else {
    new_data->object = NULL;
    new_data->object = (ModelData*)(uintptr_t)MAP_OtherState_Create(map_msg, ierr);
  };
  return new_data;    
};


MAP_EXTERNCALL MAP_InitOutputType_t* map_create_initout_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InitOutputType_t* new_type = NULL;  
  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_InitOutputType_t));
  if (new_type==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_11);    
  } else {
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_InputType_t* map_create_input_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_InputType_t* new_type = NULL;  
  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_InputType_t)); 
  if (new_type == NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_5);        
  } else { 
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_ParameterType_t* map_create_parameter_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ParameterType_t* new_type = NULL;  
  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_ParameterType_t));
  if (new_type==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_6);    
  } else {
    new_type->object = NULL;
  };
  return new_type;
};


MAP_EXTERNCALL MAP_ConstraintStateType_t* map_create_constraint_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ConstraintStateType_t* new_type = NULL;  
  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_ConstraintStateType_t)); 
  if (new_type==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_8);    
  } else {
    new_type->object = NULL;
  };
  return new_type;      
};


MAP_EXTERNCALL MAP_OutputType_t* map_create_output_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_OutputType_t* new_type = NULL;  
  *ierr = MAP_SAFE;
  new_type = malloc(sizeof(MAP_OutputType_t));
  if (new_type==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_10);    
  } else {
    new_type->object = NULL;
  };
  return new_type;    
};


MAP_EXTERNCALL MAP_ContinuousStateType_t* map_create_continuous_type(char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ContinuousStateType_t* new_data = NULL;  
  *ierr = MAP_SAFE;
  new_data = malloc(sizeof(MAP_ContinuousStateType_t)); 
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_7);    
  } else {
    new_data->object = NULL;    
  };
  return new_data;    
};
