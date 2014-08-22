/****************************************************************
 *   Copyright (C) 2014 mdm                                     *
 *   marco[dot]masciola[at]gmail                                *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#include "initialization.h"


/**
 * @file 
 */


MAP_ERROR_CODE initialize_fortran_types(MAP_InputType_t* u_type, 
                                        MAP_ParameterType_t* p_type, 
                                        MAP_ContinuousStateType_t* x_type, 
                                        MAP_ConstraintStateType_t* z_type, 
                                        MAP_OtherStateType_t* other_type, 
                                        MAP_OutputType_t* y_type, 
                                        MAP_InitOutputType_t* initout_type)
{
  /* parameters are skipped for now; they are set in fortran since depth, 
   * gravity and sea density are set by glue code 
   */

  /* inputs */
  u_type->x = NULL;     u_type->x_Len = 0;
  u_type->y = NULL;     u_type->y_Len = 0;
  u_type->z = NULL;     u_type->z_Len = 0;

  /* continuous state */
  x_type->dummy=-999.9;

  /* constraint state */  
  z_type->H = NULL;     z_type->H_Len = 0;
  z_type->V = NULL;     z_type->V_Len = 0;
  z_type->x = NULL;     z_type->x_Len = 0;
  z_type->y = NULL;     z_type->y_Len = 0;
  z_type->z = NULL;     z_type->z_Len = 0;

  /* other state */
  other_type->H = NULL;     other_type->H_Len = 0;
  other_type->V = NULL;     other_type->V_Len = 0;
  other_type->Ha = NULL;    other_type->Ha_Len = 0;
  other_type->Va = NULL;    other_type->Va_Len = 0;
  other_type->x = NULL;     other_type->x_Len = 0;
  other_type->y = NULL;     other_type->y_Len = 0;
  other_type->z = NULL;     other_type->z_Len = 0;
  other_type->xa = NULL;    other_type->xa_Len = 0;
  other_type->ya = NULL;    other_type->ya_Len = 0;
  other_type->za = NULL;    other_type->za_Len = 0;
  other_type->Fx_connect = NULL;    other_type->Fx_connect_Len = 0;
  other_type->Fy_connect = NULL;    other_type->Fy_connect_Len = 0;
  other_type->Fz_connect = NULL;    other_type->Fz_connect_Len = 0;
  other_type->Fx_anchor = NULL;    other_type->Fx_anchor_Len = 0;
  other_type->Fy_anchor = NULL;    other_type->Fy_anchor_Len = 0;
  other_type->Fz_anchor = NULL;    other_type->Fz_anchor_Len = 0;

  /* outputs */
  y_type->Fx = NULL;              y_type->Fx_Len = 0;
  y_type->Fy = NULL;              y_type->Fy_Len = 0;
  y_type->Fz = NULL;              y_type->Fz_Len = 0;
  y_type->wrtOutput = NULL;       y_type->wrtOutput_Len = 0;
  
  /* init outputs */
  initout_type->progName[0] = '\0' ;
  initout_type->version[0] = '\0';
  initout_type->compilingData[0] = '\0';
  initout_type->writeOutputHdr = NULL;     initout_type->writeOutputHdr_Len = 0;
  initout_type->writeOutputUnt = NULL;     initout_type->writeOutputUnt_Len = 0;

  return MAP_SAFE;
};


void initialize_init_type_to_null(MAP_InitInputType_t* init_type)
{
  /* initialize the native Fortran/C types */
  init_type->gravity = -999.9;
  init_type->seaDensity = -999.9;
  init_type->depth = -999.9;
  init_type->fileName[0] = '\0';
  init_type->summaryFileName[0] = '\0';
  init_type->libraryInputLine[0] = '\0';
  init_type->nodeInputLine[0] = '\0';
  init_type->elementInputLine[0] = '\0';
  init_type->optionInputLine[0] = '\0';
};



void initialize_init_data_to_null(InitializationData* init_data)
{
  /* initialize the MAP initialization internal data strcture */
  init_data->libraryInputString = bstrListCreate();
  init_data->nodeInputString = bstrListCreate();
  init_data->elementInputString = bstrListCreate();
  init_data->solverOptionsString = bstrListCreate();
  init_data->expandedNodeInputString = bstrListCreate();
  init_data->expandedElementInputString = bstrListCreate();
};


void initialize_model_data_to_null(ModelData* model_data)
{
  model_data->MAP_SOLVE_TYPE = -999;
  model_data->yList = NULL; 
  model_data->sizeOfCableLibrary = 0;
  model_data->sizeOfElements = 0;
  model_data->sizeOfNodes = 0;    
  initialize_inner_solve_data_defaults(&model_data->inner_loop);    
  initialize_outer_solve_data_defaults(&model_data->outer_loop);    
  initialize_vessel_to_null(&model_data->vessel);    
  initialize_model_option_defaults(&model_data->modelOptions);
};


size_t cable_element_meter(const void *el) 
{
  return sizeof(Element);
};


size_t cable_library_meter(const void *el) 
{
  return sizeof(CableLibrary);
};


size_t node_meter(const void *el) 
{
  return sizeof(Node);
};


MAP_ERROR_CODE allocate_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{ 
  data->yList = malloc(sizeof(OutputList)); 
  if (data->yList==NULL) {
    *ierr = map_set_universal_error(NULL, map_msg, *ierr, MAP_FATAL_46);    
    return MAP_FATAL;
  };
  return MAP_SAFE;
};


void initialize_inner_solve_data_defaults(InnerSolveAttributes* inner)
{
  inner->f_tol = 1e-6;
  inner->g_tol = 1e-6;
  inner->x_tol = 1e-6;
  inner->max_its = 500;
  inner->m = 2;
  inner->n = 2;
  inner->factor = 1.0E2;             
  inner->factor = 1.0E2;             
  inner->ldfjac = 2; 
  inner->mode = 1;             
  inner->nprint = 2;           
  inner->info = 0;             
};


void initialize_model_option_defaults(ModelOptions* options) 
{
  options->repeat_angle_size = 0;
  options->repeat_angle = NULL;
  options->integrationDt = 0.01;
  options->kbLm = 3.0E6;
  options->cbLm = 3.0E5;
  options->waveKinematics = false;
}; 


/* deallocated in free_outer_solve_data() */
void initialize_outer_solve_data_defaults(OuterSolveAttributes* outer) 
{
  outer->fd = BACKWARD_DIFFERENCE;
  outer->pg = false;
  outer->tol = 1e-6;
  outer->epsilon = 1e-3;
  outer->maxIts = 500;
  outer->jac = NULL;
  outer->x = NULL;
  outer->b = NULL;
  outer->l = NULL;
  outer->u = NULL;
  outer->y = NULL;
};


void initialize_vessel_to_null(Vessel* floater)
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
