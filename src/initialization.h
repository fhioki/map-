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


#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H


#include "map.h"


MAP_ERROR_CODE initialize_fortran_types(MAP_InputType_t* u_type, MAP_ParameterType_t* p_type, MAP_ContinuousStateType_t* x_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, MAP_InitOutputType_t* initout_type);
MAP_ERROR_CODE allocate_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
size_t cable_library_meter(const void* el);
size_t node_meter(const void* el);
size_t cable_element_meter( const void* el);


/**
 * @brief   Initializes the Fortran/C iteroperability types
 * @details This is called in the map_create_init_type routine following successful allocation 
 *          of memory. This should not be called directly by the user code. 
 * @param   init_type, Fortran/C interoperable type {@link MAP_InitInputType_t}
 * @see     map_create_init_type()
 */
void initialize_init_type_to_null(MAP_InitInputType_t* init_type);


/**
 * @brief   Initializes MAP internal initialization data
 * @details This is called in the map_create_init_type routine following successful allocation 
 *          of memory. This should not be called directly by the user code. The internal states
 *          are nullified and set to zero. 
 * @param   init_data, internal MAP initialization data {@link InitializationData}
 * @see     map_create_init_type()
 */
void initialize_init_data_to_null(InitializationData* init_data);


/**
 * @brief   Initialized all model data parameters to NULL or default values
 * @details Calls other functions to nullify and assign model data variables default values.          
 * @param   floater, the vessel strcture
 * @see     map_create_other_type()
 */
void initialize_model_data_to_null(ModelData* model_data);


/**
 * @brief   Initializes MAP internal initialization data
 * @details This is called in the map_create_init_type routine following successful allocation 
 *          of memory. This should not be called directly by the user code. The internal states
 *          are nullified and set to zero. 
 * @param   options, model variables not related to solver items/options
 * @see     map_create_init_type()
 */
void initialize_model_option_defaults(ModelOptions* options);


/**
 * @brief   Sets inner-loop solver data 
 * @details Jacobian is 2x2 and residual vectors are 2x1. Sizes are fixed because only two
 *          algebraic equations are solved for each lines. 
 *          {@link free_outer_solve_data}
 * @param   inner, inner-loop solver data
 * @see     map_create_other_type()
 */
void initialize_inner_solve_data_defaults(InnerSolveAttributes* inner);


/**
 * @brief   Sets vessel parameters to NULL
 * @details This function nullifies Vessel-associated pointers. The Vessel structure is a
 *          convenient way to displace all fairleads connected to a vessel with a single 
 *          function call. Since the number of 'VESSEL' nodes are unknown at compile time,
 *          vessel arrays are sized at run-time. Deallocation occurs in {@link free_vessel}
 * @param   floater, the vessel strcture
 * @see     map_create_other_type()
 */
void initialize_vessel_to_null(Vessel* floater);


/**
 * @brief   Sets outer solver data (Jacobian, residual vector, ect) to NULL
 * @details Unlike the inner solver data strcture, Jacobian and residual vectors sizes are not 
 *          declared at compile time. These are resized depending on the number of 'connect'
 *          node declared in the problem. Outer solve data is deallocated in 
 *          {@link free_outer_solve_data}
 * @param   outer, outer-loop solver data
 * @see     map_create_other_type()
 */
void initialize_outer_solve_data_defaults(OuterSolveAttributes* outer);


/**
 * @brief   Allocate ModelData
 * @details Called by {@link  map_create_other_type} to allocate memory for the internal 
 *          state (model) data type. 'Other States', as FAST calls them, are states not 
 *          fitting a regular role as a parameter, constraint, input, ect. Other states
 *          contain information on the line connectivity matrix, how reference to poperties
 *          for each line, buoyancy properties of the nodes, ect. Deallocated using
 *          interaction with python and C based programs. Structure is free'd by calling
 *          {@link MAP_OtherState_Delete}.
 * @param   map_msg, error message
 * @param   ierr, error code
 * @see     map_create_other_type()
 * @return  instance of the interal model struct (different from the FAST-required derived types)  
 */
ModelData* MAP_OtherState_Create(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate InitializationData
 * @details Called by {@link  map_create_init_type} to allocate memory for the iinitialization
 *          data type. The reason why a layer is added to the initialization data is due to 
 *          Fortran interactions. It is straighforward to pass 1D character arrays between
 *          Fortran and C instead of 2D arrays. 2D arrays would make more sense since multiple 
 *          lines from the MAP input file can be packed in one step. {@link MAP_InitInputType_t}
 *          in responsible for the 1D arrays. which are passed from Fortran to C. MAP then takes
 *          the 1D aray and packs it into InitializationData. This is used to subsequently 
 *          initialize the model. Structure is free'd by calling {@link MAP_InitInput_Delete}.
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  instance of the packed initialization strings (different from the FAST-required derived types)  
 */
InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr);



#endif /* _INITIALIZATION_H */
