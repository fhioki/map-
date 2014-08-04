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


#ifndef _INITIALIZATION_H
#define _INITIALIZATION_H

#include "MAP_Types.h"

MAP_EXTERNCALL void map_add_cable_library_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL void map_add_node_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL void map_add_element_input_text(MAP_InitInputType_t* init_type);
MAP_EXTERNCALL void map_add_options_input_text(MAP_InitInputType_t* init_type);


/**
 * @brief Set the water depth. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param depth water depth [m]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                                
 *    SUBROUTINE mapextern_set_depth(interf, fc_val) BIND(C,name='map_set_sea_depth')  
 *      IMPORT                            
 *      IMPLICIT NONE                     
 *      TYPE(MAP_ParameterType_C) interf
 *      REAL(C_DOUBLE), VALUE :: fc_val  
 *    END SUBROUTINE mapextern_set_depth    
 * END INTERFACE                      
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_depth(p%C_obj, depth)
 * @endcode
 */
MAP_EXTERNCALL void map_set_sea_depth(MAP_ParameterType_t* p_type, const MapReal depth);


/**
 * @brief Set the water density. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param rho water density [kg/m^3]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                                         
 *    SUBROUTINE mapextern_set_density(interf, fc_val) BIND(C,name='map_set_sea_density') 
 *      IMPORT                                            
 *      IMPLICIT NONE                                     
 *      TYPE(MAP_ParameterType_C) interf                
 *      REAL(C_DOUBLE), VALUE :: fc_val                      
 *    END SUBROUTINE mapextern_set_density                      
 * END INTERFACE                                          
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_set_density(p%C_obj, rho)
 * @endcode
 */
MAP_EXTERNCALL void map_set_sea_density(MAP_ParameterType_t* p_type, const MapReal rho);


/**
 * @brief Set the gravitational constant. Should be called before {@link map_init()}
 * @param p_type paramter type, native C struct {@link MAP_ParameterType_t}
 * @param grtavity gravitational acceleration [m/s^2]
 *
 * Example Fortran usage:
 * @code
 * INTERFACE                                                               
 *    SUBROUTINE MAP_map_set_gravity(interf, fc_val) BIND(C,name='map_set_gravity')
 *      IMPORT                                                         
 *      IMPLICIT NONE                                                  
 *      TYPE(MAP_ParameterType_C) interf                               
 *      REAL(C_DOUBLE), VALUE :: fc_val                                
 *    END SUBROUTINE MAP_map_set_gravity                                   
 * END INTERFACE                                                       
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_set_gravity(p%C_obj, g)
 * @endcode
 */
MAP_EXTERNCALL void map_set_gravity(MAP_ParameterType_t* p_type, const MapReal gravity);


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


/**
 * @brief   Allocate MAP_InitInputType_t and InitializationData
 * @details Called to allocate memory for the initialzation data for both the Fortran
 *          derived data and internal state data. Following sucessful allocation, 
 *          {@link initialize_init_type_to_null} and {@link initialize_init_data_to_null}
 *          are both called to nullify data. If not called, memory errors results. This should 
 *          the first function called when interacting with MAP. This is a necessary function for
 *          interaction with python and C based programs
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  initialization input type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_InitInputType_t* map_create_init_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_OtherStateType_t and ModelData
 * @details Called to allocate memory for the other states for both the Fortran
 *          derived data and internal state data. This is a necessary function for
 *          interaction with python and C based programs. The 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  other state type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_OtherStateType_t* map_create_other_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_InitOutputType_t 
 * @details Called to allocate memory for the initialization output type. The only obligation of
 *          this struct is to store the program version, necessary for FAST. This function is a
 *          necessary call for C and python, but can be ignored for Fortran if the MAP template 
 *          is followed (that is, ISO C Binding is followed in the mapping of Fortran type and C 
 *          structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  initialization output type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_InitOutputType_t* map_create_initout_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_InputType_t
 * @details Called to allocate memory for the input type. The program inputs are the fairlead 
 *          displacement due to the motion of the vessel the cable are attached to. This function 
 *          is a necessary call for C and python, but can be ignored for Fortran if the MAP 
 *          template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  input type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_InputType_t* map_create_input_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_ParameterType_t
 * @details Called to allocate memory for the parameter type. Parameters are time-invariant 
 *          constants, such as gravitational constant.  This function 
 *          is a necessary call for C and python, but can be ignored for Fortran if the MAP 
 *          template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  parameter type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_ParameterType_t* map_create_parameter_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_ConstraintType_t
 * @details Called to allocate memory for the constraint type. Constraints are variables solved
 *          through an algebraic equation. This is fairlead end forces (H and V) and node positions.  
 *          This function is a necessary call for C and python, but can be ignored for Fortran if
 *          the MAP template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  constraint type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_ConstraintStateType_t* map_create_constraint_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_OutputType_t
 * @details Called to allocate memory for the output type. IMPORTANT: this is different from the {@link OutList}.
 *          Output types are forces at the line fairlead only for lines connecting to the vessel. 
 *          This function is a necessary call for C and python, but can be ignored for Fortran if
 *          the MAP template is followed (that is, ISO C Binding is followed in the mapping of Fortran 
 *          type and C structures). 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  output type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_OutputType_t* map_create_output_type(char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief   Allocate MAP_ContinuousStateType_t
 * @details Called to allocate memory for the coninuous type. Not currently used, but it still is
 *          required to be allocated for FAST consistentcy. 
 * @param   map_msg, error message
 * @param   ierr, error code
 * @return  continuous type (equivalent C binding struct)  
 */
MAP_EXTERNCALL MAP_ContinuousStateType_t* map_create_continuous_type(char* map_msg, MAP_ERROR_CODE* ierr);


#endif /* _INITIALIZATION_H */
