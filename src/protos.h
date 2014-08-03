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


#ifndef _PROTOS_H
#define _PROTOS_H

#include "cminpack/cminpack.h"
#include "MAP_Types.h"

/*
 * system prototypes
 */

MAP_ERROR_CODE set_vessel(Vessel* floater, const MAP_InputType_t* uType, char* map_msg, MAP_ERROR_CODE* ierr);
int strcicmp( char const* a, char const* b );
MAP_ERROR_CODE is_numeric ( const char* s );
void end_color( char* dest, const char* src );
MAP_ERROR_CODE map_get_version(MAP_InitOutputType_t* ioType);
MAP_ERROR_CODE map_set_universal_error( char* user_str, char* stat_str, const MAP_ERROR_CODE* current, const MAP_ERROR_CODE ierr );
MAP_ERROR_CODE print_help_to_screen();
void map_reset_universal_error( char* user_str, MAP_ERROR_CODE* ierr );

/*
 *
 */
MAP_ERROR_CODE associate_vartype_ptr(VarTypePtr* type, double* arr, int index);
MAP_ERROR_CODE set_vartype( char* unit, char* alias, const int num, VarType* type, char const* property );
MAP_ERROR_CODE set_vartype_ptr( char *unit, char *alias, const int num, VarTypePtr* type, char const* property );
MAP_ERROR_CODE set_vartype_ptr_float(char* unit, char* alias, const int num, VarTypePtr* type, const MapReal value);
MAP_ERROR_CODE set_element_vartype_ptr( char *unit, char *alias, const int num, VarTypePtr* type);

const char* remove_first_character(const char* str);
void map_end_unix_color( char* dest );
void initialize_element( Element* ptr, const int ElementNum );
MAP_ERROR_CODE initialize_external_applied_force(char* unit, char* alias, const int num, VarType* type, char const* property);
MAP_ERROR_CODE initialize_node_sum_force_ptr(char* unit, char* alias, const int num, VarTypePtr* type);
void initialize_vartype( char* unit, char* alias, VarType* type, const int num );

int set_cable_library_list( ModelData* data, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr );
int set_node_list(const MAP_ParameterType_t* pType,  MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, ModelData* data, char** const nodeInputString, char* map_msg, MAP_ERROR_CODE* ierr);
int set_element_list(MAP_ConstraintStateType_t* zType, ModelData* data, char** const elementInputString, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_model_options_list(ModelData* data, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);

size_t cable_library_meter(const void* el);
size_t node_meter(const void* el);
size_t vartype_meter(const void* el);
size_t vartype_ptr_meter(const void *el);
size_t cable_element_meter( const void* el);
void set_element_null(Element* ptr);
void set_node_null(Node* ptr);
MAP_ERROR_CODE associate_element_with_cable_property(Element* in_element, ModelData* data, char* word, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE associate_element_with_anchor_node(ModelData* data, Element* new_element, const int element_num, char* word, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE associate_element_with_fairlead_node(ModelData* data, Element* new_element, const int element_num, char* word, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE repeat_nodes(ModelData* dataObj, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE repeat_elements(ModelData* dataObj, InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);

/*
 * This sets the pointers to NULL for the vessel object and gives it default properties. Only 
 * to be used in the python glue code. 
 *
 * @todo: need to associate the node with inputs
 * @acceses: set_vartype_float( )
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE set_vartype_float(char* unit, char* alias, const int num, VarType* type, MapReal const value);


/*
 * @input: Fortran types: Input, Parameter, Continuous States, Constraint State, Other State, Output, Init Output 
 * @returns: error code
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE initialize_fortran_types(MAP_InputType_t* uType, MAP_ParameterType_t* pType, MAP_ContinuousStateType_t* xType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, MAP_InitOutputType_t* initoutType);
MAP_ERROR_CODE free_fortran_types(MAP_InputType_t* uType, MAP_ParameterType_t* pType, MAP_ContinuousStateType_t* xType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType);

MAP_ERROR_CODE allocate_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE free_outlist(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);

/*
 * Sets init data to NULL or -9999
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */
MAP_ERROR_CODE initialize_init_input(InitializationData* init);

MAP_ERROR_CODE first_solve(ModelData* data, MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, MAP_OutputType_t* yType, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE allocate_outer_solve_data(MinPackDataOuter* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE free_outer_solve_data(MinPackDataOuter* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE solve_line(ModelData* data, double time, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE node_solve_sequence(ModelData* data, MAP_InputType_t* uType, MAP_ConstraintStateType_t* zType, MAP_OtherStateType_t* otherType, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE line_solve_sequence(ModelData* otherType, double t, char* map_msg, MAP_ERROR_CODE* ierr);

/*
 * Initialized omega (weight per unit length) and cross-section area of a cable. The formula is 
 *
 *   A=\pi*\frac{radius^{2}}{4}
 *   \omega=g*(\mu-A*\rho_{seawater})
 *
 * @acceses: none
 * @calledby: mapcall_msqs_init( )
 */                                                   
MAP_ERROR_CODE initialize_cable_library_variables(ModelData *data, MAP_ParameterType_t* paramFortType, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_line_variables_pre_solve(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_element_initial_guess(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_line_variables_post_solve(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE calculate_node_sum_force(ModelData* data);
void print_machine_name_to_screen( );

/*
 * stuff for nodes...
 */
MAP_ERROR_CODE reset_node_force_to_zero(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
void add_to_sum_fx(Node* node, const MapReal fx);
void add_to_sum_fy(Node* node, const MapReal fy);
void add_to_sum_fz(Node* node, const MapReal fz);

/**
 * Numeric routines
 */
int inner_function_evals(void* elementPtr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag);
MAP_ERROR_CODE call_minpack_lmder(Element* element, MinPackDataInner* mp, ModelOptions* opt, const int lineNum, const double time, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE set_psi(Element* element, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lu(MinPackDataOuter* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lu_back_substitution(MinPackDataOuter* ns, const int n, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE forward_difference_jacobian(MAP_OtherStateType_t* otherType, MAP_ConstraintStateType_t* zType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE backward_difference_jacobian(MAP_OtherStateType_t* otherType, MAP_ConstraintStateType_t* zType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);






/*
 * MARK: numeric helper functions
 */
MapReal residual_function_length_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l);
MapReal residual_function_height_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h);
MapReal jacobian_dxdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dxdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dzdh_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);
MapReal jacobian_dzdv_no_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA);

MapReal residual_function_length_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal l, const MapReal cb);
MapReal residual_function_height_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal h, const MapReal cb);
MapReal jacobian_dxdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);
MapReal jacobian_dxdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);
MapReal jacobian_dzdh_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);
MapReal jacobian_dzdv_contact(const MapReal V, const MapReal H, const MapReal w, const MapReal Lu, const MapReal EA, const MapReal cb);

MapReal set_vertical_excursion(Element* element);
MapReal set_horizontal_excursion(Element* element);
MAP_ERROR_CODE check_maximum_line_length(Element* element, const bool contactFlag, char* map_msg, MAP_ERROR_CODE* ierr);
MapReal get_maximum_line_length(Element* element);
MAP_ERROR_CODE set_element_initial_guess(ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);

MAP_ERROR_CODE increment_dof_by_delta(double* inputType, const MapReal delta, const int size);
MAP_ERROR_CODE restore_original_displacement(double* inputType, const double* initialValue, const int size);
MAP_ERROR_CODE reset_force_to_zero(double* fx, double* fy, double* fz, double* mx, double* my, double* mz, const int size);
MAP_ERROR_CODE set_force_minus(const double* inputType, double* force, const int size);
MAP_ERROR_CODE set_force_plus(const double* inputType, double* force, const int size);
MAP_ERROR_CODE calculate_stiffness(double* K, double* force, const double delta, const int size);
MAP_ERROR_CODE set_moment_minus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size);
MAP_ERROR_CODE set_moment_plus(const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size);
MAP_ERROR_CODE increment_phi_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size);
MAP_ERROR_CODE increment_the_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size);
MAP_ERROR_CODE increment_psi_dof_by_delta(MAP_InputType_t* uType, const Vessel* vessel, const double delta, const int size);
MAP_ERROR_CODE fd_x_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_y_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_z_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalPos, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_phi_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_the_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE fd_psi_sequence(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, MAP_OutputType_t* yType, MAP_ConstraintStateType_t* zType, Fd* force, const double epsilon, const int size, const double* originalX, const double* originalY, const double* originalZ, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE calculate_stiffness_2(double* K, Fd* force, const double delta, const int size);


MAP_ERROR_CODE set_moment_plus_2(const MAP_InputType_t* uType, const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size);
MAP_ERROR_CODE set_moment_minus_2(const MAP_InputType_t* uType, const MAP_OutputType_t* outputType, const Vessel* vessel, double* mx, double* my, double* mz, const int size);

#endif /* _PROTOS_H */
