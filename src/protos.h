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

MAP_EXTERNCALL void map_offset_vessel(MAP_OtherStateType_t* otherType, MAP_InputType_t* uType, double x, double y, double z, double phi, double the, double psi, char* map_msg, MAP_ERROR_CODE* ierr);
void end_color( char* dest, const char* src );
void map_reset_universal_error( char* user_str, MAP_ERROR_CODE* ierr );

/**
 * @brief   Set error code with customized message
 * @details Can be called in any routine to set the contents of map_message. The routine
 *          checks to ensure map_msg is not longer than MAX_ERROR_STRING_LENGTH. If it is,
 *          then the message is truncated. The error code returned depends on the max 
 *          level reached previously durring program execusion. A MAP_WARNING will not 
 *          override a MAP_FATAL or MAP_ERROR level. 
 * @param   user_string, custom message appending the error message
 * @param   map_msg, error message of length no greater than MAP_ERROR_STRING_LENGTH
 * @param   ierr, ierr status before function is invoked, should be MAP_SAFE, MAP_WARNING, MAP_ERROR, MAP_FATAL
 * @param   new_error, new MAP warning code; any entry in enum MAP_ERROR_CODE
 * @return  MAP error code
 */
MAP_ERROR_CODE map_set_universal_error(bstring user_msg, char* map_msg, const MAP_ERROR_CODE ierr, const MAP_ERROR_CODE new_error);


/*
 *
 */
MAP_ERROR_CODE set_element_vartype_ptr(char* unit, char* alias, const int num, VarTypePtr* type);

const char* remove_first_character(const char* str);
MAP_ERROR_CODE initialize_external_applied_force(char* unit, char* alias, const int num, VarType* type, char const* property);
MAP_ERROR_CODE initialize_node_sum_force_ptr(char* unit, char* alias, const int num, VarTypePtr* type);

size_t vartype_meter(const void* el);
size_t vartype_ptr_meter(const void *el);
MAP_ERROR_CODE free_outer_solve_data(OuterSolveAttributes* ns, const int size, char* map_msg, MAP_ERROR_CODE* ierr);
void print_machine_name_to_screen( );

/*
 * stuff for nodes...
 */

/**
 * Numeric routines
 */
int inner_function_evals(void* elementPtr, int m, int n, const __cminpack_real__* x, __cminpack_real__* fvec, __cminpack_real__* fjac, int ldfjac, int iflag);






/*
 * MARK: numeric helper functions
 */

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
