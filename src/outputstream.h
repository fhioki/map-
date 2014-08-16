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


#ifndef _OUTPUT_STREAM_H
#define _OUTPUT_STREAM_H


#include <time.h>
#include "map.h"
#include "maperror.h"
#include "MAP_Types.h"


/**
 * @brief Set the name out the MAP summary output file. Does not need to be called; the default
 *        summary file name is 'outlist.map.sum'.
 * @param init_type initalization type, native C struct {@link InitializationData_t}
 * @param map_msg MAP error message
 * @param ierr MAP error code
 *
 * Example Fortran usage:
 * @code
 * ! Interface block declaration:
 * INTERFACE 
 *    SUBROUTINE mapextern_map_set_summary_file_name(fc_init, fc_char, fc_int) BIND(C,name='map_set_summary_file_name')
 *      IMPORT                                    
 *      IMPLICIT NONE                             
 *      TYPE(MAP_InitInputType_C) fc_init        
 *      CHARACTER(KIND=C_CHAR), DIMENSION(*) :: fc_char
 *      INTEGER(KIND=C_INT) :: fc_int            
 *    END SUBROUTINE mapextern_map_set_summary_file_name
 * END INTERFACE
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_set_summary_file_name(InitInp%C_obj, ErrMsg, ErrStat)
 * @endcode
 * @todo: need to free summaryFileName. This is done in delete_all_init_data(...), should be called in Fortran routines
 */
MAP_EXTERNCALL void map_set_summary_file_name(MAP_InitInputType_t* init_type, char* map_msg, MAP_ERROR_CODE* ierr); 


/**
 * @brief Obtains the variable name array corresponding to the outputs selected in the MAP input file. For example,
 *        str_array can be:
 *        <pre>
 *        X[2]     H[1]     X[6]     H[3]     X[10]    H[5]     X[14]    H[7]
 *        </pre>
 * @param n number of header blocks. Should be proportional to the number of itms being output to the FAST output file
 * @param str_array the string being output.         
 * @param other_type Fortran other state derived type
 *
 * Example Fortran usage:
 * @code
 * ! Interface block declaration:
 * INTERFACE      
 *    SUBROUTINE mapextern_map_get_header_string(fc_int, fc_string, fc_other) BIND(C,name='map_get_header_string')   
 *      IMPORT                                 
 *      IMPLICIT NONE                          
 *      INTEGER(KIND=C_INT) :: fc_int          
 *      TYPE( MAP_OtherStateType_C ) fc_other                                                 
 *      TYPE(C_PTR), DIMENSION(FC_int) :: fc_string
 *    END SUBROUTINE mapextern_map_get_header_string
 * END INTERFACE                                                 
 *
 *   ...
 *
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_get_header_string(num_header_str, hdr_str_ptrs, other%C_obj)
 * @endcode
 * @todo this should raise and error when count!=n
 */
MAP_EXTERNCALL void map_get_header_string(int* n, char** str_array, MAP_OtherStateType_t* other_type);


/**
 * @brief Obtains the units of the outputs passed back to the calling program. str_array can be:
 *        <pre>
 *        [m]     [N]     [m]     [N]     [m]     [N]     [m]     [N]   
 *        </pre>
 * @param n number of header blocks. Should be proportional to the number of itms being output to the FAST output file
 * @param str_array the string being output.         
 * @param other_type Fortran other state derived type
 *
 * Example Fortran usage:
 * @code
 * ! Interface block declaration:
 * INTERFACE
 *    SUBROUTINE mapextern_map_get_unit_string(fc_int, fc_string, fc_other) BIND(C,name='map_get_unit_string')          
 *      IMPORT                                     
 *      IMPLICIT NONE                              
 *      INTEGER(KIND=C_INT) :: fc_int              
 *      TYPE(MAP_OtherStateType_C) fc_other                                                 
 *      TYPE(C_PTR), DIMENSION(FC_int) :: fc_string
 *    END SUBROUTINE mapextern_map_get_unit_string        
 * END INTERFACE                                                 
 *
 *   ...
 * ! access the function using this subroutine call: 
 * CALL mapextern_map_get_header_string(num_header_str, unit_str_ptrs, other%C_obj)
 * @endcode
 * @todo this should raise and error when count!=n
 */
MAP_EXTERNCALL void map_get_unit_string(int* n, char** str_array ,MAP_OtherStateType_t* other_type);


/**
 * @brief
 * @param y_type output type, native C struct {@link MAP_OutputType_t}
 * @param other_type other state type, native C struct {@link MAP_OtherStateType_t}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE get_iteration_output_stream(MAP_OutputType_t *y_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief
 * @param init_data MAP internal initialization data structure
 * @param param_type parmeter type, native C struct {@link MAP_ParameterType_t}
 * @param model_data internal state data {@link ModelData}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE write_summary_file(InitializationData* init_data, MAP_ParameterType_t* param_type, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 * @brief Writes all cable properties defined in the MAP input file:
 *        <pre>
 *        Cable Type          : {value}
 *        Diameter     [m]    : {value}
 *        Mass Density [kg/m] : {value}
 *        EA           [N]    : {value}
 *        omega        [N/m]  : {value}
 *        CB                  : {value}
 *        </pre>
 * @param file file where string is being dumped to
 * @param model_data internal state data {@link ModelData}
 * @todo  include new properties for the LM model
 */
MAP_ERROR_CODE write_cable_library_information_to_summary_file(FILE* file, ModelData* model_data);

/**
 * @brief writes the node type. Can be VESSEL, CONNECT, or FIX
 *        <pre>
 *        Type      |  {VESSEL} 
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param node_type identifies the node type as a NONE, CONNECT, FIX, or VESSEL
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_type_to_summary_file(const int num_col, const int count_to_four, const NodeType node_type, char* line_char);

/**
 * @brief header for the node block. Prints the node number. 
 *        <pre>
 *                  | Node {number} Data"
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param node_num
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_header_to_summary_file(const int num_col, const int count_to_four, const int node_num, char* line_char);

/**
 * @brief writes the node x global displacement (position) to the node output block. Units are in [m]: 
 *        <pre>
 *        X  [m]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param x_pos node global x position [m]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_x_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_pos, char* line_char);

/**
 * @brief writes the node y global displacement (position) to the node output block. Units are in [m]: 
 *        <pre>
 *        Y  [m]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param y_pos node global y position [m]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_y_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_pos, char* line_char);

/**
 * @brief writes the node z global displacement (position) to the node output block. Units are in [m]: 
 *        <pre>
 *        Z  [m]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param z_pos node global z position [m]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_z_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_pos, char* line_char);

/**
 * @brief writes the node point mass value to the output block. Units are in [kg]: 
 *        <pre>
 *        M  [kg]  | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param point_mass node point mass [kg]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_mass_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_mass, char* line_char);

/**
 * @brief writes the volumetric displacement of the buoyancy module to the node output block. Units are in [m^3]: 
 *        <pre>
 *        B  [m^3]  | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param point_buoyancy node point buoyancy [m^3]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_buoyancy_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_buoy, char* line_char);

/**
 * @brief writes the x-direction sum force to the node output block. Units are in [N]: 
 *        <pre>
 *        FX [N]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param x_sum_force node global x sum force (including external forces) [N]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_x_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_sum_force, char* line_char);

/**
 * @brief writes the y-direction sum force to the node output block. Units are in [N]: 
 *        <pre>
 *        FY [N]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param y_sum_force node global y sum force (including external forces) [N]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_y_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_sum_force, char* line_char);

/**
 * @brief writes the z-direction sum force to the node output block. Units are in [N]: 
 *        <pre>
 *        FZ [N]   | {value}
 *        </pre>
 * @param num_col column number, can be 1, 2, 3, or 4
 * @param count_to_four number of columns currently printed in the node output block 
 * @param z_sum_force node global z sum force (including external forces) [N]
 * @param line_char output string terminated with '\n'
 */
MAP_ERROR_CODE write_node_z_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_sum_force, char* line_char);

/**
 * @brief write the complete node block to the summary file
 *        <pre>
 *                  | Node 1 Data            Node 2 Data            Node 3 Data            Node 4 Data
 *                  | -------------------------------------------------------------------------------------------
 *        Type      |  FIX                    CONNECT                VESSEL                 FIX
 *        X  [m]    |   0.000                  0.000                  0.000                  0.000
 *        Y  [m]    |   0.000                  0.000                  0.000                  0.000
 *        Z  [m]    |   0.000                  0.000                  0.000                  0.000
 *        M  [kg]   |   0.000                  0.000                  0.000                  0.000
 *        B  [m^3]  |   0.000                  0.000                  0.000                  0.000
 *        FX [N]    |   0.000                  0.000                  0.000                  0.000
 *        FY [N]    |   0.000                  0.000                  0.000                  0.000
 *        FZ [N]    |   0.000                  0.000                  0.000                  0.000
 *        </pre>
 * @param file file where string is being dumped to
 * @param model_data internal state data {@link ModelData}
 * @param map_msg error message
 * @param ierr error code
 */
MAP_ERROR_CODE write_node_information_to_summary_file(FILE* file, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr);

/**
 * @brief write the complete element block to the summary file
 *        <pre>
 *                        | Element 1
 *                        | ---------------------------------------
 *        Material        |  Material
 *        Lu        [m]   |  0.000
 *        Lb        [m]   |  0.000
 *        H         [N]   |  0.000
 *        V         [N]   |  0.000
 *        T         [N]   |  0.000
 *        Alpha     [deg] |  0.000
 *        HAnch     [N]   |  0.000
 *        VAnch     [N]   |  0.000
 *        TAnch     [N]   |  0.000
 *        AlphaAnch [deg] |  0.000
 *        L^2-Norm        |  0.000
 *        Function Evals  |  0
 *        Jacobian Evals  |  0
 *        Term. criteria  |  0
 *        </pre>
 * @param file file where string is being dumped to
 * @param model_data internal state data {@link ModelData}
 */
MAP_ERROR_CODE write_element_information_to_summary_file(FILE* file, ModelData* model_data);

/**
 * @brief prints the expanded MAP input file. This can be used as a check to make sure the repeat flags 
 *        are correctly interpreted. 
 * @param file file where string is being dumped to
 * @param init_data initialization output type, native C struct {@link MAP_InitOutputType_t}
 */
MAP_ERROR_CODE write_expanded_input_file_to_summary_file(FILE* file, InitializationData* init_data);


#endif // _OUTPUT_STREAM_H
