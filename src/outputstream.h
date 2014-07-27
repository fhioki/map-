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

#include "MAP_Types.h"

/**
 *
 */
MAP_EXTERNCALL void get_header_string(int* N, char** strArr, MAP_OtherStateType_t* otherType);


/**
 *
 */
MAP_EXTERNCALL void get_unit_string(int* N, char** strArr ,MAP_OtherStateType_t* otherType);


/**
 *
 */
MAP_ERROR_CODE get_iteration_output_stream(MAP_OutputType_t *yType, MAP_OtherStateType_t* otherType, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 *
 */
MAP_ERROR_CODE set_output_list(ModelData* data, MAP_InitOutputType_t* ioType, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 *
 */
MAP_ERROR_CODE write_summary_file(InitializationData* init, MAP_ParameterType_t* paramType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 *
 */
MAP_ERROR_CODE write_cable_library_information_to_summary_file(FILE* file, ModelData* dataObj);


/**
 *
 */
MAP_ERROR_CODE write_node_type_to_summary_file(const int columnNumber, const int countToFour, const NodeType nodeType, char* line);

/**
 *
 */
MAP_ERROR_CODE write_node_header_to_summary_file(const int columnNumber, const int countToFour, const int nodeNumber, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_x_position_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* xPosition, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_y_position_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* yPosition, char* line);

/**
 *
 */
MAP_ERROR_CODE write_node_z_position_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* zPosition, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_mass_information_to_summary_file(const int columnNumber, const int countToFour, VarType* nodePointMass, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_buoyancy_information_to_summary_file(const int columnNumber, const int countToFour, VarType* nodeBuoyancyModule, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_x_sum_force_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* xSumForce, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_y_sum_force_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* ySumForce, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_z_sum_force_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* zSumForce, char* line);


/**
 *
 */
MAP_ERROR_CODE write_node_information_to_summary_file(FILE* file, ModelData* dataObj, char* map_msg, MAP_ERROR_CODE* ierr);


/**
 *
 */
MAP_ERROR_CODE write_element_information_to_summary_file(FILE* file, ModelData* dataObj);


/**
 *
 */
MAP_ERROR_CODE write_expanded_input_file_to_summary_file(FILE* file, InitializationData* initData);


#endif // _OUTPUT_STREAM_H
