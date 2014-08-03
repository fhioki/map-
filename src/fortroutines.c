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
#include "protos.h"
#include "pyprotos.h"
#include "MAP_Types.h"


MAP_EXTERNCALL void map_add_cable_library_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = (InitializationData*)init_type->object; 
  const int str_length = strlen(init_type->libraryInputLine)+1;
  const int n = init_data->librarySize;
  char** temp_string_array = NULL;
  
  temp_string_array = realloc(init_data->libraryInputString, sizeof(*init_data->libraryInputString)*(n+1));
  init_data->libraryInputString = temp_string_array;
  init_data->libraryInputString[n] = malloc(sizeof(char)*str_length);    
  strcpy(init_data->libraryInputString[n], init_type->libraryInputLine);  
  init_data->librarySize++;
};


MAP_EXTERNCALL void map_add_node_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int str_length = strlen(init_type->nodeInputLine)+1;
  const int n = init_data->nodeSize;
  char** temp_string_array = NULL;

  temp_string_array = realloc(init_data->nodeInputString, sizeof(*init_data->nodeInputString)*(n+1));
  init_data->nodeInputString = temp_string_array;
  init_data->nodeInputString[n] = malloc(sizeof(char)*str_length);    
  strcpy(init_data->nodeInputString[n], init_type->nodeInputLine);
  init_data->nodeSize++;
};


MAP_EXTERNCALL void map_add_element_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  const int n = init_data->elementSize;
  const int str_length = strlen(init_type->elementInputLine)+1;
  char** temp_string_array = NULL;

  temp_string_array = realloc(init_data->elementInputString, sizeof(*init_data->elementInputString)*(n+1));
  init_data->elementInputString = temp_string_array;
  init_data->elementInputString[n] = malloc(sizeof(char)*str_length);    
  strcpy(init_data->elementInputString[n], init_type->elementInputLine);
  init_data->elementSize++;
};


MAP_EXTERNCALL void map_add_options_input_text(MAP_InitInputType_t* init_type)
{
  InitializationData* init_data = init_type->object; 
  int n = init_data->solverOptionsSize;
  int str_length = strlen(init_type->optionInputLine)+1;
  char** temp_string_array = NULL;

  temp_string_array = realloc(init_data->solverOptionsString, sizeof(*init_data->solverOptionsString)*(n+1));
  init_data->solverOptionsString = temp_string_array;
  init_data->solverOptionsString[n] = malloc(sizeof(char)*str_length);    
  strcpy(init_data->solverOptionsString[n], init_type->optionInputLine);
  init_data->solverOptionsSize++;
};
