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


/**
 * @file 
 */


#include "outputstream.h"


/** @addtogroup FortranCall */
/* @{ */
MAP_EXTERNCALL void map_set_summary_file_name(MAP_InitInputType_t* init_type, char *map_msg, MAP_ERROR_CODE *ierr) 
{  
  InitializationData* init_data = init_type->object; 
  
  if (init_data->summaryFileName!=NULL) {
    MAPFREE(init_data->summaryFileName);
  };

  init_data->summaryFileName = bformat("%s", init_type->summaryFileName);
};


// MAP_EXTERNCALL void map_get_header_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
// { 
//   int count = 0;    
//   ModelData* model_data = other_type->object;
//   VarTypePtr* vartype_ptr = NULL;
//   VarType* vartype = NULL;
// 
//   list_iterator_start(&model_data->yList->out_list_ptr);
//   while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
//     vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr);
//     strcpy(str_array[count],vartype_ptr->name);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list_ptr);     
// 
//   list_iterator_start(&model_data->yList->out_list);
//   while (list_iterator_hasnext(&model_data->yList->out_list)) { 
//     vartype = (VarType*)list_iterator_next(&model_data->yList->out_list);
//     strcpy(str_array[count],vartype->name);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list);     
//   /* @todo this should raise and error when count != n */
// };


// MAP_EXTERNCALL void map_get_unit_string(int* n, char** str_array, MAP_OtherStateType_t* other_type)
// { 
//   int count = 0;    
//   ModelData* model_data = other_type->object;
//   VarTypePtr* vartype_ptr = NULL;
//   VarType* vartype = NULL;
// 
//   list_iterator_start(&model_data->yList->out_list_ptr);
//   while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
//     vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr );
//     strcpy(str_array[count],vartype_ptr->units);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list_ptr);     
// 
//   list_iterator_start(&model_data->yList->out_list);
//   while (list_iterator_hasnext(&model_data->yList->out_list)) { 
//     vartype = (VarType*)list_iterator_next(&model_data->yList->out_list );
//     strcpy(str_array[count],vartype->units);
//     count++;
//   };
//   list_iterator_stop(&model_data->yList->out_list);     
// };
/* @} */


MAP_ERROR_CODE get_iteration_output_stream(MAP_OutputType_t *y_type, MAP_OtherStateType_t* other_type, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int count = 0;
  ModelData* model_data = other_type->object;
  VarTypePtr* vartype_ptr = NULL;
  VarType* vartype = NULL;

  if (!y_type->wrtOutput) { /* if NULL, then wrtOutput is not initialized */
    /* first set size of out list */
    y_type->wrtOutput_Len = (int)list_size(&model_data->yList->out_list_ptr);
    y_type->wrtOutput_Len += (int)list_size(&model_data->yList->out_list);
    
    /* allocation array */
    y_type->wrtOutput = malloc(sizeof(double)*y_type->wrtOutput_Len);
  };

  list_iterator_start(&model_data->yList->out_list_ptr);
  while (list_iterator_hasnext(&model_data->yList->out_list_ptr)) { 
    vartype_ptr = (VarTypePtr*)list_iterator_next(&model_data->yList->out_list_ptr);
    y_type->wrtOutput[count] = *vartype_ptr->value;
    count++;
  };
  list_iterator_stop(&model_data->yList->out_list_ptr);     

  list_iterator_start(&model_data->yList->out_list);
  while (list_iterator_hasnext(&model_data->yList->out_list)) { 
    vartype = (VarType*)list_iterator_next(&model_data->yList->out_list);
    y_type->wrtOutput[count] = vartype->value;
    count++;
  };
  list_iterator_stop(&model_data->yList->out_list);     
  return MAP_SAFE;
};


MAP_ERROR_CODE write_summary_file(InitializationData* init, MAP_ParameterType_t* paramType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  struct tm* tm_info;
  time_t timer;
  FILE* file = fopen(init->summaryFileName->data, "w");
  char time_buffer[TIME_BUFFER_SIZE] = "\0";
  bstring user_msg = NULL;
  
  time(&timer);
  tm_info = localtime(&timer);
  
  if (file==NULL) {
    *ierr = map_set_universal_error(init->summaryFileName, map_msg, *ierr, MAP_FATAL_38);
    return MAP_FATAL;
  };
  
  //__get_machine_name(name);  
  strftime(time_buffer, TIME_BUFFER_SIZE, "%A %B %d-%Y at %H:%M:%S %p", tm_info);

  /* May want to retag version number of git between major/minor releases. A trick to include git revision number
   * can be found archived here: http://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
   *
   * $ git tag -a v1.00 -m "Release v1.00"
   */

  /* unboubtedly this will be useful on windows. Not sure how to track git version 
   * on windows with MSVC compiler, or for those who don't use version control? */
#  ifndef GITVERSION 
#    define GITVERSION PROGVERSION
#  endif

  fprintf(file,"Outputs were generated using MAP++ %s.\n\n", GITVERSION); 
  fprintf(file,"    %s\n\n", time_buffer); 
  fprintf(file,"    Gravity constant          [m/s^2]  : %1.2f\n", paramType->g); 
  fprintf(file,"    Sea density               [kg/m^3] : %1.2f\n", paramType->rhoSea);
  fprintf(file,"    Water depth               [m]      : %1.2f\n", paramType->depth);
  fprintf(file,"    Vessel reference position [m]      : %1.2f , %1.2f , %1.2f\n\n", data->vessel.refOrigin.x.value, data->vessel.refOrigin.y.value, data->vessel.refOrigin.z.value); 
  
  success = write_cable_library_information_to_summary_file(file, data); CHECKERRK(MAP_ERROR_6);
  success = write_node_information_to_summary_file(file, data, map_msg, ierr); CHECKERRK(MAP_ERROR_7);
  success = write_element_information_to_summary_file(file, data); CHECKERRK(MAP_ERROR_8);
  success = write_expanded_input_file_to_summary_file(file, init);
  fclose(file);  
  
  MAP_RETURN;
};


MAP_ERROR_CODE write_cable_library_information_to_summary_file(FILE* file, ModelData* model_data)
{
  CableLibrary* library_iter = NULL;  

  list_iterator_start(&model_data->cableLibrary);    
  while (list_iterator_hasnext(&model_data->cableLibrary)) { 
    library_iter = (CableLibrary*)list_iterator_next(&model_data->cableLibrary);    
    fprintf(file, "    Cable Type          : %s\n", library_iter->label->data);
    fprintf(file, "    Diameter     [m]    : %1.4f\n", library_iter->diam);
    fprintf(file, "    Mass Density [kg/m] : %1.2f\n", library_iter->massDensityInAir);
    fprintf(file, "    EA           [N]    : %1.2f\n", library_iter->ea);
    fprintf(file, "    omega        [N/m]  : %1.2f\n", library_iter->omega);
    fprintf(file, "    CB                  : %1.2f\n\n", library_iter->cb);
  };
  list_iterator_stop(&model_data->cableLibrary);  
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_type_to_summary_file(const int num_col, const int count_to_four, const NodeType node_type, char* line_char)
{
  int j = 0;
  int size = 0;

  if (!num_col) { /* if this is the firs column, then line the right side with information */
    if (node_type==VESSEL) {
      map_snprintf(line_char, 256, "Type      |  VESSEL");
    } else if (node_type==FIX ) { 
      map_snprintf(line_char, 256, "Type      |  FIX");
    } else if (node_type==CONNECT ) {
      map_snprintf(line_char, 256, "Type      |  CONNECT");
    } else {
      map_snprintf(line_char, 256, "Type      |  UNKNOWN");
    };
  } else { /* otherwise, just print the information */
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) { /* add white spaces */
      map_strcat(line_char, 256, " "); 
    };

    if (node_type==VESSEL) {
      map_strcat(line_char, 256, " VESSEL");        
    } else if (node_type==FIX ) { 
      map_strcat(line_char, 256, " FIX");        
    } else if(node_type==CONNECT ) {
      map_strcat(line_char, 256, " CONNECT");        
    } else {
      map_strcat(line_char, 256, " UNKNOWN");        
    };
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_header_to_summary_file(const int num_col, const int count_to_four, const int node_num, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";

  if (!num_col) { 
    map_snprintf(line_char, 256, "          | Node %d Data", node_num);        
  } else {
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " ");       
    };
    map_snprintf(buffer, 64, "%sNode %d Data", spaces, node_num);
    map_strcat(line_char, 256, buffer);       
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_x_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_pos, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal x = *x_pos->value;

  if (!num_col) { 
    if (!x_pos->isFixed) {
      if (x>=0.0) {
        map_snprintf(line_char, 256, "X  [m]    | ( %1.3f)", x);
      } else {
        map_snprintf(line_char, 256, "X  [m]    | (%1.3f)", x);        
      };
    } else {
      if (x>0.0) {
        map_snprintf(line_char, 256, "X  [m]    |   %1.3f", x);        
      } else {
        map_snprintf(line_char, 256, "X  [m]    |  %1.3f", x);        
      };
    };
  } else { 
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!x_pos->isFixed) {
      if (x>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, x);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, x);        
      };
    } else {
      if (x>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, x);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, x);        
      };
    };
    map_strcat(line_char, 256, buffer);       
  };

  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_y_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_pos, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal y = *y_pos->value;

  if (!num_col) { 
    if (!y_pos->isFixed) {
      if (y>=0.0) {
        map_snprintf(line_char, 256, "Y  [m]    | ( %1.3f)", y);
      } else {
        map_snprintf(line_char, 256, "Y  [m]    | (%1.3f)", y);        
      };
    } else {
      if (y>=0.0) {
        map_snprintf(line_char, 256, "Y  [m]    |   %1.3f", y);        
      } else {
        map_snprintf(line_char, 256, "Y  [m]    |  %1.3f", y);        
      };
    };
  } else { 
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };
    if (!y_pos->isFixed) {
      if (y>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, y);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, y);        
      };
    } else {
      if (y>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, y);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, y);        
      };
    };
    map_strcat(line_char, 256, buffer);       
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_z_position_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_pos, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal z = *z_pos->value;
  
  if (!num_col) { 
    if (!z_pos->isFixed) {
      if (z>=0.0) {
        map_snprintf(line_char, 256, "Z  [m]    | ( %1.3f)", z);
      } else {
        map_snprintf(line_char, 256, "Z  [m]    | (%1.3f)", z);        
      };
    } else {
      if (z>=0.0) {
        map_snprintf(line_char, 256, "Z  [m]    |   %1.3f", z);        
      } else {
        map_snprintf(line_char, 256, "Z  [m]    |  %1.3f", z);        
      };
    };
  } else { 
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!z_pos->isFixed) {
      if (z>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, z);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, z);        
      };
    } else {
      if (z>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, z);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, z);        
      };
    };
    map_strcat(line_char, 256, buffer);       
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_mass_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_mass, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal mass = point_mass->value;

  if (!num_col) { 
    if (!point_mass->isFixed) {
      if (mass>=0.0) {
        map_snprintf(line_char, 256, "M  [kg]   | ( %1.3f)", mass);
      } else {
        map_snprintf(line_char, 256, "M  [kg]   | (%1.3f)", mass);        
      };
    } else {
      if (mass>=0.0) {
        map_snprintf(line_char, 256, "M  [kg]   |   %1.3f", mass);        
      } else {
        map_snprintf(line_char, 256, "M  [kg]   |  %1.3f", mass);        
      };
    };

  } else { 
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!point_mass->isFixed) {
      if (mass>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, mass);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, mass);        
      };
    } else {
      if (mass>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, mass);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, mass);        
      };
    };
    map_strcat(line_char, 256, buffer);       
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_buoyancy_information_to_summary_file(const int num_col, const int count_to_four, VarType* point_buoy, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal buoyancy = point_buoy->value;

  if (!num_col) { 
    if (!point_buoy->isFixed) {
      if (buoyancy>=0.0) {
        map_snprintf(line_char, 256, "B  [m^3]  | ( %1.3f)", buoyancy);
      } else {
        map_snprintf(line_char, 256, "B  [m^3]  | (%1.3f)", buoyancy);        
      };
    } else {
      if (buoyancy>=0.0) {
        map_snprintf(line_char, 256, "B  [m^3]  |   %1.3f", buoyancy);        
      } else {
        map_snprintf(line_char, 256, "B  [m^3]  |  %1.3f", buoyancy);        
      };
    };
  } else { 
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };
    
    if (!point_buoy->isFixed) {
      if (buoyancy>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, buoyancy);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, buoyancy);        
      };
    } else {
      if (buoyancy>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, buoyancy);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, buoyancy);        
      };
    };
    map_strcat(line_char, 256, buffer);       
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_x_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* x_sum_force, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal fx = *(x_sum_force->value);

  if (!num_col) { 
    if (!x_sum_force->isFixed) {
      if (fx>=0.0) {
        map_snprintf(line_char, 256, "FX [N]    | ( %1.3f)", fx);
      } else {
        map_snprintf(line_char, 256, "FX [N]    | (%1.3f)", fx);        
      };
    } else {
      if (fx>=0.0) {
        map_snprintf(line_char, 256, "FX [N]    |   %1.3f", fx);        
      } else {
        map_snprintf(line_char, 256, "FX [N]    |  %1.3f", fx);        
      };
    };
  } else {
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    }
    if (!x_sum_force->isFixed) {
      if (fx>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, fx);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, fx);        
      };
    } else {
      if (fx>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, fx);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, fx);        
      };
    };
    map_strcat(line_char, 256, buffer);
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_y_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* y_sum_force, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal fy = *(y_sum_force->value);

  if (!num_col) { 
    if (!y_sum_force->isFixed) {
      if (fy>=0.0) {
        map_snprintf(line_char, 256, "FY [N]    | ( %1.3f)", fy);
      } else {
        map_snprintf(line_char, 256, "FY [N]    | (%1.3f)", fy);        
      };
    } else {
      if (fy>=0.0) {
        map_snprintf(line_char, 256, "FY [N]    |   %1.3f", fy);        
      } else {
        map_snprintf(line_char, 256, "FY [N]    |  %1.3f", fy);        
      };
    };
  } else {
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    }

    if (!y_sum_force->isFixed) {
      if (fy>=0.0) {
        map_snprintf(buffer, 64, "%s( %1.3f)", spaces, fy);
      } else {
        map_snprintf(buffer, 64, "%s(%1.3f)", spaces, fy);        
      };
    } else {
      if (fy>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, fy);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, fy);        
      };
    };
    map_strcat(line_char, 256, buffer);
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_z_sum_force_to_summary_file(const int num_col, const int count_to_four, VarTypePtr* z_sum_force, char* line_char)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal fz = *(z_sum_force->value);

  if (!num_col) { 
    if (!z_sum_force->isFixed) {
      if (fz>=0.0) {
        map_snprintf(line_char, 256, "FZ [N]    | ( %1.3f)", fz);
      } else {
        map_snprintf(line_char, 256, "FZ [N]    | (%1.3f)", fz);        
      };
    } else {
      if (fz>=0.0) {
        map_snprintf(line_char, 256, "FZ [N]    |   %1.3f", fz);        
      } else {
        map_snprintf(line_char, 256, "FZ [N]    |  %1.3f", fz);        
      };
    };
  } else {
    size = (int)strlen(line_char)-23*count_to_four; /* each node column is 23 character wide */ 
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    }
    if (!z_sum_force->isFixed) {
      if (fz>=0.0) {
        map_snprintf( buffer, 64, "%s( %1.3f)", spaces, fz);
      } else {
        map_snprintf( buffer, 64, "%s(%1.3f)", spaces, fz);        
      };
    } else {
      if (fz>=0.0) {
        map_snprintf(buffer, 64, "%s  %1.3f", spaces, fz);        
      } else {
        map_snprintf(buffer, 64, "%s %1.3f", spaces, fz);        
      };
    };
    map_strcat(line_char, 256, buffer);
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_node_information_to_summary_file(FILE* file, ModelData* model_data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int num = 0;
  int i = 0;
  int col_cnt = 0; /* number of columns. The MAP summary file prints 4 columns of node information. After 
                    * column 4, carriage return */
  char line0[256] = "";
  char line1[256] = "";
  char line2[256] = "";
  char line3[256] = "";
  char line4[256] = "";
  char line5[256] = "";
  char line6[256] = "";
  char line7[256] = "";
  char line8[256] = "";
  char line9[256] = "";
  Node* node_iter = NULL;  
  const int FOUR = 4;
  const unsigned int num_nodes = list_size(&model_data->node);  
  unsigned int col = 0;
  MAP_ERROR_CODE success = MAP_SAFE;
  bstring user_msg = NULL;

  do {
    for (i=0 ; i<num_nodes ; i+=FOUR) {
      if (i+FOUR>num_nodes) {
        num = num_nodes-i;
      } else {
        num = FOUR;
      };    
      for (col=i ; col<i+num ; col++) {
        node_iter = (Node*)list_get_at(&model_data->node, col);      
        success = write_node_header_to_summary_file(col-i, col_cnt, col+1, line0); CHECKERRQ(MAP_FATAL_70);
        success = write_node_type_to_summary_file(col-i, col_cnt, node_iter->type, line1); CHECKERRQ(MAP_FATAL_70);
        success = write_node_x_position_to_summary_file(col-i, col_cnt, &node_iter->positionPtr.x, line2); CHECKERRQ(MAP_FATAL_70);
        success = write_node_y_position_to_summary_file(col-i, col_cnt, &node_iter->positionPtr.y, line3); CHECKERRQ(MAP_FATAL_70);
        success = write_node_z_position_to_summary_file(col-i, col_cnt, &node_iter->positionPtr.z, line4); CHECKERRQ(MAP_FATAL_70);
        success = write_node_mass_information_to_summary_file(col-i, col_cnt, &node_iter->MApplied, line5); CHECKERRQ(MAP_FATAL_70);
        success = write_node_buoyancy_information_to_summary_file(col-i, col_cnt, &node_iter->BApplied, line6); CHECKERRQ(MAP_FATAL_70);
        success = write_node_x_sum_force_to_summary_file(col-i, col_cnt, &node_iter->sumForcePtr.fx, line7); CHECKERRQ(MAP_FATAL_70);
        success = write_node_y_sum_force_to_summary_file(col-i, col_cnt, &node_iter->sumForcePtr.fy, line8); CHECKERRQ(MAP_FATAL_70);
        success = write_node_z_sum_force_to_summary_file(col-i, col_cnt, &node_iter->sumForcePtr.fz, line9); CHECKERRQ(MAP_FATAL_70);
        col_cnt++;
      };
      col_cnt = 0;
      map_strcat(line0, 256, "\n");
      map_strcat(line1, 256, "\n");
      map_strcat(line2, 256, "\n");
      map_strcat(line3, 256, "\n");
      map_strcat(line4, 256, "\n");
      map_strcat(line5, 256, "\n");
      map_strcat(line6, 256, "\n");
      map_strcat(line7, 256, "\n");
      map_strcat(line8, 256, "\n");
      map_strcat(line9, 256, "\n");
    
      fprintf(file, "%s",line0);
      fprintf(file, "          | -------------------------------------------------------------------------------------------\n");
      fprintf(file, "%s",line1);
      fprintf(file, "%s",line2);
      fprintf(file, "%s",line3);
      fprintf(file, "%s",line4);
      fprintf(file, "%s",line5);
      fprintf(file, "%s",line6);
      fprintf(file, "%s",line7);
      fprintf(file, "%s",line8);
      fprintf(file, "%s\n\n",line9);
    
      line0[0] = 0; 
      line1[0] = 0; 
      line2[0] = 0; 
      line3[0] = 0; 
      line4[0] = 0; 
      line5[0] = 0; 
      line6[0] = 0; 
      line7[0] = 0; 
      line8[0] = 0; 
      line9[0] = 0; 
    };
  } while (0);
  MAP_RETURN;
};


MAP_ERROR_CODE write_element_information_to_summary_file(FILE* file, ModelData* model_data)
{
  const unsigned int num_elements = list_size(&model_data->element);
  bstring user_msg = NULL;
  Element* element_iter = NULL;  
  char line0[256] = "";
  char line1[256] = "";
  char line2[256] = "";
  char line3[256] = "";
  char line4[256] = "";
  char line5[256] = "";
  char line6[256] = "";
  char line7[256] = "";
  char line8[256] = "";
  char line9[256] = "";
  char line10[256] = "";
  int i = 0;  

  for (i=0 ; i<num_elements ; i++) {
    element_iter = (Element*)list_get_at(&model_data->element, i);
    
    if (element_iter->Lu.value>0.0) {
      map_strcat(line0, 256, " ");
    };
    if (element_iter->lb.value>0.0) {
      map_strcat(line1, 256, " ");
    };
    if (*(element_iter->H.value)>0.0) {
        map_strcat(line2, 256, " ");
    };
    if (*(element_iter->V.value)>0.0) {
      map_strcat(line3, 256, " ");
    };
    if (element_iter->T.value>0.0) {
      map_strcat(line4, 256, " ");
    };
    if (element_iter->alpha.value>0.0) {
      map_strcat(line5, 256, " ");
    };
    if (element_iter->HAtAnchor.value>0.0) {
      map_strcat(line6, 256, " ");
    };
    if (element_iter->VAtAnchor.value>0.0) {
      map_strcat(line7, 256, " ");
    };
    if (element_iter->TAtAnchor.value>0.0) {
      map_strcat(line8, 256, " ");
    };
    if (element_iter->alphaAtAnchor.value>0.0) {
      map_strcat(line9, 256, " ");
    };
        
    switch (element_iter->convergeReason) {
    case 0 :
      map_strcat(line10, 256, "Improver input parameter for inner loop solver.");
      break;
    case 1 :
      map_strcat(line10, 256, "SQRT( H^2 + V^2 ) <= INNER_FTOL.");
      break;
    case 2 :
      map_strcat(line10, 256, "||x_i - x_{i-1}||  <= INNER_XTOL.");
      break;
    case 3 :
      map_strcat(line10, 256, "Both SQRT( H^2 + V^2 )<=INNER_FTOL and ||x_i - x_{i-1}||  <= INNER_XTOL.");
      break;
    case 4 :
      map_strcat(line10, 256, "The cosine of the angle between [H, V] and any column of the Jacobian is at most INNER_GTOL in absolute value.");
      break;
    case 5 :
      map_strcat(line10, 256, "INNER_MAX_ITS reached. Try increasing the maximum allowable iteration count.");
      break;
    case 6 :
      map_strcat(line10, 256, "INNER_FTOL si too small. SQRT( H^2 + V^2 ) cannot be reduced further.");
      break;
    case 7 :
      map_strcat(line10, 256, "INNER_XTOL is too small. No further refinements can be made in [H, V].");
      break;
    case 8 :
      map_strcat(line10, 256, "INNER_GTOL is too small. The iterated variable vector [H, V] is orthogonal to the Jacobian columns.");
      break;
    default :
      map_strcat(line10, 256, "Uncaught error.");
      break;
    };

    fprintf(file, "                | Element %d\n", i+1);    
    fprintf(file, "                | ---------------------------------------\n");    
    fprintf(file, "Material        |  %s\n", element_iter->lineProperty->label->data);
    fprintf(file, "Lu        [m]   | %s%1.3f\n", line0, element_iter->Lu.value);
    fprintf(file, "Lb        [m]   | %s%1.3f\n", line1, element_iter->lb.value); 
    fprintf(file, "H         [N]   | %s%1.3f\n", line2, *(element_iter->H.value));
    fprintf(file, "V         [N]   | %s%1.3f\n", line3, *(element_iter->V.value));
    fprintf(file, "T         [N]   | %s%1.3f\n", line4, element_iter->T.value);
    fprintf(file, "Alpha     [deg] | %s%1.3f\n", line5, element_iter->alpha.value*RAD2DEG); 
    fprintf(file, "HAnch     [N]   | %s%1.3f\n", line6, element_iter->HAtAnchor.value);
    fprintf(file, "VAnch     [N]   | %s%1.3f\n", line7, element_iter->VAtAnchor.value);
    fprintf(file, "TAnch     [N]   | %s%1.3f\n", line8, element_iter->TAtAnchor.value);
    fprintf(file, "AlphaAnch [deg] | %s%1.3f\n", line9, element_iter->alphaAtAnchor.value*RAD2DEG);
    fprintf(file, "L^2-Norm        |  %1.7g\n", element_iter->residualNorm);
    fprintf(file, "Function Evals  |  %d\n", element_iter->numFuncEvals);
    fprintf(file, "Jacobian Evals  |  %d\n", element_iter->numJacEvals);
    fprintf(file, "Term. criteria  |  %d : %s\n", element_iter->convergeReason, line10);
    fprintf(file, "\n\n");

    line0[0] = 0;
    line1[0] = 0;
    line2[0] = 0;
    line3[0] = 0;
    line4[0] = 0;
    line5[0] = 0;
    line6[0] = 0;
    line7[0] = 0;
    line8[0] = 0;
    line9[0] = 0;
    line10[0] = 0;
  };
  return MAP_SAFE;
};


MAP_ERROR_CODE write_expanded_input_file_to_summary_file(FILE* file, InitializationData* init_data)
{
  const char* line = NULL;
  bstring user_msg = NULL;
  int len = 0;
  int i = 0;

  fprintf(file, "---------------------- LINE DICTIONARY ---------------------------------------\n");
  fprintf(file, "LineType  Diam      MassDenInAir   EA            CB\n");
  fprintf(file, "[-]       [m]       [kg/m]         [N]           [-]\n");  
  for (i=0 ; i<init_data->libraryInputString->qty ; i++) {
    line = init_data->libraryInputString->entry[i]->data;
    len = init_data->libraryInputString->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };

  fprintf(file, "---------------------- NODE PROPERTIES ---------------------------------------\n");
  fprintf(file, "Node  Type       X       Y       Z      M     B     FX    FY    FZ\n");
  fprintf(file, "[-]   [-]       [m]     [m]     [m]    [kg]  [mˆ3]  [N]   [N]   [N]\n");
  for (i=0 ; i<init_data->expandedNodeInputString->qty ; i++) {    
    line = init_data->expandedNodeInputString->entry[i]->data;
    len = init_data->expandedNodeInputString->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };
  
  fprintf(file, "---------------------- LINE PROPERTIES ---------------------------------------\n");
  fprintf(file, "Element  LineType  UnstrLen  NodeAnch  NodeFair  Flags\n");
  fprintf(file, "[-]      [-]       [m]       [-]       [-]       [-]\n");
  for (i=0 ; i<init_data->expandedElementInputString->qty ; i++) {    
    line = init_data->expandedElementInputString->entry[i]->data;
    len = init_data->expandedElementInputString->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };
  
  fprintf(file, "---------------------- SOLVER OPTIONS-----------------------------------------\n");
  fprintf(file, "Option \n");
  fprintf(file, "[-] \n");
  for (i=0 ; i<init_data->solverOptionsString->qty ; i++) {    
    line = init_data->solverOptionsString->entry[i]->data;
    len = init_data->solverOptionsString->entry[i]->slen;
    fprintf(file, "%s", line);
    if (line[len-1]!='\n') {
      fprintf(file, "\n" );
    };
  };
  return MAP_SAFE;
};
