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


#include "map.h"
#include "mapsys.h"
#include "maperror.h"
#include "protos.h"


const char* remove_first_character(const char* str)
{
  return str+1;
};



// /**
//  *
//  */
// void initialize_vartype(char* unit, char* alias, VarType* type, const int num)
// {
//   int cx = 0;
//   char buffer[64] = "";
// 
//   cx = map_snprintf(buffer, 64, "%s[%d]", alias, num); assert(cx>=0);
//   type->name = malloc(sizeof(char)*strlen(buffer)+1);
//   type->units = malloc(sizeof(char)*strlen(unit)+1);
//   strcpy(type->name, buffer);
//   strcpy(type->units, unit);
//   type->referenceCounter = 0;
//   type->id = num;  
//   type->value = -999.9;
//   type->isFixed = false;
// };


// /**
//  *
//  */
// MAP_ERROR_CODE set_vartype(char* unit, char* alias, const int num, VarType* type, char const* property)
// {
//   const char* return_str = NULL;
// 
//   type->name = malloc(sizeof(char)*strlen(alias)+1);
//   type->units = malloc(sizeof(char)*strlen(unit)+1);
//   
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);
// 
//   type->referenceCounter = 0;
//   type->id = num;  
//   type->value = 0;
// 
//   if (property[0]== '#') { 
//     /* this variable is an iterated parameter */  
//     type->isFixed = false;
//     return_str = remove_first_character(property);
//     if (!strcmp(return_str, "") || !strcmp(return_str, "\n")) {
//       /* entry in empty in the MAP input file */
//       type->value = -999.9;
//       return MAP_SAFE;
//     } else if (is_numeric(return_str)) { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else {
//       /* converted to numeric value */
//       type->value = (MapReal)atof(return_str); 
//       return MAP_SAFE;
//     };
//   } else { 
//     /* this variable is constant */
//     type->isFixed = true;
//     if (is_numeric(property)) 
//     { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else { 
//       /* converted to numeric value */
//       type->value = (MapReal)atof(property);
//       return MAP_SAFE;
//     };
//   };
//   return MAP_SAFE;
// };


/**
 *
 */
MAP_ERROR_CODE associate_vartype_ptr(VarTypePtr* type, double* arr, int index)
{
  type->value = &arr[index-1];
  return MAP_SAFE;
};



// /**
//  *
//  */
// MAP_ERROR_CODE set_vartype_ptr(char* unit, char* alias, const int num, VarTypePtr* type, char const* property)
// {
//   const char* return_str = NULL;
//   
//   type->name = malloc(sizeof(char)*strlen(alias)+1); /* @todo: this needs to be freed */
//   type->units = malloc(sizeof(char)*strlen(unit)+1); /* @todo: this needs to be freed */
//   
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);
//   
//   type->referenceCounter = 0;
//   type->id = num;  
// 
//   if (property[0]=='#') { 
//     /* this variable is an iterated parameter */      
//     type->isFixed = false;
//     return_str = remove_first_character(property);
//     if (!strcmp(return_str, "") || !strcmp(return_str, "\n")) {
//       /* entry in empty in the MAP input file */
//       *type->value = -999.9;
//       return MAP_SAFE;
//     } else if (is_numeric(return_str)) { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else {
//       /* converted to numeric value */
//       *type->value = (MapReal)atof(return_str); 
//       return MAP_SAFE;
//     };
//   } else { 
//     /* this variable is constant */
//     type->isFixed = true;
//     if (is_numeric(property)) { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else { 
//       /* converted to numeric value */
//       *type->value = (MapReal)atof(property); 
//       return MAP_SAFE;
//     };
//   };
//   return MAP_SAFE;
// };








// MAP_ERROR_CODE set_vartype_ptr_b(const bstring* unit, const bstring* alias, const int num, VarTypePtr* type, const bstring* property)
// {
//   const char* return_str = NULL;
//   
//   type->name = malloc(sizeof(char)*strlen(alias)+1); /* @todo: this needs to be freed */
//   type->units = malloc(sizeof(char)*strlen(unit)+1); /* @todo: this needs to be freed */
//   
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);
//   
//   type->referenceCounter = 0;
//   type->id = num;  
// 
//   if (property->data[0]=='#') { 
//     /* this variable is an iterated parameter */      
//     type->isFixed = false;
//     return_str = remove_first_character(property);
//     if (!strcmp(return_str, "") || !strcmp(return_str, "\n")) {
//       /* entry in empty in the MAP input file */
//       *type->value = -999.9;
//       return MAP_SAFE;
//     } else if (is_numeric(return_str)) { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else {
//       /* converted to numeric value */
//       *type->value = (MapReal)atof(return_str); 
//       return MAP_SAFE;
//     };
//   } else { 
//     /* this variable is constant */
//     type->isFixed = true;
//     if (is_numeric(property)) { 
//       /* cannot convert to numeric value */
//       return MAP_FATAL;
//     } else { 
//       /* converted to numeric value */
//       *type->value = (MapReal)atof(property); 
//       return MAP_SAFE;
//     };
//   };
//   return MAP_SAFE;
// };





// /**
//  *
//  */
// MAP_ERROR_CODE set_vartype_ptr_float(char* unit, char* alias, const int num, VarTypePtr* type, const MapReal value)
// {
//   const char* return_str = NULL;
//   
//   type->name = malloc(sizeof(char)*strlen(alias)+1); /* @todo: this needs to be freed */
//   type->units = malloc(sizeof(char)*strlen(unit)+1); /* @todo: this needs to be freed */
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);  
// 
//   type->referenceCounter = 0;
//   type->id = num;  
//   type->isFixed = true;
//   *type->value = -value;
//   return MAP_SAFE;
// };


// MAP_ERROR_CODE set_element_vartype_ptr(char* unit, char* alias, const int num, VarTypePtr* type)
// {  
//   type->name = malloc(sizeof(char)*strlen(alias)+1); 
//   type->units = malloc(sizeof(char)*strlen(unit)+1); 
//   if (type->name==NULL || type->units==NULL) {
//     return MAP_FATAL;
//   };
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);  
//   type->referenceCounter = 0;
//   type->id = num;  
//   *(type->value) = -999.9;
//   return MAP_SAFE;
// };
// 
// 
// MAP_ERROR_CODE set_vartype_float(char* unit, char* alias, const double num, VarType* type, MapReal const value)
// {
//   type->name = malloc(sizeof(char)*strlen(alias)+1);
//   type->units = malloc(sizeof(char)*strlen(unit)+1);
// 
//   strcpy(type->name, alias);
//   strcpy(type->units, unit);
// 
//   type->isFixed = false;
//   type->referenceCounter = 0;
//   type->id = num;  
//   type->value = value;
// 
//   return MAP_SAFE;
// };


/**
 *
 */
int strcicmp(char const* a, char const* b)
{
  for (; ; a++,b++) {
    int d = tolower(*a) - tolower(*b);
    if (d!=0 || !*a) {
      return d;
    };
  };
};


/**
 * Check if a string can be converted to a valid numeric value
 *
 * 0 (MAP_SAFE) is reserved for 'safe' exit out of functions. 
 *
 * @param s     must point to a user-provided memory location
 * @return      0 for success, 3 for failure
 *
 * @see MAP_ERROR_CODE
 */
MAP_ERROR_CODE is_numeric(const char* s)
{
  char* p = NULL;
  if (s==NULL || *s=='\0' || isspace(*s)) {
    return MAP_SAFE;
  };
  strtod (s, &p);
  if (*p=='\0') {
    return MAP_FATAL;
  } else {
    return MAP_SAFE;
  };
};


/**
 * done
 */
void copy_target_string(char* target, unsigned char* source)
{
  while (*source) {
    *target = *source;
    source++;
    target++;
  };
  *target = '\0';
};


/**
 * done
 */
MAP_ERROR_CODE map_get_version(MAP_InitOutputType_t* ioType)
{
  bstring out_string = NULL;
  int ret = 0;

  /* first set the program version defined in the mapsys.h header file 
   * @todo: program version should be tied to the gi revision numner
   */
  out_string = bformat("<%s>",PROGVERSION);
  if (out_string->slen>MAX_INIT_VERSION_STRING_LENGTH) { /* overflow */
    return MAP_FATAL; /* @todo: give proper error code */
  };
  copy_target_string(ioType->version, out_string->data);
  ret = bdestroy(out_string);

  /* the set the compiling date. This is #defined in the mapsys.h header */
  out_string = bformat("<%c%c%c-%c%c-%c%c%c%c>",BUILD_MONTH_CH0,BUILD_MONTH_CH1,BUILD_MONTH_CH2,BUILD_DAY_CH0,BUILD_DAY_CH1,BUILD_YEAR_CH0,BUILD_YEAR_CH1,BUILD_YEAR_CH2,BUILD_YEAR_CH3);
  if (out_string->slen>MAX_INIT_COMPILING_DATA_STRING_LENGTH) { /* overflow */
    return MAP_FATAL; /* @todo: give proper error code */
  };
  copy_target_string(ioType->compilingData, out_string->data);
  ret = bdestroy(out_string);
  return MAP_SAFE;
};


void print_machine_name_to_screen( ) {
  char name[50]="";

  // __get_machine_name(name);
  printf( "%s Ver. %s ", PROGNAME, PROGVERSION); 
  printf( "%c",BUILD_MONTH_CH0 );// build month
  printf( "%c",BUILD_MONTH_CH1 );
  printf( "%c",BUILD_MONTH_CH2 );
  printf( "-" );
  printf( "%c",BUILD_DAY_CH0 );// build day
  printf( "%c",BUILD_DAY_CH1 );
  printf( "-" );
  printf( "%c",BUILD_YEAR_CH0 ); // build year 
  printf( "%c",BUILD_YEAR_CH1 );
  printf( "%c",BUILD_YEAR_CH2 );
  printf( "%c\n",BUILD_YEAR_CH3 );
}


/**
 * @see: http://stackoverflow.com/questions/504810/how-do-i-find-the-current-machines-full-hostname-in-c-hostname-and-domain-info
 */
void __get_machine_name(char* machineName)
{
  char name[150];  
  
#if defined(_WIN32) || defined(_WIN64)
  int i = 0;
  TCHAR infoBuf[150];
  DWORD bufCharCount = 150;
  memset(name, 0, 150);
  if (GetComputerName(infoBuf, &bufCharCount)) {
    for (i=0 ; i<150 ; i++) {
      name[i] = infoBuf[i];
    };
  } else {
    strcpy(name, "Unknown_Host_Name");
  };
#else
  memset(name, 0, 150);
  gethostname(name, 150);
#endif
  strncpy(machineName ,name, 150);
};


/**
 *
 */
MAP_ERROR_CODE print_help_to_screen()
{
  print_machine_name_to_screen( );

  printf("MAP Input file section defitions:\n");
  printf("  Line dictionary definitions:\n");   
  printf("    -LineType, --User-defined name of line [-]  \n");   
  printf("    -Diam,     --Line diamater, used to calulate area and line displacement per unit length [m]  \n");   
  printf("    -MassDen,  --Mass (in air) per unit length [kg/m]  \n");   
  printf("    -EA,       --Axial stiffness [N] \n");   
  printf("    -CB,       --Cable/seabed Coulumb friction coefficient [-]  \n");   
  printf("    -CIntDamp, --Internal structural damping coefficient [Pa-s]  \n");   
  printf("    -Ca,       --Cross-flow added-mass coeficient [-]\n");   
  printf("    -Cdn,      --Cross-flow drag coefficient [-]\n");   
  printf("    -Cdt,      --Tangent (skin) drag coefficient[-]\n");   
  printf("  Node property definitions:\n");
  printf("    -Node,     --Node number; first starts at 1 [-]\n");
  printf("    -Type,     --Type of node. Must be one of: VESSEL, FIX, CONNECT [-]\n");
  printf("    -X,        --Node X position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -Y,        --Node Y position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -Z,        --Node Z position. '#' must prefix CONNECT nodes; constitutes user initial guess [m]\n");
  printf("    -M,        --Applied point mass at node [kg]\n");
  printf("    -B,        --Applied point buoyancy module at node [m^3]\n");  
  printf("    -FX,       --Applied X external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("    -FY,       --Applied Y external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("    -FZ,       --Applied Z external force at node. '#' must prefix VESSEL and FIX nodes [N]\n");
  printf("  Element property definitions:\n");
  printf("    -Element,  --Element number; first starts at 1 [-]\n");
  printf("    -LineType, --Must match property defined in 'Line Dictions'[-]\n");
  printf("    -UnstrLen, --Unstretched line length [m]\n");
  printf("    -NodeAnch, --Anchor node number corresponding to 'Node Property Definitions' section [-]\n");
  printf("    -NodeFair, --Fairlead node number corresponding to 'Node Property Definitions' section [-]\n");
  printf("    -Flags,    --User run-time flag; see below [-]\n");

  printf("    \n");
  printf("  Element run-time options definitions\n");
  printf("    Outputs:\n");
  printf("      -gx_pos,       --Fairlead posiion in global X [m]\n");
  printf("      -gy_pos,       --Fairlead posiion in global Y [m]\n");
  printf("      -gx_pos,       --Fairlead posiion in global Z [m]\n");
  printf("      -gx_a_pos,     --Anchor posiion in global X [m]\n");
  printf("      -gy_a_pos,     --Anchor posiion in global Y [m]\n");
  printf("      -gz_a_pos,     --Anchor posiion in global Z [m]\n");
  printf("      -gx_force,     --Fairlead force in global X (include applied forces) [N]\n");
  printf("      -gy_force,     --Fairlead force in global Y (include applied forces) [N]\n");
  printf("      -gz_force,     --Fairlead force in global Z (include applied forces) [N]\n");
  printf("      -h_fair,       --Horizontal force at fairlead (does NOT include applied forces) [N]\n");
  printf("      -v_fair,       --Vertical force at fairlead (does NOT include applied forces) [N]\n");
  printf("      -h_anch,       --Horizontal force at anchor (does NOT include applied forces) [N]\n");
  printf("      -v_anch,       --Vertical force at anchor (does NOT include applied forces) [N]\n");
  printf("      -tension_fair, --Line force-magnitude at fairlead (include applied loads) [N]\n");
  printf("      -tension_anch, --Line force-magnitude at anchor (include applied loads) [N]\n");
  printf("      -azimuth,      --Line lateral offset angle global X axis [deg]\n");
  printf("      -altitude,     --Line inclination angle relative to global XY plane at fiarlead [deg]\n");
  printf("      -lay_length,   --Length of line on seabed [m]\n");
  printf("      -line_tension, -- \n");
  printf("    Model features:\n");
  printf("      -omit_contact,       --Ignore cable/seabed contact\n");
  printf("      -seg_size <10>,      --Number of discrete elements in line\n");
  printf("      -damage_time <NULL>, --Line breakage occurs at specified time [s]\n");
  printf("      -diagnostic,         --Run line solver diagnostics until specified time [s]\n");
  printf("\n");
  printf("  Model option definitions\n");
  printf("    General model features:\n");
  printf("      -ref_position <0.0> <0.0> <0.0>\n");
  printf("      -repeat <NULL> ... <NULL>\n");
  printf("    MSQS solver options:\n");
  printf("      -inner_ftol <float>,\n");
  printf("      -inner_gtol <float>,\n");
  printf("      -inner_xtol <float>,\n");
  printf("      -inner_max_its <int>,\n");
  printf("      -outer_tol <float>,\n");
  printf("      -outer_max_its <int>,\n");
  printf("      -outer_epsilon <float>,\n");
  printf("      -outer_bd,\n");
  printf("      -outer_cd,\n");
  printf("      -outer_fd,\n");
  printf("      -pg_cooked <1000.0> <1.0>,\n");
  printf("      -integration_dt <0.01>,\n");
  printf("    LM model feature (not suported yet):\n");
  printf("      -kb_default      --Seabed stiffness parameter\n");
  printf("      -cb_default      --Seabed damping parameter\n");
  printf("      -wave_kinematics --Enables wave kinematics to drag interaction from surface waves\n");
  printf( "\nMAP++ Copyright (C) 2014 and GNU GPL by Marco Masciola and others\n" );
  printf( "SimCList Copyright (C) 2010 by Mij <http://mij.oltrelinux.com/devel/simclist/>\n" );
  printf( "MinPack Copyright (C) 1999 by the University of Chicago\n" );
  printf( "Modifications to MinPack by Frederic Devernay <http://devernay.free.fr/hacks/cminpack/>\n" );
  printf( "\nMAP++ is free software; see the source for copying conditions.\n" );
  printf( "This software is distributed on an \"AS IS\" BASIS, WITHOUT WARRANTIES\n" );
  printf( "OR CONDITIONS OF ANY KIND, either express or implied. See\n" );
  printf( "<https://gnu.org/licenses/gpl.html> for more details.\n" );
  printf("    \n");
  return MAP_SAFE;
};



