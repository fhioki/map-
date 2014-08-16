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


