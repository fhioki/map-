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
#include "maperror.h"
#include "protos.h"
#include "pyprotos.h"


extern const char MAP_ERROR_STRING[][1024];


/**
 *
 */
MAP_EXTERNCALL InitializationData* MAP_InitInput_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  InitializationData* new_data = NULL;
  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(InitializationData));
  if (new_data == NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_4);    
    return new_data;
  } else {
    initialize_init_data_to_null(new_data);
    *ierr = MAP_SAFE;
    return new_data;
  };
};


/**
 *
 */
MAP_EXTERNCALL void MAP_InitInput_Delete(InitializationData* init_data)
{
  MAPFREE(init_data); 
};


MAP_EXTERNCALL ModelData* MAP_OtherState_Create(char* map_msg, MAP_ERROR_CODE* ierr)
{
  ModelData* new_data = NULL;
  map_reset_universal_error(map_msg, ierr);
  new_data = malloc(sizeof(ModelData));
  if (new_data==NULL) {
    *ierr = map_set_universal_error("", map_msg, ierr, MAP_FATAL_43);    
    return new_data;
  } else {
    new_data->sizeOfCableLibrary = -9999;
    new_data->sizeOfElements = -9999;
    new_data->sizeOfNodes = -9999;    
    return new_data;
  };
};


/**
 * @todo: slated for removal; deprecated when above funtion is removed
 */
MAP_EXTERNCALL void MAP_OtherState_Delete(ModelData* model_data)
{
  MAPFREE(model_data);
};
