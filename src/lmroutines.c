/****************************************************************
 *   Copyright (C) 2014                                         *
 *                                                              *
 * Licensed to the Apache Software Foundation (ASF) under one   *
 * or more contributor license agreements.  See the NOTICE file *
 * distributed with this work for additional information        *
 * regarding copyright ownership.  The ASF licenses this file   *
 * to you under the Apache License, Version 2.0 (the            *
 * "License"); you may not use this file except in compliance   *
 * with the License.  You may obtain a copy of the License at   *
 *                                                              *
 *   http://www.apache.org/licenses/LICENSE-2.0                 *
 *                                                              *
 * Unless required by applicable law or agreed to in writing,   *
 * software distributed under the License is distributed on an  *
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       *
 * KIND, either express or implied.  See the License for the    *
 * specific language governing permissions and limitations      *      
 * under the License.                                           *  
 ****************************************************************/


#include "lmroutines.h"
#include "mapinit.h"


MAP_ERROR_CODE allocate_lumped_mass_data(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  Line* line_iter = NULL;

  checkpoint();
  list_iterator_start(&domain->line); /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);    

    list_init(&line_iter->interior_node); 
    list_init(&line_iter->element); 

    list_attributes_copy(&line_iter->interior_node, node_meter, 1);  
    list_attributes_copy(&line_iter->element, line_element_meter, 1);  
  };
  list_iterator_stop(&domain->line); /* ending the iteration session */  

  return MAP_SAFE;
}
