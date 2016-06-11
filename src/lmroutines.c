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
  MAP_ERROR_CODE success = MAP_SAFE;
  Line* line_iter = NULL;
  int N = 0;
  int i = 0;
  Element new_element;
  Element* element_iter;
  Node new_int_node;
  Node* int_node_iter;

  list_iterator_start(&domain->line); /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line_iter = (Line*)list_iterator_next(&domain->line);    
    N = line_iter->segment_size;

    /* allocate nodes */
    list_init(&line_iter->interior_node); 
    list_attributes_copy(&line_iter->interior_node, node_meter, 1);  
    for (i=0 ; i<N-1 ; i++) {
      list_append(&line_iter->interior_node, &new_int_node); /* append node to list */
      int_node_iter = (Node*)list_get_at(&line_iter->interior_node, i);
      success = reset_node(int_node_iter);
    };

    /* allocate elements */
    list_init(&line_iter->element); 
    list_attributes_copy(&line_iter->element, line_element_meter, 1);  
    for (i=0 ; i<N ; i++) {
      list_append(&line_iter->element, &new_element); /* append node to list */
      element_iter = (Element*)list_get_at(&line_iter->element, i);
      success = reset_element(element_iter);
    };

  };
  list_iterator_stop(&domain->line); /* ending the iteration session */  

  return MAP_SAFE;
};


MAP_ERROR_CODE reset_element(Element* element)
{
  element->tension = 0.0;
  element->l = 0.0;
  element->r1 = NULL; 
  element->r2 = NULL; 

  return MAP_SAFE;
};

