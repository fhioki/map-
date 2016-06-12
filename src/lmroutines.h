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


#ifndef _LMROUNTINES_H
#define _LMROUNTINES_H


#include "mapsys.h"
#include "maperror.h"

// #include "simclist/simclist.h"


MAP_ERROR_CODE reset_element(Element* element);
// MAP_ERROR_CODE initialize_interior_node(Node* node, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_set_node_list(const MAP_ParameterType_t* p_type,  MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_set_line_list(MAP_ConstraintStateType_t* z_type, Domain* domain, struct bstrList* line_input_string, char* map_msg, MAP_ERROR_CODE* ierr);
//MAP_ERROR_CODE lm_allocate_types_for_nodes(MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_set_element_list(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_allocate_input_type(MAP_InputType_t* u_type, int N, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_allocate_continuous_type(MAP_ContinuousStateType_t* x_type, int N, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_allocate_other_type(MAP_OtherStateType_t* other_type, int fix_num, int connect_num, int lm_num, char* map_msg, MAP_ERROR_CODE* ierr);
MAP_ERROR_CODE lm_allocate_output_type(MAP_OutputType_t* y_type, int N, char* map_msg, MAP_ERROR_CODE* ierr);

// #include "Line.hpp" 
// #include "Connection.hpp"

// #include <vector>
// #include <string>
// #include <stdlib.h>


// /* Using a link list for LineList and ConnectList avoids the standard
//  * library vectors. C and C++ stuff can be completely kept separate 
//  * this way. 
//  */
// struct LMContainer_t {
//   list_t LineList;    // vector<Line> LineList;          
//   list_t ConnectList; // vector<Connection> ConnectList;
// }; typedef struct LMContainer_t LMContainer;


#endif /* _LMROUTINES_H */
