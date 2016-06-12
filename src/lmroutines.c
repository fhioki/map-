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





MAP_ERROR_CODE reset_element(Element* element)
{
  element->tension = 0.0;
  element->l = 0.0;
  element->r1 = NULL; 
  element->r2 = NULL; 

  return MAP_SAFE;
};


MAP_ERROR_CODE lm_set_line_list(MAP_ConstraintStateType_t* z_type, Domain* domain, struct bstrList* line_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  const int num_lines = line_input_string->qty;
  Line new_line;
  Line* line_iter = NULL;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  bstring alias = NULL;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   
  success = reset_line(&new_line);

  for(i=0 ; i<num_lines ; i++) {         
    list_append(&domain->line, &new_line);
    line_iter = (Line*)list_get_at(&domain->line, i);
    i_parsed = 0;
    next = 0;
    parsed = bsplits(line_input_string->entry[i], &tokens);

    MAP_BEGIN_ERROR_LOG;  

    while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
      if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
        if (next==0) { /* use this first option as an opportunity to set the run-time flags to false */             
          next++;
        } else if (next==1) {
          success = associate_line_with_cable_property(line_iter, domain, (char*)parsed->entry[i_parsed]->data, map_msg, ierr); CHECKERRQ(MAP_FATAL_32);           
          next++;
        } else if (next==2) { 
          alias = bformat("Lu[%d]", i+1);
          success = set_vartype("[m]", alias, i, &line_iter->Lu, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_26);
          success = bdestroy(alias);
          next++;
        } else if (next==3) { 
          line_iter->idx_anch = (int)atoi(parsed->entry[i_parsed]->data); 
          success = associate_line_with_anchor_node(line_iter, domain, i+1, (char*)parsed->entry[i_parsed]->data,  map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
          next++;
        } else if (next==4) { 
          line_iter->idx_fair = (int)atoi(parsed->entry[i_parsed]->data);           
          success = associate_line_with_fairlead_node(line_iter, domain, i+1, (char*)parsed->entry[i_parsed]->data,  map_msg, ierr); CHECKERRQ(MAP_FATAL_32);        
          next++;
        } else { /* set the node mass */            
          success = set_line_option_flags(parsed, &i_parsed, line_iter, map_msg, ierr);
        };
      };
      i_parsed++;
    };

    MAP_END_ERROR_LOG;   

    success = bstrListDestroy(parsed);
  };  
  return MAP_SAFE;
};


MAP_ERROR_CODE lm_set_node_list(const MAP_ParameterType_t* p_type,  MAP_InputType_t* u_type, MAP_ConstraintStateType_t* z_type, MAP_OtherStateType_t* other_type, MAP_OutputType_t* y_type, Domain* domain, struct bstrList* node_input_string, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  const int num_nodes = node_input_string->qty;
  int i = 0;
  int i_parsed = 0;
  int next = 0;
  Node new_node;
  Node* node_iter = NULL;
  struct bstrList* parsed = NULL;
  struct tagbstring tokens; 
  bstring alias = NULL;
  bstring value_string = NULL;
  // const double depth = p_type->depth;
  ReferencePoint u_reference_point;

  cstr2tbstr(tokens," \t\n\r"); /* token for splitting line into indivdual words is a tab and space */   

  //success = allocate_types_for_nodes(u_type, z_type, other_type, y_type, domain, node_input_string, map_msg, ierr);
  success = reset_node(&new_node); /* create an empty node */
   
  for(i=0 ; i<num_nodes ; i++) {         
    list_append(&domain->node, &new_node); /* append node to list */
    node_iter = (Node*)list_get_at(&domain->node, i);
    // success = set_node_vartype(node_iter);
    i_parsed = 0;
    next = 0;
    parsed = bsplits(node_input_string->entry[i], &tokens);

    MAP_BEGIN_ERROR_LOG;  

    while (i_parsed<parsed->qty-1) { /* iterating through all strings */              
      if (parsed->entry[i_parsed]->slen) { /* if the string length is not 0 */
        if (next==0) {            
          next++;
        } else if (next==1) {
          if (biseqcstrcaseless(parsed->entry[i_parsed],"FIX")) {
            node_iter->type = FIX;
          } else if (biseqcstrcaseless(parsed->entry[i_parsed],"CONNECT")) {
            node_iter->type = CONNECT;
          } else if (biseqcstrcaseless(parsed->entry[i_parsed],"VESSEL")) {
            node_iter->type = VESSEL;
          } else {
            set_universal_error_with_message(map_msg, ierr, MAP_FATAL_25, "Value: <%s>", parsed->entry[i_parsed]->data);
          };
          next++;
        } else if (next==2) { /* set initial X node position values */
          node_iter->xi = atof(parsed->entry[i_parsed]->data);
          // alias = bformat("X[%d]", i+1);                          
          // success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.x, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_17);
          // bdestroy(alias);
          next++;
        } else if (next==3) { /* set initial Y node position values */
          node_iter->yi = atof(parsed->entry[i_parsed]->data);
          // alias = bformat("Y[%d]", i+1);                          
          // success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.y, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_18);
          // bdestroy(alias);
          next++;
        } else if (next==4) { /* set initial Z node position values */
          if (biseqcstrcaseless(parsed->entry[i_parsed],"DEPTH")) {         
            if (node_iter->type!=FIX) { /* can only use 'DEPTH' flag in input file for FIX (anchor) nodes */
              set_universal_error_with_message(map_msg, ierr, MAP_FATAL_71, "Value: <%s>", parsed->entry[i_parsed]->data);
            } else {
              node_iter->zi = -p_type->depth;
              // value_string = bformat("%f", -depth);                          
              // success = set_vartype_ptr("[m]", alias, i, &node_iter->position_ptr.z, value_string); CHECKERRQ(MAP_FATAL_19);
              // success = bdestroy(value_string);
            };
          } else { /* all other nodes not using the 'DEPTH' flag */
            node_iter->zi = atof(parsed->entry[i_parsed]->data);
          };                  
          next++;
        } else if (next==5) { /* set the node mass */            
          alias = bformat("M[%d]", i+1);                          
          success = set_vartype("[kg]", alias, i, &node_iter->M_applied, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_20);
          bdestroy(alias);
          next++;  
        } else if (next==6) { /* set the node buoyancy */
          alias = bformat("B[%d]", i+1);                          
          success = set_vartype("[m^3]", alias, i, &node_iter->B_applied, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_21);
          bdestroy(alias);
          next++; 
        } else if (next==7) { /* set applied X external force (or user guess) of the node */                    
          alias = bformat("FX[%d]", i+1);                          
          success = set_vartype("[N]", alias, i, &node_iter->external_force.fx, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_22);
          bdestroy(alias);
          next++;
        } else if (next==8) { /* set applied Y external force (or user guess) of the node */            
          alias = bformat("FY[%d]", i+1);                          
          success = set_vartype("[N]", alias, i, &node_iter->external_force.fy, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_23);
          bdestroy(alias);
          next++;
        } else if (next==9) { /* set applied Z external force (or user guess) of the node */
          alias = bformat("FZ[%d]", i+1);                          
          success = set_vartype("[N]", alias, i, &node_iter->external_force.fz, parsed->entry[i_parsed]); CHECKERRQ(MAP_FATAL_24);
          bdestroy(alias);
          next++;
        } else {            
          next++;
        };
      };
      i_parsed++;
    };
  
    MAP_END_ERROR_LOG;   

    success = bstrListDestroy(parsed);
  };  
  checkpoint();
  MAP_RETURN_STATUS(*ierr);
};



MAP_ERROR_CODE lm_set_element_list(Domain* domain, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  Line* line = NULL;
  int N = 0;
  int i = 0;
  int k = 0;
  // int xnum = 0;
  // int unum = 0;
  // int other_num = 0;
  // int lmnum = 0;
  int size = 0;
  int N_start = 0;
  Element new_element;
  Element* element_iter;
  Node new_node;
  Node* node_iter;

  list_iterator_start(&domain->line); /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->line)) { /* tell whether more values available */ 
    line = (Line*)list_iterator_next(&domain->line);    
    N = line->segment_size;
    N_start = 0; list_size(&domain->node) - 1;

    list_init(&line->element); 
    list_attributes_copy(&line->element, line_element_meter, 1);  

    success = reset_node(&new_node);

    for (i=0 ; i<N-1 ; i++) {
      list_append(&domain->node, &new_node); /* append node to list */
      size = list_size(&domain->node);
      node_iter = (Node*)list_get_at(&domain->node, size-1);
      node_iter->type = LM;
    };
    // lmnum += (N-1);
    domain->lm_num += (N-1);

    /* allocate elements */
    for (i=0 ; i<N ; i++) {
      k = N_start + i;
      list_append(&line->element, &new_element); /* append node to list */
      element_iter = (Element*)list_get_at(&line->element, i);
      success = reset_element(element_iter);
    
      if (i==0) {
        node_iter = (Node*)list_get_at(&domain->node, k);
        element_iter->r1 = line->anchor;
        element_iter->r2 = node_iter;
      } else if (i==N-1) {
        node_iter = (Node*)list_get_at(&domain->node, k);
        element_iter->r1 = node_iter;
        element_iter->r2 = line->fairlead;
      } else {
        element_iter->r1 = (Node*)list_get_at(&domain->node, k);
        element_iter->r2 = (Node*)list_get_at(&domain->node, k+1);      
      };
    };
  };
  list_iterator_stop(&domain->line); /* ending the iteration session */  

  list_iterator_start(&domain->node); /* starting an iteration "session" */
  while (list_iterator_hasnext(&domain->node)) { /* tell whether more values available */ 
    node_iter = (Node*)list_iterator_next(&domain->node);        
    if ((node_iter->type!=FIX) && (node_iter->type!=VESSEL)) {
      if (node_iter->type==LM) {
        (domain->n_x)++;
      } else if (node_iter->type==CONNECT) {      
        (domain->connect_num)++;
        if (node_iter->assigned==false) {
          (domain->n_x)++;
          node_iter->assigned = true;
        }; 
      } else {
        /* throw error */
      };
    } else {
      if (node_iter->type==VESSEL) {
        (domain->n_y)++;
        (domain->n_u)++;
      } else if (node_iter->type==FIX) {
        (domain->fix_num)++;
      } else {
        /* throw error */
      }
    };
  };
  list_iterator_stop(&domain->node); /* ending the iteration session */  

  // printf("number of lm force: %d\n", domain->lm_num);
  // printf("number of fix force: %d\n", domain->fix_num);
  // printf("number of connect force: %d\n", domain->connect_num);
  // printf("number of output: %d\n", domain->n_y);
  // printf("number of output: %d\n", domain->n_u);
  // printf("number of derivs: %d\n", domain->n_x);

  return MAP_SAFE;
};


MAP_ERROR_CODE lm_allocate_input_type(MAP_InputType_t* u_type, int N, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  u_type->x_Len = N;
  u_type->y_Len = N;
  u_type->z_Len = N;
  u_type->xd_Len = N;
  u_type->yd_Len = N;
  u_type->zd_Len = N;
  u_type->x = malloc(u_type->x_Len*sizeof(double));
  u_type->y = malloc(u_type->y_Len*sizeof(double));
  u_type->z = malloc(u_type->z_Len*sizeof(double));
  u_type->xd = malloc(u_type->xd_Len*sizeof(double));
  u_type->yd = malloc(u_type->yd_Len*sizeof(double));
  u_type->zd = malloc(u_type->zd_Len*sizeof(double));
  return MAP_SAFE;
};


MAP_ERROR_CODE lm_allocate_continuous_type(MAP_ContinuousStateType_t* x_type, int N, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  x_type->rd_Len = N*3; // multiplied by three because three equations per node
  x_type->rdd_Len = N*3;
  x_type->rd = malloc(x_type->rd_Len*sizeof(double));      
  x_type->rdd = malloc(x_type->rdd_Len*sizeof(double));     

  return MAP_SAFE;
};


MAP_ERROR_CODE lm_allocate_other_type(MAP_OtherStateType_t* other_type, int fix_num, int connect_num, int lm_num, char* map_msg, MAP_ERROR_CODE* ierr)
{
  printf("number of fixed states: %d\n",fix_num+connect_num+lm_num);
  other_type->x_Len = fix_num;
  other_type->y_Len = fix_num;
  other_type->z_Len = fix_num;
  other_type->x = malloc(other_type->x_Len*sizeof(double));
  other_type->y = malloc(other_type->y_Len*sizeof(double));
  other_type->z = malloc(other_type->z_Len*sizeof(double));

  /* If the node is VESSEL, then the applied force is an output state. Otherwise, 
   * it has to be an other state because it can't be associated with any other type. 
   * This is what is done below.
   */
  other_type->Fx_connect_Len = connect_num;
  other_type->Fy_connect_Len = connect_num;
  other_type->Fz_connect_Len = connect_num;
  other_type->Fx_connect = malloc(other_type->Fx_connect_Len*sizeof(double));
  other_type->Fy_connect = malloc(other_type->Fy_connect_Len*sizeof(double));
  other_type->Fz_connect = malloc(other_type->Fz_connect_Len*sizeof(double));
   
  other_type->Fx_anchor_Len = fix_num;
  other_type->Fy_anchor_Len = fix_num;
  other_type->Fz_anchor_Len = fix_num;
  other_type->Fx_anchor = malloc(other_type->Fx_anchor_Len*sizeof(double));
  other_type->Fy_anchor = malloc(other_type->Fy_anchor_Len*sizeof(double));
  other_type->Fz_anchor = malloc(other_type->Fz_anchor_Len*sizeof(double));

  other_type->Fx_lm_Len = lm_num;
  other_type->Fy_lm_Len = lm_num;
  other_type->Fz_lm_Len = lm_num;
  other_type->Fx_lm = malloc(other_type->Fx_lm_Len*sizeof(double));
  other_type->Fy_lm = malloc(other_type->Fy_lm_Len*sizeof(double));
  other_type->Fz_lm = malloc(other_type->Fz_lm_Len*sizeof(double));

  return MAP_SAFE;
};


MAP_ERROR_CODE lm_allocate_output_type(MAP_OutputType_t* y_type, int N, char* map_msg, MAP_ERROR_CODE* ierr) 
{
  y_type->Fx_Len = N;
  y_type->Fy_Len = N;
  y_type->Fz_Len = N;
  y_type->Fx = malloc(y_type->Fx_Len*sizeof(double));
  y_type->Fy = malloc(y_type->Fy_Len*sizeof(double));
  y_type->Fz = malloc(y_type->Fz_Len*sizeof(double));

  return MAP_SAFE;
};

