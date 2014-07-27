/*
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


#include <time.h>
#include "map.h"
#include "maperror.h"
#include "protos.h"


/**
 *
 */
MAP_ERROR_CODE write_summary_file(InitializationData* init, MAP_ParameterType_t* paramType, ModelData* data, char* map_msg, MAP_ERROR_CODE* ierr)
{
  MAP_ERROR_CODE success = MAP_SAFE;
  struct tm* tmInfo;
  time_t timer;
  FILE* file = fopen(init->summaryFileName, "w");
  char timeBuffer[64] = "";
  char name[256] = "";
  
  time(&timer);
  tmInfo = localtime(&timer);
  
  if (file==NULL) {
    *ierr = map_set_universal_error(init->summaryFileName, map_msg, ierr, MAP_FATAL_38);
    return MAP_FATAL;
  };
  
  __get_machine_name(name);  
  strftime(timeBuffer, 64, "%A %B %d-%Y at %H:%M:%S %p", tmInfo);
  fprintf(file,"Outputs were generated using %s %s-%s on %s.\n\n\n", PROGNAME, PROGVERSION, name, timeBuffer ); 
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


/**
 *
 */
MAP_ERROR_CODE write_cable_library_information_to_summary_file(FILE* file, ModelData* dataObj)
{
  CableLibrary* cableLibraryIter = NULL;  
  list_iterator_start(&dataObj->cableLibrary);    
  while (list_iterator_hasnext(&dataObj->cableLibrary)) { 
    cableLibraryIter = (CableLibrary*)list_iterator_next(&dataObj->cableLibrary);    
    fprintf(file, "    Cable Type          : %s\n", cableLibraryIter->label);
    fprintf(file, "    Diameter     [m]    : %1.4f\n", cableLibraryIter->diam);
    fprintf(file, "    Mass Density [kg/m] : %1.2f\n", cableLibraryIter->massDensityInAir);
    fprintf(file, "    EA           [kN]   : %1.2f\n", cableLibraryIter->ea);
    fprintf(file, "    omega        [N/m]  : %1.2f\n", cableLibraryIter->omega);
    fprintf(file, "    CB                  : %1.2f\n\n", cableLibraryIter->cb);
  };
  list_iterator_stop(&dataObj->cableLibrary);  
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_type_to_summary_file(const int columnNumber, const int countToFour, const NodeType nodeType, char* line)
{
  int j = 0;
  int size = 0;

  if (!columnNumber) { /* if this is the firs column, then line the right side with information */
    if (nodeType==VESSEL) {
      map_snprintf(line, 256, "Type      |  VESSEL");
    } else if (nodeType==FIX ) { 
      map_snprintf(line, 256, "Type      |  FIX");
    } else if (nodeType==CONNECT ) {
      map_snprintf(line, 256, "Type      |  CONNECT");
    } else {
      map_snprintf(line, 256, "Type      |  UNKNOWN");
    };
  } else { /* otherwise, just print the information */
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) { /* add white spaces */
      map_strcat(line, 256, " "); 
    };

    if (nodeType==VESSEL) {
      map_strcat(line, 256, " VESSEL");        
    } else if (nodeType==FIX ) { 
      map_strcat(line, 256, " FIX");        
    } else if(nodeType==CONNECT ) {
      map_strcat(line, 256, " CONNECT");        
    } else {
      map_strcat(line, 256, " UNKNOWN");        
    };
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_header_to_summary_file(const int columnNumber, const int countToFour, const int nodeNumber, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";

  if (!columnNumber) { 
    map_snprintf(line, 256, "          | Node %d Data", nodeNumber);        
  } else {
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " ");       
    };
    map_snprintf(buffer, 64, "%sNode %d Data", spaces, nodeNumber);
    map_strcat(line, 256, buffer);       
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_x_position_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* xPosition, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal x = *xPosition->value;

  if (!columnNumber) { 
    if (!xPosition->isFixed) {
      if (x>=0.0) {
        map_snprintf(line, 256, "X  [m]    | ( %1.3f)", x);
      } else {
        map_snprintf(line, 256, "X  [m]    | (%1.3f)", x);        
      };
    } else {
      if (x>0.0) {
        map_snprintf(line, 256, "X  [m]    |   %1.3f", x);        
      } else {
        map_snprintf(line, 256, "X  [m]    |  %1.3f", x);        
      };
    };
  } else { 
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!xPosition->isFixed) {
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
    map_strcat(line, 256, buffer);       
  };

  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_y_position_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* yPosition, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal y = *yPosition->value;

  if (!columnNumber) { 
    if (!yPosition->isFixed) {
      if (y>=0.0) {
        map_snprintf(line, 256, "Y  [m]    | ( %1.3f)", y);
      } else {
        map_snprintf(line, 256, "Y  [m]    | (%1.3f)", y);        
      };
    } else {
      if (y>=0.0) {
        map_snprintf(line, 256, "Y  [m]    |   %1.3f", y);        
      } else {
        map_snprintf(line, 256, "Y  [m]    |  %1.3f", y);        
      };
    };
  } else { 
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!yPosition->isFixed) {
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
    map_strcat(line, 256, buffer);       
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_z_position_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* zPosition, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal z = *zPosition->value;
  
  if (!columnNumber) { 
    if (!zPosition->isFixed) {
      if (z>=0.0) {
        map_snprintf(line, 256, "Z  [m]    | ( %1.3f)", z);
      } else {
        map_snprintf(line, 256, "Z  [m]    | (%1.3f)", z);        
      };
    } else {
      if (z>=0.0) {
        map_snprintf(line, 256, "Z  [m]    |   %1.3f", z);        
      } else {
        map_snprintf(line, 256, "Z  [m]    |  %1.3f", z);        
      };
    };
  } else { 
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!zPosition->isFixed) {
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
    map_strcat(line, 256, buffer);       
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_mass_information_to_summary_file(const int columnNumber, const int countToFour, VarType* nodePointMass, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal mass = nodePointMass->value;

  if (!columnNumber) { 
    if (!nodePointMass->isFixed) {
      if (mass>=0.0) {
        map_snprintf(line, 256, "M  [kg]   | ( %1.3f)", mass);
      } else {
        map_snprintf(line, 256, "M  [kg]   | (%1.3f)", mass);        
      };
    } else {
      if (mass>=0.0) {
        map_snprintf(line, 256, "M  [kg]   |   %1.3f", mass);        
      } else {
        map_snprintf(line, 256, "M  [kg]   |  %1.3f", mass);        
      };
    };

  } else { 
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };

    if (!nodePointMass->isFixed) {
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
    map_strcat(line, 256, buffer);       
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_buoyancy_information_to_summary_file(const int columnNumber, const int countToFour, VarType* nodeBuoyancyModule, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal buoyancy = nodeBuoyancyModule->value;

  if (!columnNumber) { 
    if (!nodeBuoyancyModule->isFixed) {
      if (buoyancy>=0.0) {
        map_snprintf(line, 256, "B  [m^3]  | ( %1.3f)", buoyancy);
      } else {
        map_snprintf(line, 256, "B  [m^3]  | (%1.3f)", buoyancy);        
      };
    } else {
      if (buoyancy>=0.0) {
        map_snprintf(line, 256, "B  [m^3]  |   %1.3f", buoyancy);        
      } else {
        map_snprintf(line, 256, "B  [m^3]  |  %1.3f", buoyancy);        
      };
    };
  } else { 
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    };
    
    if (!nodeBuoyancyModule->isFixed) {
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
    map_strcat(line, 256, buffer);       
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_x_sum_force_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* xSumForce, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal fx = *(xSumForce->value);

  if (!columnNumber) { 
    if (!xSumForce->isFixed) {
      if (fx>=0.0) {
        map_snprintf(line, 256, "FX [kN]   | ( %1.3f)", fx);
      } else {
        map_snprintf(line, 256, "FX [kN]   | (%1.3f)", fx);        
      };
    } else {
      if (fx>=0.0) {
        map_snprintf(line, 256, "FX [kN]   |   %1.3f", fx);        
      } else {
        map_snprintf(line, 256, "FX [kN]   |  %1.3f", fx);        
      };
    };
  } else {
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    }

    if (!xSumForce->isFixed) {
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
    map_strcat(line, 256, buffer);
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_y_sum_force_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* ySumForce, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal fy = *(ySumForce->value);

  if (!columnNumber) { 
    if (!ySumForce->isFixed) {
      if (fy>=0.0) {
        map_snprintf(line, 256, "FY [kN]   | ( %1.3f)", fy);
      } else {
        map_snprintf(line, 256, "FY [kN]   | (%1.3f)", fy);        
      };
    } else {
      if (fy>=0.0) {
        map_snprintf(line, 256, "FY [kN]   |   %1.3f", fy);        
      } else {
        map_snprintf(line, 256, "FY [kN]   |  %1.3f", fy);        
      };
    };
  } else {
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    }

    if (!ySumForce->isFixed) {
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
    map_strcat(line, 256, buffer);
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_z_sum_force_to_summary_file(const int columnNumber, const int countToFour, VarTypePtr* zSumForce, char* line)
{
  int j = 0;
  int size = 0;
  char spaces[64] = "";
  char buffer[64] = "";
  MapReal fz = *(zSumForce->value);

  if (!columnNumber) { 
    if (!zSumForce->isFixed) {
      if (fz>=0.0) {
        map_snprintf(line, 256, "FZ [kN]   | ( %1.3f)", fz);
      } else {
        map_snprintf(line, 256, "FZ [kN]   | (%1.3f)", fz);        
      };
    } else {
      if (fz>=0.0) {
        map_snprintf(line, 256, "FZ [kN]   |   %1.3f", fz);        
      } else {
        map_snprintf(line, 256, "FZ [kN]   |  %1.3f", fz);        
      };
    };
  } else {
    size = (int)strlen(line)-23*countToFour;
    for (j=0 ; j<SPACE_LENGTH-size ; j++) {
      map_strcat(spaces, 64, " "); 
    }

    if (!zSumForce->isFixed) {
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
    map_strcat(line, 256, buffer);
  };
  return MAP_SAFE;
};


/**
 *
 */
MAP_ERROR_CODE write_node_information_to_summary_file(FILE* file, ModelData* dataObj, char* map_msg, MAP_ERROR_CODE* ierr)
{
  int num =0;
  int i = 0;
  int cnt = 0;
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
  Node* nodeIter = NULL;  
  const int FOUR = 4;
  unsigned int numberOfNodes = 0;
  unsigned int col = 0;
  MAP_ERROR_CODE success = MAP_SAFE;

  numberOfNodes = list_size(&dataObj->node);

  do {
    for (i=0 ; i<numberOfNodes ; i+=FOUR) {
      if (i+FOUR>numberOfNodes) {
        num = numberOfNodes-i;
      } else {
        num = FOUR;
      };
    
      for (col=i ; col<i+num ; col++) {
        nodeIter = (Node*)list_get_at(&dataObj->node, col);      
        success = write_node_header_to_summary_file(col-i, cnt, col+1, line0); CHECKERRQ(MAP_FATAL_70);
        success = write_node_type_to_summary_file(col-i, cnt, nodeIter->type, line1); CHECKERRQ(MAP_FATAL_70);
        success = write_node_x_position_to_summary_file(col-i, cnt, &nodeIter->positionPtr.x, line2); CHECKERRQ(MAP_FATAL_70);
        success = write_node_y_position_to_summary_file(col-i, cnt, &nodeIter->positionPtr.y, line3); CHECKERRQ(MAP_FATAL_70);
        success = write_node_z_position_to_summary_file(col-i, cnt, &nodeIter->positionPtr.z, line4); CHECKERRQ(MAP_FATAL_70);
        success = write_node_mass_information_to_summary_file(col-i, cnt, &nodeIter->MApplied, line5); CHECKERRQ(MAP_FATAL_70);
        success = write_node_buoyancy_information_to_summary_file(col-i, cnt, &nodeIter->BApplied, line6); CHECKERRQ(MAP_FATAL_70);
        success = write_node_x_sum_force_to_summary_file(col-i, cnt, &nodeIter->sumForcePtr.fx, line7); CHECKERRQ(MAP_FATAL_70);
        success = write_node_y_sum_force_to_summary_file(col-i, cnt, &nodeIter->sumForcePtr.fy, line8); CHECKERRQ(MAP_FATAL_70);
        success = write_node_z_sum_force_to_summary_file(col-i, cnt, &nodeIter->sumForcePtr.fz, line9); CHECKERRQ(MAP_FATAL_70);
        cnt++;
      };
      cnt=0;
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


/**
 *
 */
MAP_ERROR_CODE write_element_information_to_summary_file( FILE *file, ModelData *dataObj )
{
  unsigned int numberOfElements = 0;
  Element* elementIter = NULL;  
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

  numberOfElements = list_size(&dataObj->element);
  
  for (i=0 ; i<numberOfElements ; i++) {
    elementIter = (Element*)list_get_at(&dataObj->element, i);
    
    if (elementIter->Lu.value>0.0) {
      map_strcat(line0, 256, " ");
    };
    if (elementIter->lb.value>0.0) {
      map_strcat(line1, 256, " ");
    };
    if (*(elementIter->H.value)>0.0) {
        map_strcat(line2, 256, " ");
    };
    if (*(elementIter->V.value)>0.0) {
      map_strcat(line3, 256, " ");
    };
    if (elementIter->T.value>0.0) {
      map_strcat(line4, 256, " ");
    };
    if (elementIter->alpha.value>0.0) {
      map_strcat(line5, 256, " ");
    };
    if (elementIter->HAtAnchor.value>0.0) {
      map_strcat(line6, 256, " ");
    };
    if (elementIter->VAtAnchor.value>0.0) {
      map_strcat(line7, 256, " ");
    };
    if (elementIter->TAtAnchor.value>0.0) {
      map_strcat(line8, 256, " ");
    };
    if (elementIter->alphaAtAnchor.value>0.0) {
      map_strcat(line9, 256, " ");
    };
        
    switch (elementIter->convergeReason) {
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
    fprintf(file, "Material        |  %s\n", elementIter->lineProperty->label);
    fprintf(file, "Lu        [m]   | %s%1.3f\n", line0, elementIter->Lu.value);
    fprintf(file, "Lb        [m]   | %s%1.3f\n", line1, elementIter->lb.value); 
    fprintf(file, "H         [N]   | %s%1.3f\n", line2, *(elementIter->H.value));
    fprintf(file, "V         [N]   | %s%1.3f\n", line3, *(elementIter->V.value));
    fprintf(file, "T         [N]   | %s%1.3f\n", line4, elementIter->T.value);
    fprintf(file, "Alpha     [deg] | %s%1.3f\n", line5, elementIter->alpha.value*RAD2DEG); 
    fprintf(file, "HAnch     [N]   | %s%1.3f\n", line6, elementIter->HAtAnchor.value);
    fprintf(file, "VAnch     [N]   | %s%1.3f\n", line7, elementIter->VAtAnchor.value);
    fprintf(file, "TAnch     [N]   | %s%1.3f\n", line8, elementIter->TAtAnchor.value);
    fprintf(file, "AlphaAnch [deg] | %s%1.3f\n", line9, elementIter->alphaAtAnchor.value*RAD2DEG);
    fprintf(file, "L^2-Norm        |  %1.7g\n", elementIter->residualNorm);
    fprintf(file, "Function Evals  |  %d\n", elementIter->numFuncEvals);
    fprintf(file, "Jacobian Evals  |  %d\n", elementIter->numJacEvals);
    fprintf(file, "Term. criteria  |  %d : %s\n", elementIter->convergeReason, line10);
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


/**
 *
 */
MAP_ERROR_CODE write_expanded_input_file_to_summary_file(FILE* file, InitializationData* initData)
{
  int i = 0;

  fprintf(file, "---------------------- LINE DICTIONARY ---------------------------------------\n");
  fprintf(file, "LineType  Diam      MassDenInAir   EA            CB\n");
  fprintf(file, "(-)       (m)       (kg/m)        (kN)           (-)\n");  
  for (i=0 ; i<initData->librarySize ; i++) {
    fprintf(file, "%s", initData->libraryInputString[i]);
    if (initData->libraryInputString[i][strlen(initData->libraryInputString[i])-1]!='\n') {
      fprintf(file, "\n" );
    };
  };

  fprintf(file, "---------------------- NODE PROPERTIES ---------------------------------------\n");
  fprintf(file, "Node  Type       X       Y       Z      M     B     FX    FY    FZ\n");
  fprintf(file, "(-)   (-)       (m)     (m)     (m)    (kg)  (mˆ3) (kN)  (kN)  (kN)\n");
  for (i=0 ; i<initData->sizeOfFullNodeString ; i++) {    
    fprintf(file, "%s", initData->expandedNodeInputString[i]);
    if (initData->expandedNodeInputString[i][strlen(initData->expandedNodeInputString[i])-1]!='\n') {
      fprintf(file, "\n" );
    };
  };

  fprintf(file, "---------------------- LINEP ROPERTIES ---------------------------------------\n");
  fprintf(file, "Element  LineType  UnstrLen  NodeAnch  NodeFair  Flags\n");
  fprintf(file, "(-)      (-)       (m)       (-)       (-)       (-)\n");
  for (i=0 ; i<initData->sizeOfFullElementString ; i++) {    
    fprintf(file, "%s", initData->expandedElementInputString[i]);
    if (initData->expandedElementInputString[i][strlen(initData->expandedElementInputString[i])-1]!='\n') {
      fprintf(file, "\n" );
    };
  };

  fprintf(file, "---------------------- SOLVER OPTIONS-----------------------------------------\n");
  fprintf(file, "Option \n");
  fprintf(file, "(-) \n");
  for (i=0 ; i<initData->solverOptionsSize ; i++) {    
    fprintf(file, "%s", initData->solverOptionsString[i]);
    if (initData->solverOptionsString[i][strlen(initData->solverOptionsString[i])-1]!='\n') {
      fprintf(file, "\n");
    };
  };
  return MAP_SAFE;
};
