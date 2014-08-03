/***************************************************************************
 *   Copyright (C) 2014 mdm                                                *
 *   marco.masciola@gmail.com                                              *
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


#ifndef _FREE_DATA_H
#define _FREE_DATA_H


/**
 * @brief     Deallocates the memory space for the init structure. Should be called immediately after map_init()
 * @param     init MAP-native initialization data structure. This is distinct from the FAST-framework data structure 
 * @param     map_msg error string
 * @param     ierr error code
 * @return    MAP_SAFE if it completes successfully
 * @see       {@link map_init()}
 */
MAP_EXTERNCALL int free_init_data (InitializationData* init, char* map_msg, MAP_ERROR_CODE* ierr);

/**
 * @brief     Deallocates all elements. Function loops through the elemenet link list and frees allocated data. Pointers
 *            are nullified.  
 * @param     element the element link list
 * @return    MAP_SAFE if it completes successfully
 * @see       {@link Element_t()}
 */
MAP_ERROR_CODE free_element(list_t *restrict element);

/**
 * @brief     Deallocates all nodes. Function loops through the elemenet link list and frees allocated data. Pointers
 *            are nullified.  
 * @param     node the node link list
 * @return    MAP_SAFE if it completes successfully
 * @see       {@link Element_t()}
 */
MAP_ERROR_CODE free_node(list_t *restrict node);


void MAP_InitInput_Delete(InitializationData* init_data);


void MAP_OtherState_Delete(ModelData* model_data);


MAP_ERROR_CODE free_vessel(Vessel* floater);


#endif // _FREE_DATA_H
