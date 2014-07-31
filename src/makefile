# =========== MAP Makefile        =============================================
#	     
#  Copyright Jan. 21, 2014
#  Author: Marco D. Masciola
# 	
#  MAP is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
#
#  MAP is distributed in the hope that it will be useful, but WITHOUT ANY 
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
#  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along 
#  with MAP. If not, see:
#  
#  < http://www.gnu.org/licenses/>
#  ============================================================================

CMINPACK_DIR = cminpack
SIMCLIST_DIR = simclist

VPATH=cminpack:simclist

CC_TOOLS  = gcc 
CFLAGS    = -g -fPIC -D DEBUG -Icminpack -Isimclist
#CFLAGS    = -g -O2 -fPIC -Wstrict-aliasing
LDFLAGS   = -shared 
DEBUG     = 
OBJ       = lmder.o \
	dpmpar.o \
	lmpar.o \
	qrsolv.o \
	enorm.o \
	enorm_u.o \
	qrfac.o \
	simclist.o \
	nwtcfunctions.o \
	fortdatamanager.o \
	mapinit.o \
	maperror.o \
	freedata.o \
	initialization.o \
	supportfiles.o \
	lineupdateroutines.o \
	numeric.o \
	outputstream.o \
	linegenerator.o \
	MAP_Types.o \
	fortroutines.o \
	pyroutines.o

all : $(OBJ)
	$(CC_TOOLS) $(LDFLAGS) -o map.so $(DEBUG) $(OBJ)

.c.o :
	$(CC_TOOLS) -c $(CFLAGS) $<

clean:
	rm -rf *.so *.o *~

memcheck:	
	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py	



# CC_TOOLS  = gcc 
# CFLAGS    = -g -fPIC -D DEBUG 
# #CFLAGS    = -g -O2 -fPIC -Wstrict-aliasing
# LDFLAGS   = -shared 
# DEBUG     = 
# OBJ       = $(CMINPACK_DIR)/lmder.o \
# 	$(CMINPACK_DIR)/dpmpar.o \
# 	$(CMINPACK_DIR)/lmpar.o \
# 	$(CMINPACK_DIR)/qrsolv.o \
# 	$(CMINPACK_DIR)/enorm.o \
# 	$(CMINPACK_DIR)/enorm_u.o \
# 	$(CMINPACK_DIR)/qrfac.o \
# 	$(SIMCLIST_DIR)/simclist.o \
# 	nwtcfunctions.o \
# 	fortdatamanager.o \
# 	mapinit.o \
# 	maperror.o \
# 	freedata.o \
# 	initialization.o \
# 	supportfiles.o \
# 	lineupdateroutines.o \
# 	numeric.o \
# 	outputstream.o \
# 	linegenerator.o \
# 	MAP_Types.o \
# 	fortroutines.o \
# 	pyroutines.o
# 
# all : $(OBJ)
# 	$(CC_TOOLS) $(LDFLAGS) -o map.so $(DEBUG) $(OBJ)
# 
# .c.o :
# 	$(CC_TOOLS) -c $(CFLAGS) $<
# 
# clean:
# 	rm -rf *.so *.o *~
# 
# memcheck:	
# 	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py	
# 	#valgrind --tool=callgrind python ./main.py		
# 	#valgrind --tool=memcheck --log-file="memcheck.out" python -E -tt ./main.py		
