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

PLATFORM = $(shell uname -s)
CMINPACK_DIR = cminpack
BSTRING_DIR = bstring
SIMCLIST_DIR = simclist

VPATH = cminpack:simclist:bstring


ifeq ($(PLATFORM),Darwin)
CC_TOOLS  = clang
CFLAGS    = -g -O1 -fsanitize=address -fno-omit-frame-pointer -fPIC -D DEBUG -Icminpack -Isimclist
LDFLAGS   = -g -fsanitize=address -fno-omit-frame-pointer -dynamiclib
else ifeq ($(PLATFORM),Linux)
CC_TOOLS  = gcc
#CFLAGS    = -g -O1 -fsanitize=address -fno-omit-frame-pointer -fPIC -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" -Icminpack -Isimclist
CFLAGS    = -g -fPIC -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" 
LDFLAGS   = -g -shared 
endif





GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
DEBUG     = 
OBJ       = lmder.o \
	dpmpar.o \
	lmpar.o \
	qrsolv.o \
	enorm.o \
	enorm_u.o \
	qrfac.o \
	simclist.o \
	bstrlib.o \
	bstraux.o \
	supportfiles.o \
	freedata.o \
	initialization.o \
	mapinit.o \
	nwtcfunctions.o \
	maperror.o \
	MAP_Types.o \
	lineroutines.o \
	numeric.o

all : $(OBJ)
	$(CC_TOOLS) $(LDFLAGS) -o libmap-1.00.01.so $(DEBUG) $(OBJ) -lm
	@echo $(CFLAGS)

.c.o :
	$(CC_TOOLS) -c $(CFLAGS) $<

clean:
	rm -rf *.so *.o *~

memcheck:	
	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py	
