CC = gcc
#CCFLAGS = -g -O0 
CFLAGS = -g -O2 -Wstrict-aliasing -Wuninitialized  
OBJ = test.o

all : $(OBJ)
	$(CC) $(CCFLAGS) -o c_driver $(OBJ) -L. ../src/libmap-1.00.01.so

clean:		
	rm -f *.o c_driver *.*~ *.*# makefile~ *.mod *.valgrind


memcheck:	
	valgrind  --tool=exp-sgcheck  ./c_driver
#	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py	
