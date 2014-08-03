CC = gcc
CCFLAGS = -g -O0 
OBJ = test.o

all : $(OBJ)
	$(CC) $(CCFLAGS) -o c_driver $(OBJ) -L. ../src/libmap-1.00.01.so

clean:		
	rm -f *.o c_driver *.*~ *.*# makefile~ *.mod *.valgrind


memcheck:	
	valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./c_driver
