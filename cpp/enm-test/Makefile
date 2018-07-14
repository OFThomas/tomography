CC = g++
CFLAGS = -std=c++11 -c -I../
LFLAGS = -std=c++11 -I../
IFLAGS = 

enm-test: enm-test.o 
	$(CC) $(LFLAGS) -o enm-test enm-test.o
	rm -f *.o *~

enm-test.o: enm-test.cpp 
	$(CC) $(CFLAGS) enm-test.cpp
