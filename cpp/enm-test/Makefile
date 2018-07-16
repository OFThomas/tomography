CC = g++
CFLAGS = -std=c++11 -O3 -c -I../ 
LFLAGS = -std=c++11 -O3 -I../


enm-test: enm-test.o simulation.o estimation.o stats.o proj.o enm-test.h
	$(CC) $(LFLAGS) -o enm-test enm-test.o simulation.o estimation.o stats.o proj.o

enm-test.o: enm-test.cpp simulation.h 
	$(CC) $(CFLAGS) enm-test.cpp

simulation.o: simulation.cpp simulation.h
	$(CC) $(CFLAGS) simulation.cpp

estimation.o: estimation.cpp estimation.h
	$(CC) $(CFLAGS) estimation.cpp

stats.o: stats.cpp stats.h
	$(CC) $(CFLAGS) stats.cpp

proj.o: proj.cpp proj.h
	$(CC) $(CFLAGS) proj.cpp
