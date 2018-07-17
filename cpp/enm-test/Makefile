CC = g++
CFLAGS = -std=c++11 -O3 -c -I../ 
LFLAGS = -std=c++11 -O3 -I../
BINARY = enm-test

enm-test: enm-test.o simulation.o estimation.o stats.o proj.o progress.o enm-test.h
	$(CC) $(LFLAGS) -o $(BINARY) enm-test.o simulation.o estimation.o stats.o proj.o progress.o

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

progress.o: progress.cpp progress.h
	$(CC) $(CFLAGS) progress.cpp


clean:
	rm -f $(BINARY) *.o 
