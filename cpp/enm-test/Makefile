CC = g++
LDLIBS = -I../ -L/usr/local/lib/ -lnlopt
CFLAGS = -std=c++11 -O3 -c 
LFLAGS = -std=c++11 -O3
BINARY = enm-test

enm-test: enm-test.o simulation.o estimation.o stats.o proj.o progress.o enm-test.h
	$(CC) $(LFLAGS) -o $(BINARY) enm-test.o simulation.o estimation.o stats.o proj.o progress.o $(LDLIBS)

enm-test.o: enm-test.cpp simulation.h 
	$(CC) $(CFLAGS) enm-test.cpp $(LDLIBS)

simulation.o: simulation.cpp simulation.h
	$(CC) $(CFLAGS) simulation.cpp $(LDLIBS)

estimation.o: estimation.cpp estimation.h
	$(CC) $(CFLAGS) estimation.cpp $(LDLIBS)

stats.o: stats.cpp stats.h
	$(CC) $(CFLAGS) stats.cpp $(LDLIBS)

proj.o: proj.cpp proj.h
	$(CC) $(CFLAGS) proj.cpp $(LDLIBS)

progress.o: progress.cpp progress.h
	$(CC) $(CFLAGS) progress.cpp $(LDLIBS)


clean:
	rm -f $(BINARY) *.o 
