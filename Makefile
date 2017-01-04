CC = g++
CFLAGS = -O3 -c
LFLAGS = -O3 -o QuantumMetrology
INCLUDE = -I /home/jake/Documents/EIGEN
OBJS = main.o MeritFunction.o BFGS_Optimization.o LOTransform.o

QuantumMetrology: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS)

main.o: main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) main.cpp

MeritFunction.o: MeritFunction.cpp
	$(CC) $(CFLAGS) $(INCLUDE) MeritFunction.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	$(CC) $(CFLAGS) $(INCLUDE) BFGS_Optimization.cpp

LOTransform.o: LOTransform.cpp
	$(CC) $(CFLAGS) $(INCLUDE) LOTransform.cpp

clean:
	rm *.o QuantumMetrology
