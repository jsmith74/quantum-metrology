CC = g++
CFLAGS = -O3 -c
LFLAGS = -O3 -o QuantumMetrology
OMPFLAGS = -fopenmp
INCLUDE = -I /home/jake/Documents/EIGEN
OBJS = main.o MeritFunction.o BFGS_Optimization.o LOTransform.o MZIMeas.o BranchMeasStruct.o Integration.o

QuantumMetrology: $(OBJS)
	$(CC) $(OMPFLAGS) $(LFLAGS) $(OBJS)

main.o: main.cpp
	$(CC) $(OMPFLAGS) $(CFLAGS) $(INCLUDE) main.cpp

MeritFunction.o: MeritFunction.cpp
	$(CC) $(CFLAGS) $(INCLUDE) MeritFunction.cpp

BranchMeasStruct.o: BranchMeasStruct.cpp
	$(CC) $(CFLAGS) $(INCLUDE) BranchMeasStruct.cpp

MZIMeas.o: MZIMeas.cpp
	$(CC) $(CFLAGS) $(INCLUDE) MZIMeas.cpp

Integration.o: Integration.cpp
	$(CC) $(CFLAGS) $(INCLUDE) Integration.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	$(CC) $(CFLAGS) $(INCLUDE) BFGS_Optimization.cpp

LOTransform.o: LOTransform.cpp
	$(CC) $(CFLAGS) $(INCLUDE) LOTransform.cpp

clean:
	rm *.o QuantumMetrology *.out *.dat
