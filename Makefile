all: QuantumMetrology

QuantumMetrology: main.o MeritFunction.o BFGS_Optimization.o PUA.o AncillaAugment.o LOTransform.o
	g++ -O3 -o QuantumMetrology main.o MeritFunction.o BFGS_Optimization.o PUA.o AncillaAugment.o LOTransform.o

main.o: main.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN main.cpp

MeritFunction.o: MeritFunction.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN MeritFunction.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN BFGS_Optimization.cpp

PUA.o: PUA.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN PUA.cpp

AncillaAugment.o: AncillaAugment.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN AncillaAugment.cpp

LOTransform.o: LOTransform.cpp
	g++ -O3 -c -I /home/jake/Documents/EIGEN LOTransform.cpp

clean:
	rm *.o QuantumMetrology
