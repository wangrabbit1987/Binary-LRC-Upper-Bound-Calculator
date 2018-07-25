
start: upperbound 

upperbound: main.cpp
	g++ -o main.o -c main.cpp 
	g++ -o upperbound main.o

clean:
	rm -rf *.o *.exe *.stackdump upperbound

run: upperbound
	./upperbound
