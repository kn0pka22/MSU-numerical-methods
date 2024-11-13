OPTS = -O3 -g -Wall -Wextra -Wunused -Werror -pedantic -std=c++11
all: prog
prog: main.o solution.o pvalue.o
	g++ $(OPTS) main.o solution.o pvalue.o -o prog
main.o: main.cpp solution.hpp 
	g++ -c $(OPTS) main.cpp
solution.o: solution.cpp solution.hpp
	g++ -c $(OPTS) solution.cpp
p.o: p.cpp
	g++ -c $(OPTS)  pvalue.cpp
clean:
	rm -f main.o solution.o prog

