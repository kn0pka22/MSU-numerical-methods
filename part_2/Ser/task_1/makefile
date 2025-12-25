# OPTS = -O3 -g -Wall -Wextra -Wunused -Werror -pedantic -std=c++11
OPTS = -O3 -g 


all: prog

prog: main.o solution.o
	g++ $(OPTS) main.o solution.o -o prog
main.o: main.cpp solution.hpp
	g++ -c $(OPTS) main.cpp
solution.o: solution.cpp solution.hpp
	g++ -c $(OPTS) solution.cpp
clean:
	rm -f main.o solution.o  out.txt p.txt prog

