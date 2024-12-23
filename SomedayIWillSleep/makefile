CFLAGS = -c -Wall
TARGET = prog
ARGS = 20



all: $(TARGET)

$(TARGET): main.o sol.o 
	g++ main.o sol.o  -o $(TARGET)

main.o: main.cpp sol.hpp
	g++ $(CFLAGS) main.cpp

sol.o: sol.cpp sol.hpp
	g++ $(CFLAGS) sol.cpp

clean:
	rm -rf *.o $(TARGET) 
