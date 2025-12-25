CFLAGS = -c -Wall -std=c++11
TARGET = proga
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o funcs.o 
	g++ main.o funcs.o -o $(TARGET)

main.o: main.cpp funcs.hpp
	g++ $(CFLAGS) main.cpp

funcs.o: funcs.cpp funcs.hpp
	g++ $(CFLAGS) funcs.cpp

clean:
	rm -rf *.o $(TARGET)

run: $(TARGET)
	./$(TARGET) $(ARGS)