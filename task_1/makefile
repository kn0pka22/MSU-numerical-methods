CFLAGS = -c -Wall
TARGET = pcalculate.cpp
ARGS = 10 
#.PHONY: all clean

all: task_1 graphs

task_1: main.o task_1.o pcalculate.o
	g++ main.o task_1.o pcalculate.o -o task_1

main.o: main.cpp task_1.hpp
	g++ $(CFLAGS) main.cpp

task_1.o: task_1.cpp task_1.hpp
	g++ $(CFLAGS) task_1.cpp

pcalculate.o: pcalculate.cpp
	g++ $(CFLAGS) pcalculate.cpp

graphs: graph.ipynb
	python3 graph.ipynb

clean:
	rm -rf *.o task_1

# run: $(TARGET)
# 	./$(TARGET) $(ARGS)