CFLAGS = -c -Wall
TARGET = task_2
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o task_2.o
	g++ main.o task_2.o -o $(TARGET)  

main.o: main.cpp task_2.hpp
	g++ $(CFLAGS) main.cpp

task_1.o: task_2.cpp task_2.hpp
	g++ $(CFLAGS) task_1.cpp

# pcalculate.o: pcalculate.cpp
# 	g++ $(CFLAGS) pcalculate.cpp

run: $(TARGET)
	./$(TARGET) $(ARGS) 0 out.txt

# graphs: graph.ipynb
# 	python3 graph.ipynb

clean:
	rm -rf *.o $(TARGET) 
