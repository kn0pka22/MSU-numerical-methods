CFLAGS = -c -Wall
TARGET = task_1
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o task_1.o pcalculate.o
	g++ main.o task_1.o pcalculate.o -o $(TARGET)

main.o: main.cpp task_1.hpp
	g++ $(CFLAGS) main.cpp

task_1.o: task_1.cpp task_1.hpp
	g++ $(CFLAGS) task_1.cpp

pcalculate.o: pcalculate.cpp
	g++ $(CFLAGS) pcalculate.cpp

run: $(TARGET)
	./$(TARGET) $(ARGS) 0 out.txt
#jupyter nbconvert --to notebook --execute graph.ipynb --output graph_output.ipynb

graphs: graph.ipynb
	python3 graph.ipynb

clean:
	rm -rf *.o $(TARGET) graph_output.ipynb
