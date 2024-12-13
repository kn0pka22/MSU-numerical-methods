CFLAGS = -c -Wall
TARGET = task_5
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o task_5.o 
	g++ main.o task_5.o  -o $(TARGET)

main.o: main.cpp task_5.hpp
	g++ $(CFLAGS) main.cpp

task_5.o: task_5.cpp task_5.hpp
	g++ $(CFLAGS) task_5.cpp

run: $(TARGET)
	./$(TARGET) $(ARGS) 0 out.txt
#jupyter nbconvert --to notebook --execute graph.ipynb --output graph_output.ipynb

graphs: graph.ipynb
	python3 graph.ipynb

clean:
	rm -rf *.o $(TARGET) graph_output.ipynb
