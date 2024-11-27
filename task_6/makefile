CFLAGS = -c #-Wall
TARGET = task_6
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o task_6.o 
	g++ main.o task_6.o  -o $(TARGET)

main.o: main.cpp task_6.hpp
	g++ $(CFLAGS) main.cpp

task_6.o: task_6.cpp task_6.hpp
	g++ $(CFLAGS) task_6.cpp

# run: $(TARGET)
# 	./$(TARGET) $(ARGS) 0 out.txt
#jupyter nbconvert --to notebook --execute graph.ipynb --output graph_output.ipynb

# graphs: graph.ipynb
# 	python3 graph.ipynb

clean:
	rm -rf *.o $(TARGET) graph_output.ipynb
