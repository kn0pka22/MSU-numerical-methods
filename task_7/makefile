CFLAGS = -c -Wall
TARGET = task_7
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o task_7.o 
	g++ main.o task_7.o  -o $(TARGET)

main.o: main.cpp task_7.hpp
	g++ $(CFLAGS) main.cpp

task_7.o: task_7.cpp task_7.hpp
	g++ $(CFLAGS) task_7.cpp

# run: $(TARGET)
# 	./$(TARGET) $(ARGS) 0 out.txt
#jupyter nbconvert --to notebook --execute graph.ipynb --output graph_output.ipynb

# graphs: graph.ipynb
# 	python3 graph.ipynb

clean:
	rm -rf *.o $(TARGET) 
