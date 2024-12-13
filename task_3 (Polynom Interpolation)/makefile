CXX = g++
CXXFLAGS = -c -Wall
TARGET = task_3
ARGS = 20

.PHONY: all clean run

all: $(TARGET)

$(TARGET): main.o task_3.o linearSystemSolver.o
	$(CXX) $^ -o $@

main.o: main.cpp task_3.hpp linearSystemSolver.hpp
	$(CXX) $(CXXFLAGS) main.cpp

task_3.o: task_3.cpp task_3.hpp linearSystemSolver.hpp
	$(CXX) $(CXXFLAGS) task_3.cpp

linearSystemSolver.o: linearSystemSolver.cpp linearSystemSolver.hpp
	$(CXX) $(CXXFLAGS) linearSystemSolver.cpp

run: $(TARGET)
	./$(TARGET) $(ARGS) 0 out.txt

clean:
	rm -rf *.o $(TARGET)
