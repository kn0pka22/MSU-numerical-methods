CFLAGS = -c -Wall
TARGET = task_7
ARGS = 20

.PHONY: all clean run graphs

all: $(TARGET)

$(TARGET): main.o task_7.o 
	g++ main.o task_7.o -o $(TARGET)

main.o: main.cpp task_7.hpp
	g++ $(CFLAGS) main.cpp

task_7.o: task_7.cpp task_7.hpp
	g++ $(CFLAGS) task_7.cpp  # Уберите task_7.hpp из этой строки

clean:
	rm -rf *.o $(TARGET)