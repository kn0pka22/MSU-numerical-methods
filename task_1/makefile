CFLAGS = -c -Wall


task_1: main.o task_1.o 
	g++ main.o task_1.o -o task_1

main.o: main.cpp task_1.cpp task_1.hpp 
	g++ $(CFLAGS) main.cpp

task_1.o: task_1.cpp task_1.hpp
	g++ $(CFLAGS) task_1.cpp
clean:
	rm -rf *.o task_1

gr:
	chmod +x printAll.gpi
	./printAll.gpi