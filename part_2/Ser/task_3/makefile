CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = proga

SRCS = main.cpp funcs.cpp
OBJS = $(SRCS:.cpp=.o)


all: $(TARGET)

# Правило для сборки исполняемого файла:
# $(TARGET) зависит от $(OBJS) - объектных файлов
# Команда: $(CXX) - компилятор g++
#          $(CXXFLAGS) - флаги компиляции
#          -o $(TARGET) - имя выходного файла
#          $(OBJS) - список объектных файлов
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)
#          -c $< - компилировать исходный файл (первая зависимость)
#          -o $@ - выходной файл (имя цели)
%.o: %.cpp funcs.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)