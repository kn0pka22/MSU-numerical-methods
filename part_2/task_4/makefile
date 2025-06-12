CXX = g++
# CXXFLAGS = -Wall -Wextra -std=c++11
TARGET = proga

SRCS = main.cpp funcs.cpp
OBJS = $(SRCS:.cpp=.o)


all: $(TARGET)


$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)
%.o: %.cpp funcs.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)