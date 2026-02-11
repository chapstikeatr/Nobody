# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17

# Target executable name
TARGET = nbody

# Source files
SRC = nbody.cpp

# Default target
all: $(TARGET)

# Build rule
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Clean rule
clean:
	rm -f $(TARGET)
