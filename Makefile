CC = gcc
CFLAGS = -O3 -march=native -Wall
LDFLAGS = -lm
TARGET = tsptw

all: $(TARGET)

$(TARGET): tsptw.c
	$(CC) $(CFLAGS) -o $(TARGET) tsptw.c $(LDFLAGS)

clean:
	del /Q $(TARGET).exe 2>nul

.PHONY: all clean
