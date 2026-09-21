CC = gcc
CFLAGS = -O3 -fopenmp -Wall -Wextra
LDFLAGS = -lm -fopenmp

TARGET = openmp
SRCS = openmp.c

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
