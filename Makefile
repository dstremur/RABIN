# Compiler
CC = gcc
CFLAGS = -Iinclude -Wall -Wextra -g -O3 -fopenmp -flto -march=native 

LDFLAGS = -fopenmp -lm -flto

TARGET = bignum

 SRCS = $(wildcard src/*.c)

# Convert to object files
OBJS = $(SRCS:.c=.o)

# The default rule (runs when you just type 'make')
all: $(TARGET)

# Link the object files into the final executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) 

$(OBJS): include/bignum.h

# Rule to compile .c files into .o files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

