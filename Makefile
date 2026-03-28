# Compiler
CC = gcc
CFLAGS = -Iinclude -Wall -Wextra -g -O3 -fopenmp -flto -march=native 

LDFLAGS = -fopenmp -lm -flto

TARGET = bignum

SRCS = $(wildcard src/*.c)

OBJS = $(SRCS:.c=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) 

$(OBJS): include/bignum.h

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

