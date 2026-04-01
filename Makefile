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

test_mul: src/bigmul.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o tests/test_mul.c
	$(CC) $(CFLAGS) $^ -o test_mul $(LDFLAGS)
	./test_mul

test_primes: src/bigmul.o src/bigmod.o src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c tests/test_primes.c
	$(CC) $(CFLAGS) $^ -o test_primes $(LDFLAGS)
	./test_primes

clean:
	rm -f $(OBJS) $(TARGET)

