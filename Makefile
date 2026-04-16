# Compiler
CC = gcc
AS = nasm
CFLAGS = -Iinclude -Wall -Wextra -g -O3 -fopenmp -flto=auto -march=native 

ASFLAGS = -f elf64

LDFLAGS = -fopenmp -lm -flto

TARGET = bignum

SRCS = $(wildcard src/*.c)
ASMS = $(wildcard src/*.asm)

OBJS = $(SRCS:.c=.o) src/bn_mul_inner.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS) 

$(OBJS): include/bignum.h

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule for Assembly files
%.o: %.asm
	$(AS) $(ASFLAGS) $< -o $@

test_mul: src/bigmul.o src/bn_mul_inner.o src/bigrand.c src/bn_avx512.c src/bigprime.c src/bignum.o src/bigadd.o src/bigsub.o src/bigdiv.c src/bigshift.o tests/test_mul.c
	$(CC) $(CFLAGS) $^ -o test_mul $(LDFLAGS)
	./test_mul

test_div: src/bigmul.o src/bn_mul_inner.o src/bignum.o src/bigprime.c src/bigadd.o src/bigsub.o src/bigdiv.c src/bigshift.o tests/test_div.c src/bigrand.c
	$(CC) $(CFLAGS) $^ -o test_div $(LDFLAGS)
	./test_div

test_primes: src/bigmul.o src/bigrand.c src/bigmod.o src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c src/biglucas.c src/bigprime.c src/bigdiv.c src/bigmath.c tests/test_primes.c src/bn_mul_inner.o
	$(CC) $(CFLAGS) $^ -o test_primes $(LDFLAGS)
	./test_primes

test_primes_parallel: src/bigmul.o src/bn_mul_inner.o src/bigrand.c src/bigmod.o src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c src/biglucas.c src/bigprime.c src/bigdiv.c src/bigmath.c tests/test_primes_parallel.c 
	$(CC) $(CFLAGS) $^ -o test_primes_p $(LDFLAGS)
	./test_primes_p

test_lucas: src/bigmul.o src/bigmod.o src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c src/biglucas.c tests/test_lucas.c
	$(CC) $(CFLAGS) $^ -o test_lucas $(LDFLAGS)
	./test_lucas


test_bpsw: src/bigmul.o src/bigmod.o src/bn_mul_inner.o src/bigrand.c src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c src/biglucas.c src/bigprime.c src/bigdiv.c src/bigmath.c tests/test_bpsw.c
	$(CC) $(CFLAGS) $^ -o test_bpsw $(LDFLAGS)
	./test_bpsw

test_openssl: src/bigmul.o src/bigmod.o src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c src/biglucas.c src/bigprime.c src/bigdiv.c src/bigmath.c src/bigrand.c tests/test_openssl.c src/bn_mul_inner.o
	$(CC) $(CFLAGS) $^ -o test_openssl $(LDFLAGS)
	./test_openssl

test_matrix: src/bigmul.o src/bigrand.c src/bigmod.o src/bigmont.o src/bigexp.o src/bignum.o src/bigadd.o src/bigsub.o src/bigshift.o src/bigrabin.c src/bigmatrix.c src/biglucas.c src/bigprime.c src/bigdiv.c src/bigmath.c tests/test_matrix.c src/bn_mul_inner.o
	$(CC) $(CFLAGS) $^ -o test_matrix $(LDFLAGS)
	./test_matrix
clean:
	rm -f $(OBJS) $(TARGET)

