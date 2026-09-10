# Compiler
CC = gcc
AS = nasm
AR = ar
# Flags 
CFLAGS = -Iinclude -Wall -Wextra -g -O3 -fopenmp  -funroll-loops -fopenmp 
LDFLAGS = -fopenmp -lm -flto -lgmp
ASFLAGS = -f elf64

# directories 
SRC_DIR = src
INC_DIR = include
BUILD_DIR = build
CMD_DIR = cmd
TEST_DIR = tests

SRCS = $(shell find $(SRC_DIR) -name '*.c')
ASMS = $(shell find $(SRC_DIR) -name '*.asm' -o -name '*.s')

OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/$(SRC_DIR)/%.o, $(SRCS))
ASM_OBJS = $(patsubst $(SRC_DIR)/%.asm, $(BUILD_DIR)/$(SRC_DIR)/%.o, $(filter %.asm, $(ASMS))) \
           $(patsubst $(SRC_DIR)/%.s, $(BUILD_DIR)/$(SRC_DIR)/%.o, $(filter %.s, $(ASMS)))

LIB = $(BUILD_DIR)/libbignum.a
TARGET = $(BUILD_DIR)/bignum

TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)
TEST_BINS = $(patsubst $(TEST_DIR)/%.c, $(BUILD_DIR)/%, $(TEST_SRCS))

#.PHONY: all clean tests $(TEST_BINS)
.PHONY: $(patsubst $(BUILD_DIR)/%, test_%, $(TEST_BINS))
.SECONDARY: $(OBJS) $(ASM_OBJS) $(TEST_BINS)
.PRECIOUS: $(BUILD_DIR)/tests/%.o $(BUILD_DIR)/%


all: $(TARGET)

# 1. Compile the static library (Contains all logic from src/)
$(LIB): $(OBJS) $(ASM_OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

# 2. Compile the Main Application
$(BUILD_DIR)/$(CMD_DIR)/main.o: $(CMD_DIR)/main.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(TARGET): $(BUILD_DIR)/$(CMD_DIR)/main.o $(LIB)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# 3. Generic rules for compiling C and ASM files into build/
$(BUILD_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.asm
	@mkdir -p $(@D)
	$(AS) $(ASFLAGS) $< -o $@

$(BUILD_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.s
	@mkdir -p $(@D)
	$(AS) $(ASFLAGS) $< -o $@

# 4. Rules for compiling and linking Tests
$(BUILD_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/test_%: $(BUILD_DIR)/$(TEST_DIR)/test_%.o $(LIB)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_test_$*)

# 5. Convenience targets to run individual tests (e.g., `make test_mul`)
test_%: $(BUILD_DIR)/test_%
	./$<

# 6. Build all tests without running them
tests: $(TEST_BINS)

# Clean up all build artifacts
clean:
	rm -rf $(BUILD_DIR)
