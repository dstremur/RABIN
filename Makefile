# Compiler
CC = gcc
AS = nasm
AR = ar
# Flags
CFLAGS = -Iinclude -Wall -Wextra -g -O3 -fopenmp  -funroll-loops -fopenmp 
LDFLAGS = -fopenmp -lm -flto -lgmp
ASFLAGS = -f elf64

# Build-time AVX-512 detection: enable -mavx512f only when the compiler
 #accepts the flag and the build machine's CPU exposes avx512f.
 HAVE_AVX512 := $(shell grep -qw avx512f /proc/cpuinfo 2>/dev/null && \
     printf 'int main(void){return 0;}' | $(CC) -mavx512f -x c - -o /dev/null \
     2>/dev/null && echo yes)
 ifeq ($(HAVE_AVX512),yes)
 CFLAGS += -mavx512f
 endif

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

# test_flint requires libflint; skip it from the build when it is absent
FLINT_OK := $(shell echo '#include <flint/flint.h>' | $(CC) -E -x c - \
    -o /dev/null 2>/dev/null && echo yes)
ifneq ($(FLINT_OK),yes)
TEST_SRCS := $(filter-out $(TEST_DIR)/test_flint.c, $(TEST_SRCS))
endif

TEST_BINS = $(patsubst $(TEST_DIR)/%.c, $(BUILD_DIR)/%, $(TEST_SRCS))
EXTRA_test_flint = -lflint

#.PHONY: all clean tests $(TEST_BINS)
.PHONY: test-valgrind $(patsubst $(BUILD_DIR)/%, test_%, $(TEST_BINS))
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

TESTS = test_add test_add_u64 test_bpsw test_cmp test_div test_div_gcd \
        test_divmod test_divmod_u64 test_field test_gcd \
        test_gcd_extended test_gmp test_isqrt test_log_2 \
        test_lshift test_mod test_mod_inverse test_mod_u64

test_complete: $(TESTS)
	@echo "All specified tests completed successfully!"

# 8. Run one test (`make test-valgrind TEST=test_mul`) or the fast core set
#    under valgrind (definite/possible leaks + error abort). Valgrind cannot
#    emulate AVX-512, so tests run from a parallel no-AVX-512 build (build/vg).
VALGRIND ?= valgrind
VALGRIND_FLAGS = --error-exitcode=99 --leak-check=full \
                 --show-leak-kinds=definite,possible
VALGRIND_TESTS = test_add test_add_u64 test_cmp test_div test_divmod \
                 test_divmod_u64 test_field test_gcd test_gcd_extended \
                 test_isqrt test_log_2 test_lshift test_mod \
                 test_mod_inverse test_mod_u64

VG_DIR = $(BUILD_DIR)/vg
VG_CFLAGS = $(CFLAGS) -mno-avx512f
VG_OBJS = $(patsubst $(BUILD_DIR)/%.o, $(VG_DIR)/%.o, $(OBJS))
VG_ASM_OBJS = $(patsubst $(BUILD_DIR)/%.o, $(VG_DIR)/%.o, $(ASM_OBJS))
VG_LIB = $(VG_DIR)/libbignum.a

$(VG_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(VG_CFLAGS) -c $< -o $@

$(VG_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.asm
	@mkdir -p $(@D)
	$(AS) $(ASFLAGS) $< -o $@

$(VG_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.s
	@mkdir -p $(@D)
	$(AS) $(ASFLAGS) $< -o $@

$(VG_LIB): $(VG_OBJS) $(VG_ASM_OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

$(VG_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(VG_CFLAGS) -c $< -o $@

$(VG_DIR)/test_%: $(VG_DIR)/$(TEST_DIR)/test_%.o $(VG_LIB)
	@mkdir -p $(@D)
	$(CC) $(VG_CFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_test_$*)

test-valgrind:
	@target="$(or $(TEST),$(VALGRIND_TESTS))"; \
	for t in $$target; do \
		echo "== valgrind $$t =="; \
		$(MAKE) $(VG_DIR)/$$t || exit 1; \
		$(VALGRIND) $(VALGRIND_FLAGS) ./$(VG_DIR)/$$t || exit 1; \
	done

# Clean up all build artifacts
clean:
	rm -rf $(BUILD_DIR)
