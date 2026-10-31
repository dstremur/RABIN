# Compiler
CC = gcc
AS = nasm
AR = gcc-ar
# Flags
CFLAGS = -Iinclude -Wall -Wextra -g -march=native -O3 -fopenmp -flto=auto -funroll-loops
LDFLAGS = -fopenmp -lm -flto
ASFLAGS = -f elf64


# Build-time AVX-512 detection: enable -mavx512f only when the compiler
# accepts the flag and the build machine's CPU exposes avx512f.
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

# Tests that reference <gmp.h> link -lgmp; skip them when GMP is absent so the
# library itself keeps zero external dependencies
GMP_OK := $(shell echo '#include <gmp.h>' | $(CC) -E -x c - \
    -o /dev/null 2>/dev/null && echo yes)
GMP_TEST_NAMES := $(notdir $(basename $(shell grep -l '#include.*gmp' \
    $(TEST_DIR)/*.c 2>/dev/null)))
# test_rat only uses GMP when the header is available (guarded include), so
# it must survive the auto-skip above and link -lgmp conditionally.
GMP_TEST_NAMES := $(filter-out test_rat, $(GMP_TEST_NAMES))
EXTRA_test_rat = $(if $(filter yes,$(GMP_OK)),-lgmp)
ifneq ($(GMP_OK),yes)
TEST_SRCS := $(filter-out $(addprefix $(TEST_DIR)/,$(addsuffix .c,$(GMP_TEST_NAMES))), \
    $(TEST_SRCS))
endif

# test_openssl cross-checks against the `openssl prime` CLI (no libssl link);
# skip it when the binary is absent
OPENSSL_OK := $(shell command -v openssl >/dev/null 2>&1 && echo yes)
ifneq ($(OPENSSL_OK),yes)
TEST_SRCS := $(filter-out $(TEST_DIR)/test_openssl.c, $(TEST_SRCS))
endif

TEST_BINS = $(patsubst $(TEST_DIR)/%.c, $(BUILD_DIR)/%, $(TEST_SRCS))
EXTRA_test_flint = -lflint
# test_poly_div cross-checks against fmpz_poly_gcd when FLINT headers exist
# (the suite itself gates the FLINT code on __has_include)
EXTRA_test_poly_div = $(if $(filter yes,$(FLINT_OK)),-lflint)

#.PHONY: all clean tests $(TEST_BINS)
.PHONY: test-valgrind test-san test_quick test_complete check-api tests \
        install uninstall $(patsubst $(BUILD_DIR)/%, test_%, $(TEST_BINS))
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
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_test_$*) \
	  $(if $(filter test_$*,$(GMP_TEST_NAMES)),-lgmp)

# 5. Convenience targets to run individual tests (e.g., `make test_mul`)
test_%: $(BUILD_DIR)/test_%
	./$<

# 6. Build all tests without running them
tests: $(TEST_BINS)

# 7. Test sets: TESTS is every available suite (already dep-filtered by the
#    FLINT/GMP/OpenSSL probes above); TESTS_QUICK is the fast, dependency-free
#    correctness subset suitable for CI smoke runs.
TESTS = $(notdir $(patsubst $(TEST_DIR)/%.c, %, $(TEST_SRCS)))

TESTS_QUICK = test_add test_add_u64 test_all test_bpsw test_cert test_cmp \
              test_div test_div_gcd test_divmod test_divmod_u64 test_field \
              test_gcd_extended test_isqrt test_kernel test_kron test_lll \
              test_log_2 test_logic test_lshift test_matrix test_matrix_basic \
              test_matrix_elim test_matrix_kernel test_matrixq test_mod \
              test_mod_inverse test_mod_u64 test_modexp test_mul test_mul_i64 \
              test_neg test_ntt test_poly test_poly_div test_rat test_rns \
              test_rshift \
              test_sqr test_sub test_u64 test_vector

test_quick: $(TESTS_QUICK)
	@echo "Quick test suite completed successfully!"

test_complete: $(TESTS)
	@echo "All specified tests completed successfully!"

# 8. Detect API drift: every public symbol declared in include/*.h must be
#    defined in the built library unless allowlisted in tests/api_allow.txt.
check-api: $(LIB)
	@sh tests/check_api.sh $(LIB) tests/api_allow.txt

# 9. Address+UndefinedBehavior sanitizer build of the quick suite. Override
#    the suite with `make test-san TEST=test_mul` or TESTS="$(TESTS)".
SAN_DIR = $(BUILD_DIR)/san
SAN_CFLAGS = $(CFLAGS) -g -fsanitize=address,undefined \
             -fno-omit-frame-pointer -mno-avx512f
SAN_OBJS = $(patsubst $(BUILD_DIR)/%.o, $(SAN_DIR)/%.o, $(OBJS))
SAN_ASM_OBJS = $(patsubst $(BUILD_DIR)/%.o, $(SAN_DIR)/%.o, $(ASM_OBJS))
SAN_LIB = $(SAN_DIR)/libbignum.a

$(SAN_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(SAN_CFLAGS) -c $< -o $@

$(SAN_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.asm
	@mkdir -p $(@D)
	$(AS) $(ASFLAGS) $< -o $@

$(SAN_DIR)/$(SRC_DIR)/%.o: $(SRC_DIR)/%.s
	@mkdir -p $(@D)
	$(AS) $(ASFLAGS) $< -o $@

$(SAN_LIB): $(SAN_OBJS) $(SAN_ASM_OBJS)
	@mkdir -p $(@D)
	$(AR) rcs $@ $^

$(SAN_DIR)/$(TEST_DIR)/%.o: $(TEST_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(SAN_CFLAGS) -c $< -o $@

$(SAN_DIR)/test_%: $(SAN_DIR)/$(TEST_DIR)/test_%.o $(SAN_LIB)
	@mkdir -p $(@D)
	$(CC) $(SAN_CFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_test_$*) \
	  $(if $(filter test_$*,$(GMP_TEST_NAMES)),-lgmp)

test-san:
	@target="$(or $(TEST),$(TESTS_QUICK))"; \
	for t in $$target; do \
		echo "== sanitizer $$t =="; \
		$(MAKE) $(SAN_DIR)/$$t || exit 1; \
		UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 \
		ASAN_OPTIONS=detect_leaks=1 ./$(SAN_DIR)/$$t || exit 1; \
	done

# 8. Run one test (`make test-valgrind TEST=test_mul`) or the fast core set
#    under valgrind (definite/possible leaks + error abort). Valgrind cannot
#    emulate AVX-512, so tests run from a parallel no-AVX-512 build (build/vg).
VALGRIND ?= valgrind
VALGRIND_FLAGS = --error-exitcode=99 --leak-check=full \
                 --show-leak-kinds=definite,possible
VALGRIND_TESTS = test_add test_add_u64 test_cholesky test_cmp test_div \
                 test_divmod \
                 test_divmod_u64 test_field test_fincke_pohst test_gcd test_gcd_extended \
                 test_isqrt test_log_2 test_lshift test_mod \
                 test_mod_inverse test_mod_u64

ifneq ($(GMP_OK),yes)
VALGRIND_TESTS := $(filter-out $(GMP_TEST_NAMES), $(VALGRIND_TESTS))
endif

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
	$(CC) $(VG_CFLAGS) $^ -o $@ $(LDFLAGS) $(EXTRA_test_$*) \
	  $(if $(filter test_$*,$(GMP_TEST_NAMES)),-lgmp)

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

# --- Installation ---
PREFIX ?= /usr/local
INSTALL ?= install

rabin.pc: rabin.pc.in
	sed 's|@PREFIX@|$(PREFIX)|' $< > $@

install: $(LIB) rabin.pc
	$(INSTALL) -d $(DESTDIR)$(PREFIX)/lib $(DESTDIR)$(PREFIX)/lib/pkgconfig \
	    $(DESTDIR)$(PREFIX)/include/rabin
	$(INSTALL) -m 644 $(LIB) $(DESTDIR)$(PREFIX)/lib/librabin.a
	$(INSTALL) -m 644 include/*.h $(DESTDIR)$(PREFIX)/include/rabin/
	$(INSTALL) -m 644 rabin.pc $(DESTDIR)$(PREFIX)/lib/pkgconfig/

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/lib/librabin.a \
	    $(DESTDIR)$(PREFIX)/lib/pkgconfig/rabin.pc
	rm -rf $(DESTDIR)$(PREFIX)/include/rabin


# --- Profile-Guided Optimization (PGO) ---
.PHONY: pgo pgo-clean

pgo: clean
	@echo "==> Step 1/3: Compiling instrumented binaries..."
	$(MAKE) all CFLAGS="$(CFLAGS) -fprofile-generate" LDFLAGS="$(LDFLAGS) -fprofile-generate"
	@echo "==> Step 2/3: Generating profile data from workloads..."
	$(MAKE) test_complete CFLAGS="$(CFLAGS) -fprofile-generate" LDFLAGS="$(LDFLAGS) -fprofile-generate"
	@echo "==> Step 3/3: Rebuilding optimized library with profile feedback..."
	@find . -type f -name "*.o" -delete
	@find . -type f -name "*.a" -delete
	@rm -f $(BUILD_DIR)/bignum
	$(MAKE) all CFLAGS="$(CFLAGS) -fprofile-use -fprofile-correction -Wno-missing-profile" LDFLAGS="$(LDFLAGS) -fprofile-use"

pgo-clean: clean
	@find . -name "*.gcda" -o -name "*.gcno" -delete
