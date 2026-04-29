# Explicitly set default target - must be before all includes
.DEFAULT_GOAL := default

# ============================================================
# Include generated configuration (run ./configure if missing)
# ============================================================
ifeq ($(wildcard config.mk),)
$(info config.mk not found, running ./configure...)
_DUMMY := $(shell chmod +x configure && ./configure >&2)
ifeq ($(wildcard config.mk),)
$(error ./configure did not generate config.mk)
endif
endif
include config.mk

# ============================================================
# Compiler flags (configured for the current toolchain)
# ============================================================
ARCH_CFLAGS ?=
LTO_CFLAGS ?=
LTO_LDFLAGS ?=
OPENMP_AVAILABLE ?= no
OPENMP_CFLAGS ?=
OPENMP_LDFLAGS ?=

CFLAGS = -O3 $(ARCH_CFLAGS) -fPIC $(LTO_CFLAGS) $(OPENMP_CFLAGS) $(ASAN_CFLAGS)

# Linker flags
LDFLAGS = $(LTO_LDFLAGS) $(OPENMP_LDFLAGS) $(ASAN_LDFLAGS)

# ============================================================
# Local directories
# ============================================================
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

# ============================================================
# Bundled PML determinant subset (force use of pml_det over any system PML)
# ============================================================
PML_SOURCE_DIR ?= pml_det
PML_BUILD_DIR ?= $(BUILD_DIR)/pml_det
PML_OBJ_DIR ?= $(PML_BUILD_DIR)/obj
PML_BUNDLED_LIB ?= $(PML_BUILD_DIR)/libpml.a
PML_BUNDLED_AVAILABLE := $(if $(wildcard $(PML_SOURCE_DIR)/src/pml.h),yes,no)
PML_ORIGIN := configure

ifeq ($(PML_BUNDLED_AVAILABLE),yes)
PML_ORIGIN := bundled-pml_det
PML_INCLUDE_ROOT := $(PML_SOURCE_DIR)
PML_INCLUDE_PATH := $(PML_SOURCE_DIR)/src
PML_LIB_PATH := $(PML_BUILD_DIR)
PML_HEADER_CHECK := yes
PML_HEADER_PREFIX :=
PML_HEADER_DIR := $(PML_INCLUDE_PATH)
NMOD_POLY_MAT_UTILS_CHECK := yes
PML_DYNAMIC_LIB_CHECK := no
PML_STATIC_LIB_CHECK := yes
PML_SO_PATH :=
PML_A_PATH := $(PML_BUNDLED_LIB)
PML_DIR_EXISTS := yes
PML_AVAILABLE := yes
PML_FLAGS := -DHAVE_PML -DPML_HAVE_MACHINE_VECTORS=0
PML_LIBS := $(PML_BUNDLED_LIB)
PML_STATIC_LIBS := $(PML_BUNDLED_LIB)
INCLUDE_FLAGS := $(filter-out -I$(PML_INCLUDE_ROOT) -I$(PML_INCLUDE_PATH),$(INCLUDE_FLAGS))
INCLUDE_FLAGS := -I$(PML_INCLUDE_ROOT) -I$(PML_INCLUDE_PATH) $(INCLUDE_FLAGS)
RPATH_FLAGS := -Wl,-rpath,.
ifneq ($(strip $(FLINT_LIB_PATH)),)
RPATH_FLAGS += -Wl,-rpath,$(FLINT_LIB_PATH)
endif
PML_REL_SOURCES := \
	nmod_poly_mat_extra/nmod_poly_mat_det.c \
	nmod_poly_mat_extra/nmod_poly_mat_weak_popov.c \
	nmod_poly_mat_extra/nmod_poly_mat_utils.c \
	nmod_poly_mat_extra/degrees_pivots_leadingmatrix.c \
	nmod_poly_mat_extra/kernel.c \
	nmod_poly_mat_extra/approximant_basis.c \
	nmod_poly_mat_extra/middle_product.c \
	nmod_mat_poly_extra/nmod_mat_poly_mbasis.c \
	nmod_mat_poly_extra/nmod_mat_poly_mem.c \
	nmod_mat_poly_extra/nmod_mat_poly_arith.c \
	nmod_mat_poly_extra/nmod_mat_poly_set_from.c \
	nmod_mat_poly_extra/nmod_mat_poly_shift.c \
	nmod_mat_extra/left_nullspace.c \
	nmod_extra/nmod_find_root.c
PML_SOURCES := $(addprefix $(PML_SOURCE_DIR)/src/,$(PML_REL_SOURCES))
PML_OBJECTS := $(addprefix $(PML_OBJ_DIR)/,$(PML_REL_SOURCES:.c=.o))
PML_BUILD_PREREQS := $(PML_BUNDLED_LIB)
else
PML_BUILD_PREREQS :=
endif

# ============================================================
# Install directories
# Set by ./configure (--prefix, --bindir, etc.).
# You can still override at install time: make install PREFIX=/custom
# ============================================================

# Install tool
INSTALL         ?= install
INSTALL_PROGRAM ?= $(INSTALL) -m 755
INSTALL_DATA    ?= $(INSTALL) -m 644
INSTALL_DIR     ?= $(INSTALL) -d -m 755

# ============================================================
# Combined CFLAGS
# ============================================================
ALL_CFLAGS = $(CFLAGS) $(INCLUDE_FLAGS) $(FLINT_FLAGS) $(PML_FLAGS)

# ============================================================
# External library sets
# ============================================================
EXTERNAL_LIBS = $(PML_STATIC_LIBS) $(FLINT_LIBS) $(SYSTEM_LIBS)
EXTERNAL_STATIC_ALL_LIBS = $(PML_STATIC_LIBS) $(FLINT_STATIC_LIBS) \
                           -Wl,-Bdynamic $(SYSTEM_LIBS) -Wl,--allow-multiple-definition

# ============================================================
# Source files for the math library (in src directory)
# ============================================================
MATH_SOURCES = $(SRC_DIR)/dixon_complexity.c \
               $(SRC_DIR)/dixon_flint.c \
               $(SRC_DIR)/dixon_interface_flint.c \
               $(SRC_DIR)/dixon_test.c \
               $(SRC_DIR)/dixon_with_ideal_reduction.c \
               $(SRC_DIR)/fq_mat_det.c \
               $(SRC_DIR)/macaulay_flint.c \
               $(SRC_DIR)/fq_mpoly_mat_det.c \
               $(SRC_DIR)/fq_multivariate_interpolation.c \
               $(SRC_DIR)/fq_mvpoly.c \
               $(SRC_DIR)/fq_nmod_roots.c \
               $(SRC_DIR)/fq_poly_mat_det.c \
               $(SRC_DIR)/fq_sparse_interpolation.c \
               $(SRC_DIR)/fq_unified_interface.c \
               $(SRC_DIR)/gf2n_mpoly.c \
               $(SRC_DIR)/gf2n_field.c \
               $(SRC_DIR)/gf2n_poly.c \
               $(SRC_DIR)/large_prime_system_solver.c \
               $(SRC_DIR)/polynomial_system_solver.c \
               $(SRC_DIR)/resultant_with_ideal_reduction.c \
               $(SRC_DIR)/unified_mpoly_det.c \
               $(SRC_DIR)/unified_mpoly_interface.c \
               $(SRC_DIR)/unified_mpoly_resultant.c \
               $(SRC_DIR)/fmpq_acb_roots.c \
               $(SRC_DIR)/rational_system_solver.c
# Object files (in build directory)
MATH_OBJECTS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(MATH_SOURCES))

# Main source file (in current directory)
DIXON_SRC = drsolve.c

# All source files (for LTO compilation)
ALL_SOURCES = $(DIXON_SRC) $(MATH_SOURCES)

# Library names (in current directory)
DIXON_STATIC_LIB = libdrsolve.a
DIXON_SHARED_LIB = libdrsolve.$(SHARED_LIB_EXT)

# Output executable (in current directory)
DIXON_TARGET = drsolve

# ============================================================
# Attack programs directory and files
# ============================================================
ATTACK_DIR = ../Attack
# Find all C files recursively in Attack directory and subdirectories
ATTACK_C_FILES := $(shell find $(ATTACK_DIR) -name "*.c" 2>/dev/null | grep -v ".ipynb_checkpoints" || echo "")
ATTACK_EXECUTABLES := $(patsubst %.c,%,$(ATTACK_C_FILES))

# ============================================================
# Create build directory
# ============================================================
$(BUILD_DIR):
	@echo "Creating build directory..."
	mkdir -p $(BUILD_DIR)

ifeq ($(PML_BUNDLED_AVAILABLE),yes)
$(PML_OBJ_DIR)/%.o: $(PML_SOURCE_DIR)/src/%.c | $(BUILD_DIR)
	@echo "Compiling bundled PML source $<..."
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDE_FLAGS) $(FLINT_FLAGS) $(PML_FLAGS) -c -o $@ $<

$(PML_BUNDLED_LIB): $(PML_OBJECTS) | $(BUILD_DIR)
	@echo "Building bundled PML static library..."
	@mkdir -p $(dir $@)
	ar rcs $@ $^
	@echo "Bundled PML library built: $@"
endif

# ============================================================
# OS-specific flags for static linking with duplicate symbols
# (FLINT and PML share some symbol names intentionally)
# Linux GNU ld : --allow-multiple-definition, first definition wins
# macOS Apple ld: no equivalent flag; duplicates are resolved before
#   linking by making the FLINT copies local via nmedit (done in CI).
# ============================================================
ifeq ($(HOST_OS),Darwin)
  STATIC_ALLOW_MULTI =
  # Direct-path linking: avoids -Wl,-Bstatic which is GNU ld only.
  # PML listed before FLINT so PML's versions of shared symbols win.
  STATIC_LIBS_ORDERED = $(PML_STATIC_LIBS) $(FLINT_STATIC_LIBS)
else
  STATIC_ALLOW_MULTI = -Wl,--allow-multiple-definition -no-pie
  STATIC_LIBS_ORDERED = $(PML_STATIC_LIBS) $(FLINT_STATIC_LIBS)
endif

# ============================================================
# Default target
# Link mode is set by ./configure (--enable-static / --disable-static).
# You can also override directly: make static-all
# ============================================================
ifeq ($(DEFAULT_LINK_MODE),static-all)
# Static-all: only libdrsolve.a needed as prereq.
# Do NOT depend on $(DIXON_SHARED_LIB) — building a dylib against two
# static archives that share symbols fails on macOS (Apple ld has no
# --allow-multiple-definition).
default: $(DIXON_STATIC_LIB) $(PML_BUILD_PREREQS)
	@echo "Building $(DIXON_TARGET) fully static (configured by --enable-static)..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $(ALL_SOURCES) \
		$(DIXON_STATIC_LIB) \
		$(STATIC_LIBS_ORDERED) \
		$(SYSTEM_LIBS) \
		$(STATIC_ALLOW_MULTI) \
		$(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (fully static)"
else
default: $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB) $(PML_BUILD_PREREQS)
	@echo "Building $(DIXON_TARGET) with LTO (Link Time Optimization)..."
	@echo "Libraries built, now compiling all sources together for maximum inlining..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $(ALL_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (LTO optimized with libraries available)"
endif
	@echo ""
	@echo "=== Build Configuration ==="
	@echo "  Prefix  : $(PREFIX)"
	@echo "  Link mode: $(DEFAULT_LINK_MODE)"
ifeq ($(PML_AVAILABLE),yes)
	@echo "  PML support: ENABLED"
else
	@echo "  PML support: DISABLED"
endif
ifeq ($(ENABLE_ASAN),yes)
	@echo "  AddressSanitizer: ENABLED"
endif
	@echo "==========================="
ifneq ($(ATTACK_C_FILES),)
	@echo ""
	@echo "Now building Attack programs..."
	@$(MAKE) attack-programs-verbose
endif

# Also build libraries with LTO for better performance
all: default
	@echo "Built drsolve executable, libraries, and attack programs"
	@echo "  Link mode: $(DEFAULT_LINK_MODE)"
ifeq ($(PML_AVAILABLE),yes)
	@echo "  PML support: ENABLED"
else
	@echo "  PML support: DISABLED"
endif
ifeq ($(ENABLE_ASAN),yes)
	@echo "  AddressSanitizer: ENABLED"
endif

# LTO target - compile all sources together for maximum optimization (same as default now)
$(DIXON_TARGET)-lto: $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB) $(PML_BUILD_PREREQS)
	@echo "Building $(DIXON_TARGET) with LTO (Link Time Optimization)..."
	@echo "Libraries built, now compiling all sources together for maximum inlining..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $(ALL_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (LTO optimized)"

# Traditional dynamic library target
$(DIXON_TARGET)-dynamic: $(DIXON_SRC) $(DIXON_SHARED_LIB) $(PML_BUILD_PREREQS)
	@echo "Building $(DIXON_TARGET) with dynamic drsolve library..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< -L. -ldrsolve $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (dynamic drsolve, bundled static PML, dynamic FLINT)"

# ============================================================
# Library targets
# ============================================================

# Build dynamic drsolve library
dynamic-lib: $(DIXON_SHARED_LIB)

$(DIXON_SHARED_LIB): $(MATH_OBJECTS) $(PML_BUILD_PREREQS)
	@echo "Building dynamic drsolve library..."
	$(CC) $(SHARED_LDFLAGS) -o $@ $(MATH_OBJECTS) $(EXTERNAL_LIBS) $(LDFLAGS)
	@echo "Dynamic library built: $(DIXON_SHARED_LIB)"

# Build static drsolve library
static-lib: $(DIXON_STATIC_LIB)

$(DIXON_STATIC_LIB): $(MATH_OBJECTS)
	@echo "Building static drsolve library..."
	ar rcs $@ $^
	@echo "Static library built: $(DIXON_STATIC_LIB)"

# ============================================================
# Static linking variants
# ============================================================

# Build with static drsolve library (bundled static PML, dynamic FLINT)
static: $(DIXON_TARGET)-static
	@echo "Built drsolve with static library"

$(DIXON_TARGET)-static: $(DIXON_SRC) $(DIXON_STATIC_LIB) $(PML_BUILD_PREREQS)
	@echo "Building $(DIXON_TARGET) with static drsolve library (bundled static PML, dynamic FLINT)..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< $(DIXON_STATIC_LIB) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (static drsolve, bundled static PML, dynamic FLINT)"

# Build with all static libraries (drsolve + PML + FLINT)
static-all: static-lib $(DIXON_TARGET)-static-all

$(DIXON_TARGET)-static-all: $(DIXON_SRC) $(DIXON_STATIC_LIB) $(PML_BUILD_PREREQS)
	@echo "Building $(DIXON_TARGET) with all static libraries..."
	$(CC) $(ALL_CFLAGS) -o $(DIXON_TARGET) $< $(DIXON_STATIC_LIB) \
		$(STATIC_LIBS_ORDERED) \
		$(SYSTEM_LIBS) \
		$(STATIC_ALLOW_MULTI) \
		$(LDFLAGS)
	@echo "Build complete: $(DIXON_TARGET) (fully static)"

# ============================================================
# Install / Uninstall
# ============================================================

# install: copies binary, libraries, and headers to $(PREFIX)
# Usage:
#   make install              -> installs to /usr/local  (default)
#   make install PREFIX=~/.local
#   sudo make install PREFIX=/usr
install: $(DIXON_TARGET) $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	@echo "Installing to PREFIX=$(PREFIX) ..."
	@echo ""
	@echo "--- Creating directories ---"
	$(INSTALL_DIR) "$(DESTDIR)$(BINDIR)"
	$(INSTALL_DIR) "$(DESTDIR)$(LIBDIR)"
	$(INSTALL_DIR) "$(DESTDIR)$(INCLUDEDIR)"
	@echo ""
	@echo "--- Installing executable: $(DIXON_TARGET) -> $(DESTDIR)$(BINDIR)/ ---"
	$(INSTALL_PROGRAM) $(DIXON_TARGET) "$(DESTDIR)$(BINDIR)/$(DIXON_TARGET)"
	@echo ""
	@echo "--- Installing shared library: $(DIXON_SHARED_LIB) -> $(DESTDIR)$(LIBDIR)/ ---"
	$(INSTALL_DATA) $(DIXON_SHARED_LIB) "$(DESTDIR)$(LIBDIR)/$(DIXON_SHARED_LIB)"
	@# Create versioned symlink (libdrsolve.so.1 -> libdrsolve.so) if ldconfig available
	@if command -v ldconfig >/dev/null 2>&1; then \
		echo "Running ldconfig to update shared library cache..."; \
		ldconfig "$(DESTDIR)$(LIBDIR)" 2>/dev/null || true; \
	fi
	@echo ""
	@echo "--- Installing static library: $(DIXON_STATIC_LIB) -> $(DESTDIR)$(LIBDIR)/ ---"
	$(INSTALL_DATA) $(DIXON_STATIC_LIB) "$(DESTDIR)$(LIBDIR)/$(DIXON_STATIC_LIB)"
	@if command -v ranlib >/dev/null 2>&1; then \
		ranlib "$(DESTDIR)$(LIBDIR)/$(DIXON_STATIC_LIB)"; \
	fi
	@echo ""
	@echo "--- Installing headers: $(INCLUDE_DIR)/*.h -> $(DESTDIR)$(INCLUDEDIR)/ ---"
	@for h in $(INCLUDE_DIR)/*.h; do \
		if [ -f "$$h" ]; then \
			echo "  $$h -> $(DESTDIR)$(INCLUDEDIR)/"; \
			$(INSTALL_DATA) "$$h" "$(DESTDIR)$(INCLUDEDIR)/"; \
		fi; \
	done
	@echo ""
	@echo "=== Installation complete ==="
	@echo "  Binary  : $(DESTDIR)$(BINDIR)/$(DIXON_TARGET)"
	@echo "  Shared  : $(DESTDIR)$(LIBDIR)/$(DIXON_SHARED_LIB)"
	@echo "  Static  : $(DESTDIR)$(LIBDIR)/$(DIXON_STATIC_LIB)"
	@echo "  Headers : $(DESTDIR)$(INCLUDEDIR)/"

# install-strip: same as install but strips debug symbols from binary and shared lib
install-strip: $(DIXON_TARGET) $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB)
	@$(MAKE) install INSTALL_PROGRAM='$(INSTALL_PROGRAM) -s' \
	                 INSTALL_DATA='$(INSTALL_DATA)'
	@echo "Stripping installed shared library..."
	@strip --strip-unneeded "$(DESTDIR)$(LIBDIR)/$(DIXON_SHARED_LIB)" 2>/dev/null || true

# install-headers: install only the header files (useful for dev packages)
install-headers:
	@echo "Installing headers only..."
	$(INSTALL_DIR) "$(DESTDIR)$(INCLUDEDIR)"
	@for h in $(INCLUDE_DIR)/*.h; do \
		if [ -f "$$h" ]; then \
			echo "  $$h -> $(DESTDIR)$(INCLUDEDIR)/"; \
			$(INSTALL_DATA) "$$h" "$(DESTDIR)$(INCLUDEDIR)/"; \
		fi; \
	done
	@echo "Headers installed to $(DESTDIR)$(INCLUDEDIR)/"

# uninstall: remove everything that 'make install' put in place
uninstall:
	@echo "Uninstalling from PREFIX=$(PREFIX) ..."
	@echo ""
	@echo "--- Removing executable ---"
	rm -f "$(DESTDIR)$(BINDIR)/$(DIXON_TARGET)"
	@echo "--- Removing libraries ---"
	rm -f "$(DESTDIR)$(LIBDIR)/$(DIXON_SHARED_LIB)"
	rm -f "$(DESTDIR)$(LIBDIR)/$(DIXON_STATIC_LIB)"
	@if command -v ldconfig >/dev/null 2>&1; then \
		ldconfig "$(DESTDIR)$(LIBDIR)" 2>/dev/null || true; \
	fi
	@echo "--- Removing headers ---"
	rm -rf "$(DESTDIR)$(INCLUDEDIR)"
	@echo ""
	@echo "=== Uninstall complete ==="

# ============================================================
# Attack programs
# ============================================================

# Attack programs compilation with verbose output (LTO with all sources)
attack-programs-verbose: $(PML_BUILD_PREREQS)
	@echo "Building Attack programs with LTO optimization (compiling with all sources)..."
	@if [ -z "$(ATTACK_C_FILES)" ]; then \
		echo "No C files found in $(ATTACK_DIR)"; \
		echo "Checked directory: $(ATTACK_DIR)"; \
		if [ ! -d "$(ATTACK_DIR)" ]; then \
			echo "Directory $(ATTACK_DIR) does not exist!"; \
		fi; \
	else \
		echo "Found C files in $(ATTACK_DIR):"; \
		for cfile in $(ATTACK_C_FILES); do \
			echo "  $$cfile"; \
		done; \
		echo ""; \
		success_count=0; \
		fail_count=0; \
		for cfile in $(ATTACK_C_FILES); do \
			if [ -f "$$cfile" ]; then \
				executable="$${cfile%.c}"; \
				echo "Compiling $$cfile -> $$executable (LTO with all sources)"; \
				echo "Command: $(CC) $(ALL_CFLAGS) -o \"$$executable\" \"$$cfile\" $(MATH_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)"; \
				if $(CC) $(ALL_CFLAGS) -o "$$executable" "$$cfile" $(MATH_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS) 2>&1; then \
					echo "Successfully compiled: $$executable"; \
					ls -la "$$executable" | sed 's/^/  /'; \
					success_count=$$((success_count + 1)); \
				else \
					echo "Failed to compile $$cfile"; \
					fail_count=$$((fail_count + 1)); \
				fi; \
				echo ""; \
			else \
				echo "File not found: $$cfile"; \
				fail_count=$$((fail_count + 1)); \
			fi; \
		done; \
		echo "Attack programs compilation summary:"; \
		echo "  Success: $$success_count"; \
		echo "  Failed: $$fail_count"; \
		total_files=`echo "$(ATTACK_C_FILES)" | wc -w`; \
		echo "  Total: $$total_files"; \
	fi

# Attack programs compilation (silent version, LTO with all sources)
attack-programs: $(PML_BUILD_PREREQS)
	@echo "Building Attack programs with LTO optimization..."
	@if [ -z "$(ATTACK_C_FILES)" ]; then \
		echo "No C files found in $(ATTACK_DIR)"; \
	else \
		echo "Found C files: $(ATTACK_C_FILES)"; \
		for cfile in $(ATTACK_C_FILES); do \
			if [ -f "$$cfile" ]; then \
				executable="$${cfile%.c}"; \
				echo ""; \
				echo "Compiling $$cfile -> $$executable (LTO)"; \
				echo "Command: $(CC) $(ALL_CFLAGS) -o \"$$executable\" \"$$cfile\" $(MATH_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS)"; \
				if $(CC) $(ALL_CFLAGS) -o "$$executable" "$$cfile" $(MATH_SOURCES) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS); then \
					echo "Successfully compiled: $$executable"; \
					ls -la "$$executable" | sed 's/^/  /'; \
				else \
					echo "Failed to compile $$cfile"; \
				fi; \
			fi; \
		done; \
		echo ""; \
		echo "Attack programs compilation complete"; \
	fi

# Attack programs with static drsolve library (kept for compatibility)
attack-static: $(DIXON_STATIC_LIB) $(PML_BUILD_PREREQS)
	@echo "Building Attack programs with static drsolve library..."
	@if [ -z "$(ATTACK_C_FILES)" ]; then \
		echo "No C files found in $(ATTACK_DIR)"; \
	else \
		echo "Found C files: $(ATTACK_C_FILES)"; \
		for cfile in $(ATTACK_C_FILES); do \
			if [ -f "$$cfile" ]; then \
				executable="$${cfile%.c}"; \
				echo "Compiling $$cfile -> $$executable (static drsolve)"; \
				$(CC) $(ALL_CFLAGS) -o "$$executable" "$$cfile" $(DIXON_STATIC_LIB) $(EXTERNAL_LIBS) $(RPATH_FLAGS) $(LDFLAGS) || echo "Failed to compile $$cfile"; \
			fi; \
		done; \
		echo "Attack programs compilation complete (static drsolve)"; \
	fi

# Clean attack programs
clean-attack:
	@echo "Cleaning Attack programs..."
	@if [ -n "$(ATTACK_EXECUTABLES)" ]; then \
		for exe in $(ATTACK_EXECUTABLES); do \
			if [ -f "$$exe" ]; then \
				echo "Removing $$exe"; \
				rm -f "$$exe"; \
			fi; \
		done; \
	fi
	@echo "Attack programs cleaned"

# ============================================================
# Object file compilation (src/*.c -> build/*.o)
# ============================================================
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR) $(PML_BUILD_PREREQS)
	@echo "Compiling $<..."
	$(CC) $(ALL_CFLAGS) -c -o $@ $<

# ============================================================
# Clean
# ============================================================
clean: clean-attack
	rm -f $(DIXON_TARGET) $(DIXON_STATIC_LIB) $(DIXON_SHARED_LIB) dixon libdixon.a libdixon.so libdixon.dylib
	rm -rf $(BUILD_DIR)
	@echo "Cleaned all build artifacts"

# Clean only build directory (keep executables and libraries)
clean-build:
	rm -rf $(BUILD_DIR)
	@echo "Cleaned build directory"

# Clean everything including config.mk
distclean: clean
	rm -f config.mk
	@echo "Cleaned all build artifacts and configuration"

# ============================================================
# Test / debug / info targets
# ============================================================

# Test library detection
test-paths:
	@echo "Testing library path detection..."
	@echo "LD_LIBRARY_PATH: $$LD_LIBRARY_PATH"
	@echo "SYSTEM_LIB_PATHS: $(SYSTEM_LIB_PATHS)"
	@echo "PML_DYNAMIC_LIB_CHECK: $(PML_DYNAMIC_LIB_CHECK)"
	@echo "PML_STATIC_LIB_CHECK: $(PML_STATIC_LIB_CHECK)"
	@echo "PML_SO_PATH: $(PML_SO_PATH)"
	@echo "PML_A_PATH: $(PML_A_PATH)"
	@echo "PML_AVAILABLE: $(PML_AVAILABLE)"

# Test Attack directory detection
test-attack:
	@echo "Testing attack detection..."
	@echo "ATTACK_DIR: $(ATTACK_DIR)"
	@echo "ATTACK_C_FILES: $(ATTACK_C_FILES)"
	@echo "ATTACK_EXECUTABLES: $(ATTACK_EXECUTABLES)"
	@echo ""
	@echo "Manual find test:"
	@find $(ATTACK_DIR) -name "*.c" 2>/dev/null | grep -v ".ipynb_checkpoints" || echo "No files found or directory doesn't exist"

# Show configuration
info:
	@echo "=== Build Configuration ==="
	@echo "CC: $(CC)"
	@echo "CFLAGS: $(ALL_CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "RPATH_FLAGS: $(RPATH_FLAGS)"
ifeq ($(ENABLE_ASAN),yes)
	@echo "AddressSanitizer: ENABLED"
endif
	@echo ""
	@echo "=== Install Paths ==="
	@echo "PREFIX:     $(PREFIX)"
	@echo "BINDIR:     $(BINDIR)"
	@echo "LIBDIR:     $(LIBDIR)"
	@echo "INCLUDEDIR: $(INCLUDEDIR)"
	@echo ""
	@echo "=== Directory Structure ==="
	@echo "Source directory: $(SRC_DIR)/"
	@echo "Include directory: $(INCLUDE_DIR)/"
	@echo "Build directory: $(BUILD_DIR)/"
	@echo "Output directory: ./"
	@echo ""
	@echo "=== System Library Paths ==="
	@echo "LD_LIBRARY_PATH: $$LD_LIBRARY_PATH"
	@echo "Library search paths: $(SYSTEM_LIB_PATHS)"
	@echo ""
	@echo "=== Library Status ==="
	@echo "PML origin: $(PML_ORIGIN)"
	@echo "FLINT headers found: $(FLINT_HEADER_CHECK)"
	@echo "FLINT directory exists: $(FLINT_DIR_EXISTS) at $(FLINT_INCLUDE_PATH)"
	@echo "PML headers found: $(PML_HEADER_CHECK)"
	@echo "nmod_poly_mat_utils.h found: $(NMOD_POLY_MAT_UTILS_CHECK)"
	@echo "PML dynamic library found: $(PML_DYNAMIC_LIB_CHECK)"
	@echo "PML static library found: $(PML_STATIC_LIB_CHECK)"
	@echo "PML available (all headers + libraries): $(PML_AVAILABLE)"
	@echo "PML directory exists: $(PML_DIR_EXISTS) at $(PML_INCLUDE_PATH)"
	@echo ""
	@echo "=== Found Library Paths ==="
	@echo "PML dynamic library path: $(PML_SO_PATH)"
	@echo "PML static library path: $(PML_A_PATH)"
	@echo ""
	@echo "=== Library Paths ==="
	@echo "FLINT lib: $(FLINT_LIB_PATH)"
	@echo "PML lib: $(PML_LIB_PATH)"
	@echo ""
	@echo "=== Library Files ==="
	@echo "FLINT dynamic: $(FLINT_LIBS)"
	@echo "FLINT static: $(FLINT_STATIC_LIBS)"
	@echo "PML dynamic: $(PML_LIBS)"
	@echo "PML static: $(PML_STATIC_LIBS)"
	@echo ""
	@echo "=== Attack Directory Structure ==="
	@echo "Attack directory: $(ATTACK_DIR)"
	@echo -n "Directory exists: "
	@if [ -d "$(ATTACK_DIR)" ]; then \
		echo "YES"; \
		echo "C files found:"; \
		if [ -n "$(ATTACK_C_FILES)" ]; then \
			for cfile in $(ATTACK_C_FILES); do \
				echo "  $$cfile"; \
			done; \
		else \
			echo "  No C files found"; \
		fi; \
	else \
		echo "NO"; \
	fi
	@echo "EXTERNAL_LIBS (dynamic): $(EXTERNAL_LIBS)"
	@echo "EXTERNAL_STATIC_ALL_LIBS: $(EXTERNAL_STATIC_ALL_LIBS)"

# Debug header file detection
debug-headers:
	@echo "=== Header File Detection Debug ==="
	@echo ""
	@echo "PML origin: $(PML_ORIGIN)"
	@echo ""
	@echo "=== Compiler Search Paths ==="
	@echo "Getting GCC include search paths..."
	@$(CC) -E -v -x c /dev/null 2>&1 | sed -n '/#include <...> search starts here:/,/End of search list./p' | sed 's/^/ /'
	@echo ""
	@echo "=== Environment Variables ==="
	@echo "C_INCLUDE_PATH: $(C_INCLUDE_PATH)"
	@echo "CPLUS_INCLUDE_PATH: $(CPLUS_INCLUDE_PATH)"
	@echo ""
	@echo "=== Header File Tests ==="
	@echo -n "FLINT headers (flint/flint.h): $(FLINT_HEADER_CHECK)"
	@echo ""
	@echo -n "PML headers (pml.h): $(PML_HEADER_CHECK)"
	@echo ""
	@echo -n "nmod_poly_mat_utils.h: $(NMOD_POLY_MAT_UTILS_CHECK)"
	@echo ""
	@echo "PML Available (all required headers + libraries): $(PML_AVAILABLE)"
	@echo ""
	@echo "=== Manual Path Search ==="
	@echo "Searching for FLINT headers in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/flint/flint.h" ]; then \
			echo "  FOUND: $$path/flint/flint.h"; \
		fi; \
	done
	@echo "Searching for PML headers in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/pml.h" ]; then \
			echo "  FOUND: $$path/pml.h"; \
		fi; \
	done
	@echo "Searching for nmod_poly_mat_utils.h in common locations..."
	@for path in /usr/include /usr/local/include ~/.local/include $(subst :, ,$(C_INCLUDE_PATH)); do \
		if [ -f "$$path/nmod_poly_mat_utils.h" ]; then \
			echo "  FOUND: $$path/nmod_poly_mat_utils.h"; \
		fi; \
	done
# Debug library detection
debug-libs:
	@echo "=== Library Detection Debug ==="
	@echo ""
	@echo "=== System Library Paths ==="
	@echo "LD_LIBRARY_PATH: $$LD_LIBRARY_PATH"
	@echo "Detected paths: $(SYSTEM_LIB_PATHS)"
	@echo ""
	@echo "=== PML Library Search ==="
	@echo "PML origin: $(PML_ORIGIN)"
	@echo "Searching for PML libraries..."
	@echo -n "Dynamic libraries (libpml.so*): "
	@found=no; for path in $(SYSTEM_LIB_PATHS); do \
		if [ -n "$$path" ] && [ -d "$$path" ]; then \
			if ls "$$path"/libpml.so* >/dev/null 2>&1; then \
				echo "FOUND"; \
				ls "$$path"/libpml.so* 2>/dev/null | sed 's/^/  /'; \
				found=yes; break; \
			fi; \
		fi; \
	done; if [ "$$found" = "no" ]; then echo "NOT FOUND"; fi
	@echo -n "Static libraries (libpml.a): "
	@found=no; for path in $(SYSTEM_LIB_PATHS); do \
		if [ -n "$$path" ] && [ -f "$$path/libpml.a" ]; then \
			echo "FOUND"; \
			echo "  $$path/libpml.a"; \
			found=yes; break; \
		fi; \
	done; if [ "$$found" = "no" ]; then echo "NOT FOUND"; fi
	@echo ""
	@echo "=== Detection Results ==="
	@echo "PML headers found: $(PML_HEADER_CHECK)"
	@echo "nmod_poly_mat_utils.h found: $(NMOD_POLY_MAT_UTILS_CHECK)"
	@echo "PML dynamic library found: $(PML_DYNAMIC_LIB_CHECK)"
	@echo "PML static library found: $(PML_STATIC_LIB_CHECK)"
	@echo "PML available (all requirements met): $(PML_AVAILABLE)"
	@echo "Selected PML SO path: $(PML_SO_PATH)"
	@echo "Selected PML A path: $(PML_A_PATH)"

# Debug attack directory structure
debug-attack:
	@echo "=== Attack Directory Debug ==="
	@echo ""
	@echo "=== Attack Directory Analysis ==="
	@echo "Attack directory: $(ATTACK_DIR)"
	@echo -n "Directory exists: "
	@if [ -d "$(ATTACK_DIR)" ]; then \
		echo "YES"; \
		echo "Directory contents:"; \
		find $(ATTACK_DIR) -type f -name "*.c" 2>/dev/null | grep -v ".ipynb_checkpoints" | sed 's/^/  /' || echo "  No C files found"; \
		echo ""; \
		echo "Subdirectories:"; \
		find $(ATTACK_DIR) -type d 2>/dev/null | sed 's/^/  /' || echo "  No subdirectories"; \
		echo ""; \
		echo "All files (including checkpoints):"; \
		find $(ATTACK_DIR) -type f -name "*.c" 2>/dev/null | sed 's/^/  /' || echo "  No C files found"; \
	else \
		echo "NO"; \
	fi
	@echo ""
	@echo "=== Detected Variables ==="
	@echo "ATTACK_C_FILES: $(ATTACK_C_FILES)"
	@echo "ATTACK_EXECUTABLES: $(ATTACK_EXECUTABLES)"
	@echo ""
	@echo "=== Test Compilation Check ==="
	@echo "Current directory drsolve library status:"
	@if [ -f "$(DIXON_STATIC_LIB)" ]; then \
		echo "  $(DIXON_STATIC_LIB): EXISTS"; \
	else \
		echo "  $(DIXON_STATIC_LIB): MISSING (run 'make static-lib' first)"; \
	fi
	@if [ -f "$(DIXON_SHARED_LIB)" ]; then \
		echo "  $(DIXON_SHARED_LIB): EXISTS"; \
	else \
		echo "  $(DIXON_SHARED_LIB): MISSING (run 'make dynamic-lib' first)"; \
	fi
	@echo ""
	@echo "=== Compiled Attack Programs Status ==="
	@if [ -n "$(ATTACK_EXECUTABLES)" ]; then \
		for exe in $(ATTACK_EXECUTABLES); do \
			if [ -f "$$exe" ]; then \
				echo "  $$exe: EXISTS"; \
				ls -la "$$exe" | sed 's/^/    /'; \
			else \
				echo "  $$exe: NOT COMPILED"; \
			fi; \
		done; \
	else \
		echo "  No attack executables expected"; \
	fi

debug-structure:
	@echo "=== Local Directory Structure Debug ==="
	@echo ""
	@echo "=== Current Directory ==="
	@echo "PWD: $(shell pwd)"
	@echo "Contents:"
	@ls -la . | sed 's/^/  /'
	@echo ""
	@echo "=== Source Directory ($(SRC_DIR)) ==="
	@echo -n "Directory exists: "
	@if [ -d "$(SRC_DIR)" ]; then \
		echo "YES"; \
		echo "Contents:"; \
		ls -la $(SRC_DIR) | sed 's/^/  /'; \
	else \
		echo "NO"; \
	fi
	@echo ""
	@echo "=== Include Directory ($(INCLUDE_DIR)) ==="
	@echo -n "Directory exists: "
	@if [ -d "$(INCLUDE_DIR)" ]; then \
		echo "YES"; \
		echo "Contents:"; \
		ls -la $(INCLUDE_DIR) | sed 's/^/  /'; \
	else \
		echo "NO"; \
	fi
	@echo ""
	@echo "=== Build Directory ($(BUILD_DIR)) ==="
	@echo -n "Directory exists: "
	@if [ -d "$(BUILD_DIR)" ]; then \
		echo "YES"; \
		echo "Contents:"; \
		ls -la $(BUILD_DIR) | sed 's/^/  /'; \
	else \
		echo "NO (will be created during build)"; \
	fi

# ============================================================
# Help
# ============================================================
help:
	@echo "Available targets:"
	@echo "  make (default)       - Build libraries first, then drsolve with LTO (all sources compiled together)"
	@echo "  make all             - Same as default"
	@echo "  make lto             - Same as default - Build with Link Time Optimization"
	@echo "  make dynamic         - Build drsolve with dynamic drsolve library"
	@echo "  make static          - Build drsolve with static drsolve library (bundled static PML, dynamic FLINT)"
	@echo "  make static-all      - Build drsolve with all static libraries (fully static)"
	@echo "  make dynamic-lib     - Build dynamic drsolve library only"
	@echo "  make static-lib      - Build static drsolve library only"
	@echo "  make attack-programs - Build all C programs in ../Attack directory (LTO with all sources)"
	@echo "  make attack-programs-verbose - Build Attack programs with detailed output (LTO)"
	@echo "  make attack-static   - Build all C programs in ../Attack directory (using static drsolve library)"
	@echo "  make clean-attack    - Clean all compiled Attack programs"
	@echo "  make test-paths      - Test library path detection"
	@echo "  make test-attack     - Test Attack directory detection"
	@echo "  make info            - Show build configuration (including install paths)"
	@echo "  make debug-headers   - Debug header file detection"
	@echo "  make debug-libs      - Debug external library detection"
	@echo "  make debug-structure - Debug local directory structure"
	@echo "  make debug-attack    - Debug Attack directory structure and C files"
	@echo "  make clean           - Clean all build artifacts (including Attack programs)"
	@echo "  make check           - Run test suite (pass/fail summary)"
	@echo "  make check-verbose   - Run tests with full output"
	@echo "  make clean-build     - Clean only build directory"
	@echo "  make distclean       - Clean all build artifacts and config.mk"
	@echo ""
	@echo "Install targets:"
	@echo "  make install         - Install binary, libraries, and headers to PREFIX (default: /usr/local)"
	@echo "  make install-strip   - Same as install but strips debug symbols"
	@echo "  make install-headers - Install header files only"
	@echo "  make uninstall       - Remove all installed files"
	@echo ""
	@echo "Install path overrides (all optional):"
	@echo "  make install PREFIX=~/.local"
	@echo "  make install PREFIX=/usr LIBDIR=/usr/lib/x86_64-linux-gnu"
	@echo "  make install DESTDIR=/tmp/staging PREFIX=/usr  # for package staging"
	@echo ""
	@echo "Build workflow:"
	@echo "  1. ./configure       - Detect libraries, generate config.mk"
	@echo "  2. make              - Build everything"
	@echo "  3. sudo make install - Install to /usr/local"
	@echo ""
	@echo "Directory structure:"
	@echo "  $(SRC_DIR)/          - Source files (.c)"
	@echo "  $(INCLUDE_DIR)/      - Header files (.h)"
	@echo "  $(BUILD_DIR)/        - Object files (.o) [created during build]"
	@echo "  $(ATTACK_DIR)/       - Attack programs (.c files with main())"
	@echo "  ./               - Executables and libraries"
	@echo ""
	@echo "Attack programs workflow:"
	@echo "  1. Run 'make' to build drsolve libraries AND all Attack programs automatically"
	@echo "  2. Or run 'make attack-programs' to build only Attack programs with LTO optimization"
	@echo "  3. Or run 'make attack-programs-verbose' for detailed compilation output"
	@echo "  4. Or run 'make attack-static' to build with static drsolve library (legacy)"
	@echo "  5. Each .c file in ../Attack becomes an executable with the same name"
	@echo "  6. Use 'make debug-attack' to check compilation status"
	@echo "  7. .ipynb_checkpoints directories are automatically excluded"
	@echo ""
	@echo "Attack programs compilation strategy:"
	@echo "  default - Compile Attack programs with all drsolve sources using LTO (best performance)"
	@echo "  attack-static - Use pre-built static drsolve library (legacy compatibility)"
	@echo ""
	@echo "Compilation strategy:"
	@echo "  default - Build bundled PML first, then drsolve libraries, then compile all sources with LTO"
	@echo "  dynamic - Traditional library-based compilation using pre-built drsolve library"
	@echo "  static  - Static drsolve + bundled static PML + dynamic FLINT (needs rpath)"
	@echo "  static-all  - Fully static (no runtime dependencies)"
	@echo ""
	@echo "Library structure:"
	@echo "  drsolve library: $(words $(MATH_SOURCES)) math source files"
	@echo "  Main program: drsolve.c links against drsolve library OR compiles with all sources"
	@echo "  Attack programs: Each .c in ../Attack compiles with all drsolve sources (LTO optimization)"
	@echo "  External deps: FLINT (required), PML determinant subset (bundled from pml_det when present)"
	@echo ""
	@echo "PML Detection:"
	@echo "  PML is forced to use the bundled pml_det copy when it exists"
	@echo "  The bundled library is built as $(PML_BUNDLED_LIB)"
	@echo "  Re-run ./configure for FLINT detection, then 'make test-paths' to verify"

# ============================================================
# Aliases for convenience
# ============================================================
lto: $(DIXON_TARGET)-lto
dynamic: $(DIXON_TARGET)-dynamic

# ============================================================
# Check / test targets
# ============================================================

# Colour helpers (fall back gracefully if terminal doesn't support it)
_GREEN  := \033[0;32m
_RED    := \033[0;31m
_YELLOW := \033[0;33m
_NC     := \033[0m   # No Colour

# Run a single test.
# Usage: $(call RUN_TEST, description, command)
# Increments PASS/FAIL counters; prints coloured result.
define RUN_TEST
	@printf "  %-60s" "$(1)"; \
	if $(2) >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; \
		PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; \
		FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    $(1)"; \
	fi
endef

check: $(DIXON_TARGET)
	@echo ""
	@echo "╔══════════════════════════════════════════════════════════════╗"
	@echo "║                    Dixon Test Suite                          ║"
	@echo "╚══════════════════════════════════════════════════════════════╝"
	@echo ""
	@PASS=0; FAIL=0; FAILED_TESTS=""; \
	\
	echo "--- Basic Dixon resultant ---"; \
	\
	printf "  %-60s" "Dixon: x+y+z, x*y+y*z+z*x, x*y*z+1 over F_257"; \
	if ./$(DIXON_TARGET) "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Dixon: x+y+z, x*y+y*z+z*x, x*y*z+1 over F_257"; \
	fi; \
	\
	printf "  %-60s" "Dixon: x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1 over F_257"; \
	if ./$(DIXON_TARGET) "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" "x,y" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Dixon: x^2+y^2+z^2-6 over F_257"; \
	fi; \
	\
	printf "  %-60s" "Dixon: extension field 2^8 (silent)"; \
	if ./$(DIXON_TARGET) --silent "x+y^2+t, x*y+t*y+1" "x" "2^8" >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Dixon: extension field 2^8 (silent)"; \
	fi; \
	\
	echo ""; \
	echo "--- Complexity analysis (--comp / -c) ---"; \
	\
	printf "  %-60s" "Comp: x+y+z, x*y+y*z+z*x, x*y*z+1 over F_257"; \
	if ./$(DIXON_TARGET) --comp "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Comp: x+y+z over F_257"; \
	fi; \
	\
	printf "  %-60s" "Comp -c: x^2+y^2+1, x*y+z, x+y+z^2 over F_257"; \
	if ./$(DIXON_TARGET) -c "x^2+y^2+1, x*y+z, x+y+z^2" "x,y" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Comp -c: x^2+y^2+1 over F_257"; \
	fi; \
	\
	printf "  %-60s" "Comp --omega: [4]*4 over F_65537"; \
	if ./$(DIXON_TARGET) --comp --omega 2.81 "x^4+y^4+z^4+w^4+1, x^3*y+z+1, x+y^3+z^2+w, x*y*z*w+1" "x,y,z" 65537 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Comp --omega [4]*4"; \
	fi; \
	\
	echo ""; \
	echo "--- Polynomial solver (--solve) ---"; \
	\
	printf "  %-60s" "Solve: x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1 over F_257"; \
	if ./$(DIXON_TARGET) --solve "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Solve: x^2+y^2+z^2-6 over F_257"; \
	fi; \
	\
	printf "  %-60s" "Solve: simple linear x+y-3, x-y+1 over F_257"; \
	if ./$(DIXON_TARGET) --solve "x+y-3, x-y+1" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Solve: x+y-3 linear over F_257"; \
	fi; \
	\
	echo ""; \
	echo "--- Random mode (--random / -r) ---"; \
	\
	printf "  %-60s" "Random Dixon: [3,3,2] over F_257"; \
	if ./$(DIXON_TARGET) --random "[3,3,2]" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Random Dixon: [3,3,2] over F_257"; \
	fi; \
	\
	printf "  %-60s" "Random solve -r: [2]*3 over F_257"; \
	if ./$(DIXON_TARGET) -r --solve "[2]*3" 257 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Random solve: [2]*3 over F_257"; \
	fi; \
	\
	printf "  %-60s" "Random comp -r: [3,2]*2 over F_65537"; \
	if ./$(DIXON_TARGET) -r --comp "[3,2]*2" 65537 >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Random comp: [3,2]*2 over F_65537"; \
	fi; \
	\
	echo ""; \
	echo "--- File input ---"; \
	\
	printf "  %-60s" "File: basic Dixon from generated file"; \
	printf "x,y\n257\nx+y+z, x*y+y*z+z*x, x*y*z+1\n" > /tmp/dixon_check_test.dat; \
	if ./$(DIXON_TARGET) /tmp/dixon_check_test.dat >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    File: basic Dixon from file"; \
	fi; \
	rm -f /tmp/dixon_check_test.dat; \
	\
	printf "  %-60s" "File: solver from generated file"; \
	printf "257\nx^2+y^2-5\nx+y-3\n" > /tmp/dixon_check_solver.dat; \
	if ./$(DIXON_TARGET) /tmp/dixon_check_solver.dat >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    File: solver from file"; \
	fi; \
	rm -f /tmp/dixon_check_solver.dat; \
	\
	echo ""; \
	echo "--- Error / edge cases ---"; \
	\
	printf "  %-60s" "Help flag exits cleanly"; \
	if ./$(DIXON_TARGET) --help >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_RED)[FAIL]$(_NC)\n"; FAIL=$$((FAIL+1)); \
		FAILED_TESTS="$$FAILED_TESTS\n    Help flag"; \
	fi; \
	\
	printf "  %-60s" "No args prints usage (exit 0)"; \
	if ./$(DIXON_TARGET) >/dev/null 2>&1; then \
		printf "$(_GREEN)[PASS]$(_NC)\n"; PASS=$$((PASS+1)); \
	else \
		printf "$(_YELLOW)[SKIP/WARN]$(_NC) (non-zero exit with no args)\n"; \
	fi; \
	\
	echo ""; \
	echo "════════════════════════════════════════════════════════════════"; \
	TOTAL=$$((PASS+FAIL)); \
	if [ $$FAIL -eq 0 ]; then \
		printf "Result: $(_GREEN)All $$TOTAL tests passed.$(_NC)\n"; \
	else \
		printf "Result: $(_RED)$$FAIL of $$TOTAL tests FAILED.$(_NC)\n"; \
		printf "Failed tests:$$FAILED_TESTS\n"; \
	fi; \
	echo "════════════════════════════════════════════════════════════════"; \
	echo ""; \
	rm -f solution_*.dat comp_*.dat; \
	exit $$FAIL

# Verbose check: same tests but shows full output of each command
check-verbose: $(DIXON_TARGET)
	@echo ""
	@echo "=== Verbose Test Run ==="
	@echo ""
	@set -e; \
	echo "--- Test 1: Basic Dixon ---"; \
	./$(DIXON_TARGET) "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257; \
	echo ""; \
	echo "--- Test 2: Complexity analysis ---"; \
	./$(DIXON_TARGET) --comp "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257; \
	echo ""; \
	echo "--- Test 3: Solver ---"; \
	./$(DIXON_TARGET) --solve "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" 257; \
	echo ""; \
	echo "--- Test 4: Random Dixon [3,3,2] ---"; \
	./$(DIXON_TARGET) --random "[3,3,2]" 257; \
	echo ""; \
	echo "--- Test 5: Random solver [2]*3 ---"; \
	./$(DIXON_TARGET) -r --solve "[2]*3" 257; \
	echo ""; \
	echo "--- Test 6: Extension field 2^8 (silent) ---"; \
	./$(DIXON_TARGET) --silent "x+y^2+t, x*y+t*y+1" "y" "2^8"; \
	echo ""; \
	echo "=== All verbose tests passed ==="; \
	rm -f solution_*.dat comp_*.dat

.PHONY: default all lto dynamic static static-all dynamic-lib static-lib \
        attack-programs attack-programs-verbose attack-static clean-attack \
        clean clean-build distclean test-paths test-attack info \
        debug-headers debug-libs debug-structure debug-attack help \
        install install-strip install-headers uninstall \
        check check-verbose
