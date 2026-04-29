# Building DRsolve with CMake

CMake ≥ 3.16 is required.  
The build system produces a shared library (`libdrsolve.so` / `.dylib` / `-1.dll`),
a static library (`libdrsolve.a` / `-1.a`), and the `drsolve` CLI executable.

---

## Linux / macOS  (quick start)

```bash
cmake -B build
cmake --build build -j$(nproc)
ctest --test-dir build               # run tests
sudo cmake --install build           # installs to /usr/local
```

---

## Common CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `CMAKE_BUILD_TYPE` | `Release` | `Release` / `Debug` / `RelWithDebInfo` |
| `CMAKE_INSTALL_PREFIX` | `/usr/local` | Install prefix |
| `FLINT_ROOT` | *(auto)* | Root of FLINT installation (contains `include/` and `lib/`) |
| `PML_ROOT` | *(auto)* | Root of PML installation |
| `DRSOLVE_ENABLE_PML` | `ON` | Use PML if found |
| `DRSOLVE_ENABLE_OPENMP` | `ON` | Use OpenMP if available |
| `DRSOLVE_ENABLE_LTO` | `ON` | Link-time optimisation |
| `DRSOLVE_ENABLE_ASAN` | `OFF` | AddressSanitizer |
| `DRSOLVE_STATIC` | `OFF` | Link FLINT/PML statically into `drsolve` |
| `DRSOLVE_BUILD_SHARED_LIB` | `ON` | Build `libdrsolve.so` / `.dylib` |
| `DRSOLVE_BUILD_STATIC_LIB` | `ON` | Build `libdrsolve.a` |
| `DRSOLVE_BUILD_GUI` | `ON` | Build Windows GUI targets |
| `DRSOLVE_BUILD_ATTACK` | `ON` | Build `../Attack/*.c` programs |
| `DRSOLVE_USE_BUNDLED_DEPS` | `ON` | Use bundled third-party dependencies (cross-compile only) |

Examples:

```bash
# Non-standard FLINT location
cmake -B build -DFLINT_ROOT=$HOME/.local

# Fully static binary (no .so deps at runtime)
cmake -B build -DDRSOLVE_STATIC=ON

# Debug build with AddressSanitizer
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DDRSOLVE_ENABLE_ASAN=ON

# Shared library only (skip static library)
cmake -B build -DDRSOLVE_BUILD_STATIC_LIB=OFF

# Install to custom prefix
cmake -B build -DCMAKE_INSTALL_PREFIX=/opt/drsolve
cmake --build build && cmake --install build
```

---

## Build outputs

### Linux / macOS

```
build/
  drsolve                    ← CLI executable
  lib/
    libdrsolve.so.1.0.0        ← shared library (Linux)
    libdrsolve.so → .so.1    ← symlinks
    libdrsolve.dylib           ← shared library (macOS)
    libdrsolve.a               ← static library
```

### Windows

```
build-win/
  drsolve.exe                  ← launcher (double-click or from Explorer)
  drsolve_win_gui.exe          ← GUI frontend
  bin/drsolve_cli_real.exe     ← actual CLI
  dll/libdrsolve-1.dll         ← Dixon shared library
  dll/*.dll                  ← runtime DLL dependencies
  lib/libdrsolve-1.a           ← static archive
  lib/libdrsolve-1.dll.a       ← import library
```

---

## Windows — cross-compile from Linux/macOS (MinGW-w64)

Install the cross-compiler first:

```bash
# Ubuntu/Debian
sudo apt install gcc-mingw-w64-x86-64

# Homebrew (macOS)
brew install mingw-w64
```

### Option A: Use bundled third-party dependencies (default)

The bundled `mingw/`, `pml_det/`, and `runtime/` directories are used automatically:

```bash
cmake -B build-win \
      -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake"
cmake --build build-win -j$(nproc)
```

### Option B: Auto-download from MSYS2 (recommended for clean Git repos)

CMake will automatically download required DLLs and libraries directly from MSYS2
during the build:

```bash
cmake -B build-win \
      -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake" \
      -DDRSOLVE_USE_BUNDLED_DEPS=OFF
cmake --build build-win -j$(nproc)
```

**How it works:**
- Packages are downloaded from https://packages.msys2.org/
- DLLs are automatically extracted and copied to `build-win/dll/`
- No need to commit large binary files in Git

**Updating package versions/SHA256:**

1. Visit https://packages.msys2.org/ and find the latest versions:
   - FLINT: https://packages.msys2.org/package/mingw-w64-ucrt-x86_64-flint
   - GMP: https://packages.msys2.org/package/mingw-w64-ucrt-x86_64-gmp
   - MPFR: https://packages.msys2.org/package/mingw-w64-ucrt-x86_64-mpfr
   - OpenBLAS: https://packages.msys2.org/package/mingw-w64-ucrt-x86_64-openblas
   - GCC Libs: https://packages.msys2.org/package/mingw-w64-ucrt-x86_64-gcc-libs
   - WinPthread: https://packages.msys2.org/package/mingw-w64-ucrt-x86_64-libwinpthread

2. Copy the "Version" and "SHA256" values from each page and update them in `CMakeLists.txt` (lines 94-138)

**Note about PML:**
PML library is not available in MSYS2 repositories. You still need to provide `PML_ROOT` explicitly or use the bundled determinant subset in `pml_det/`.

### Option C: Use system cross-compiler libraries

If you have MinGW-w64 libraries installed system-wide via your package manager:

```bash
cmake -B build-win \
      -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake" \
      -DDRSOLVE_USE_BUNDLED_DEPS=OFF
cmake --build build-win -j$(nproc)
```

Note: Ensure FLINT and other dependencies are installed for the MinGW-w64 cross-compiler.

---

## Windows — native build (MSYS2/UCRT64 or Visual Studio)

For MSYS2/UCRT64:

```bash
pacman -S mingw-w64-ucrt-x86_64-cmake \
          mingw-w64-ucrt-x86_64-gcc \
          mingw-w64-ucrt-x86_64-flint

cmake -B build -G "MinGW Makefiles"
cmake --build build -j$(nproc)
```

For Visual Studio: FLINT does not currently ship MSVC-compatible libs;
use the MinGW approach above.

---

## Running tests

```bash
ctest --test-dir build --output-on-failure
```

---

## Uninstall

CMake does not provide a built-in uninstall target.
Use the generated install manifest:

```bash
xargs rm -f < build/install_manifest.txt
```
