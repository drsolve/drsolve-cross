# Building DixonRes with CMake

CMake ≥ 3.16 is required.  
The build system produces a shared library (`libdixon.so` / `.dylib` / `-1.dll`),
a static library (`libdixon.a` / `-1.a`), and the `dixon` CLI executable.

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
| `DIXONRES_ENABLE_PML` | `ON` | Use PML if found |
| `DIXONRES_ENABLE_OPENMP` | `ON` | Use OpenMP if available |
| `DIXONRES_ENABLE_LTO` | `ON` | Link-time optimisation |
| `DIXONRES_ENABLE_ASAN` | `OFF` | AddressSanitizer |
| `DIXONRES_STATIC` | `OFF` | Link FLINT/PML statically into `dixon` |
| `DIXONRES_BUILD_SHARED_LIB` | `ON` | Build `libdixon.so` / `.dylib` |
| `DIXONRES_BUILD_STATIC_LIB` | `ON` | Build `libdixon.a` |
| `DIXONRES_BUILD_GUI` | `ON` | Build Windows GUI targets |
| `DIXONRES_BUILD_ATTACK` | `ON` | Build `../Attack/*.c` programs |

Examples:

```bash
# Non-standard FLINT location
cmake -B build -DFLINT_ROOT=$HOME/.local

# Fully static binary (no .so deps at runtime)
cmake -B build -DDIXONRES_STATIC=ON

# Debug build with AddressSanitizer
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DDIXONRES_ENABLE_ASAN=ON

# Shared library only (skip static library)
cmake -B build -DDIXONRES_BUILD_STATIC_LIB=OFF

# Install to custom prefix
cmake -B build -DCMAKE_INSTALL_PREFIX=/opt/dixonres
cmake --build build && cmake --install build
```

---

## Build outputs

### Linux / macOS

```
build/
  dixon                      ← CLI executable
  lib/
    libdixon.so.1.0.0        ← shared library (Linux)
    libdixon.so → .so.1      ← symlinks
    libdixon.dylib           ← shared library (macOS)
    libdixon.a               ← static library
```

### Windows

```
build-win/
  dixon.exe                  ← launcher (double-click or from Explorer)
  dixon_win_gui.exe          ← GUI frontend
  bin/dixon_cli_real.exe     ← actual CLI
  dll/libdixon-1.dll         ← Dixon shared library
  dll/*.dll                  ← runtime DLL dependencies
  lib/libdixon-1.a           ← static archive
  lib/libdixon-1.dll.a       ← import library
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

Then build (note: `-DCMAKE_TOOLCHAIN_FILE` must be an **absolute path**,
because CMake resolves relative paths against the build directory, not the source):

```bash
cmake -B build-win \
      -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake"
cmake --build build-win -j$(nproc)
```

The bundled `third_party/` directory is used automatically for FLINT and PML.

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