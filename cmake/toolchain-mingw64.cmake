# cmake/toolchain-mingw64.cmake
#
# CMake toolchain for cross-compiling DRsolve on Linux/macOS
# targeting 64-bit Windows using the MinGW-w64 toolchain.
#
# Usage:
#   cmake -B build-win \
#         -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-mingw64.cmake \
#         [-DFLINT_ROOT=/path/to/flint-win] \
#         [-DPML_ROOT=/path/to/pml-win]
#
# The cross-compiler (x86_64-w64-mingw32-gcc) must be installed,
# e.g. via:
#   Ubuntu:  sudo apt install gcc-mingw-w64-x86-64
#   Fedora:  sudo dnf install mingw64-gcc
#   Homebrew: brew install mingw-w64

set(CMAKE_SYSTEM_NAME    Windows)
set(CMAKE_SYSTEM_VERSION 6.1)        # Windows 7+
set(CMAKE_SYSTEM_PROCESSOR x86_64)

# Cross-compiler executables
set(CMAKE_C_COMPILER   x86_64-w64-mingw32-gcc)
set(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++)
set(CMAKE_RC_COMPILER  x86_64-w64-mingw32-windres)
set(CMAKE_AR           x86_64-w64-mingw32-ar)
set(CMAKE_RANLIB       x86_64-w64-mingw32-ranlib)
set(CMAKE_STRIP        x86_64-w64-mingw32-strip)

# Tell CMake where to find target-platform libraries.
# Ubuntu/Debian installs MinGW sysroot under /usr/x86_64-w64-mingw32;
# Homebrew on macOS uses a different layout under its prefix.
if(EXISTS /usr/x86_64-w64-mingw32)
  # Linux (Ubuntu/Debian/Fedora)
  set(CMAKE_FIND_ROOT_PATH /usr/x86_64-w64-mingw32)
elseif(EXISTS /opt/homebrew/opt/mingw-w64/toolchain-x86_64)
  # Homebrew arm64 Mac
  set(CMAKE_FIND_ROOT_PATH /opt/homebrew/opt/mingw-w64/toolchain-x86_64)
elseif(EXISTS /usr/local/opt/mingw-w64/toolchain-x86_64)
  # Homebrew Intel Mac
  set(CMAKE_FIND_ROOT_PATH /usr/local/opt/mingw-w64/toolchain-x86_64)
endif()

# Search headers/libs in the cross-root; build tools from host
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH) 

# The project bundles Windows third-party deps in third_party/
# These paths are set automatically by CMakeLists.txt when
# DRSOLVE_WINDOWS is true and the user has not set explicit roots.
# Override here if needed:
#
#   set(FLINT_ROOT "${CMAKE_SOURCE_DIR}/third_party/mingw")
#   set(PML_ROOT   "${CMAKE_SOURCE_DIR}/third_party/pml")