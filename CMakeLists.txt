cmake_minimum_required(VERSION 3.16)

project(drsolve
  VERSION 0.2.4
  DESCRIPTION "DRsolve: Dixon Resultant & Polynomial System Solver"
  LANGUAGES C
)

# ============================================================
# Build options
# ============================================================
option(DRSOLVE_ENABLE_PML       "Use bundled pml_det when available" ON)
option(DRSOLVE_ENABLE_OPENMP    "Enable OpenMP if available"              ON)
option(DRSOLVE_ENABLE_LTO       "Enable Link Time Optimization (-flto)"   ON)
option(DRSOLVE_ENABLE_ASAN      "Enable AddressSanitizer"                 OFF)
option(DRSOLVE_STATIC           "Link FLINT/PML statically into binary"   OFF)
option(DRSOLVE_BUILD_SHARED_LIB "Build libdrsolve as a shared library"      ON)
option(DRSOLVE_BUILD_STATIC_LIB "Build libdrsolve as a static library"      ON)
option(DRSOLVE_BUILD_GUI        "Build Windows GUI targets (Windows only)" ON)
option(DRSOLVE_BUILD_ATTACK     "Build programs in ../Attack if present"  ON)
option(DRSOLVE_USE_BUNDLED_DEPS "Use bundled third-party dependencies (cross-compile only)" ON)

# Hint paths exposed to the user  (override via -DFLINT_ROOT=... etc.)
set(FLINT_ROOT ""   CACHE PATH "Root directory of FLINT installation")
set(PML_ROOT   ""   CACHE PATH "Root directory of PML installation")

# ============================================================
# Default build type
# ============================================================
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    Release Debug RelWithDebInfo MinSizeRel)
endif()

# ============================================================
# Standard modules
# ============================================================
include(GNUInstallDirs)
include(CheckIPOSupported)
include(CheckCCompilerFlag)
include(CTest)

# Include our custom MSYS2 package fetcher
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ============================================================
# Platform detection
# ============================================================
if(WIN32 OR MINGW OR CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(DRSOLVE_WINDOWS TRUE)
else()
  set(DRSOLVE_WINDOWS FALSE)
endif()

# Detect if we're cross-compiling vs native Windows build
if(DRSOLVE_WINDOWS AND CMAKE_CROSSCOMPILING)
  set(DRSOLVE_WINDOWS_CROSS TRUE)
  message(STATUS "Windows: cross-compiling from ${CMAKE_HOST_SYSTEM_NAME}")
elseif(DRSOLVE_WINDOWS)
  set(DRSOLVE_WINDOWS_NATIVE TRUE)
  message(STATUS "Windows: native build (MSYS2/MINGW64)")
else()
  set(DRSOLVE_WINDOWS_CROSS FALSE)
  set(DRSOLVE_WINDOWS_NATIVE FALSE)
endif()

# When cross-compiling for Windows from Linux/macOS:
#   - Optionally use bundled mingw/ and pml_det/ directories, OR
#   - Use system-installed cross-compiler libraries
# Native Windows builds (MSYS2) should use pacman-installed packages.
if(DRSOLVE_WINDOWS_CROSS AND DRSOLVE_USE_BUNDLED_DEPS AND NOT FLINT_ROOT)
  if(EXISTS "${CMAKE_SOURCE_DIR}/mingw/include/flint/flint.h")
    set(FLINT_ROOT "${CMAKE_SOURCE_DIR}/mingw")
  elseif(EXISTS "${CMAKE_SOURCE_DIR}/flint/include/flint/flint.h")
    set(FLINT_ROOT "${CMAKE_SOURCE_DIR}/flint")
  endif()
endif()


# ============================================================
# Fetch MSYS2 packages for cross-compilation (when not using bundled deps)
# ============================================================
if(DRSOLVE_WINDOWS_CROSS AND NOT DRSOLVE_USE_BUNDLED_DEPS)
  include(FetchMSYS2Package)
  
  message(STATUS "Fetching MSYS2 packages for cross-compilation...")
  
  # MSYS2 package definitions - update these as needed
  # You can find latest versions and SHA256 at https://packages.msys2.org/
  set(MSYS2_PACKAGES)
  
  # FLINT (Fast Library for Number Theory)
  list(APPEND MSYS2_PACKAGES
    "mingw-w64-x86_64-flint"
    "3.5.0-1"
    "9194082fcafbdf6e7a3f6c722f428f7e9948927d3e2b440525e0582a42ab0ad9"
  )
  
  # GMP (GNU Multiple Precision Arithmetic Library)
  list(APPEND MSYS2_PACKAGES
    "mingw-w64-x86_64-gmp"
    "6.3.0-2"
    "8924433974c4add46cb46ea4f6ef283b5c5139d3f552375115b5580f855015cc"
  )
  
  # MPFR (GNU Multiple Precision Floating-Point Library)
  list(APPEND MSYS2_PACKAGES
    "mingw-w64-x86_64-mpfr"
    "4.2.2-1"
    "d774caec6c1655970fb9908b1b9521ce398ea394da62aab48fb2f050845dcd3c"
  )
  
  # OpenBLAS (Optimized BLAS library)
  list(APPEND MSYS2_PACKAGES
    "mingw-w64-x86_64-openblas"
    "0.3.31-1"
    "744ec919b71ce7e148678bc67426af03a78f20a383390330a6d94a4fc4f7675d"
  )
  
  # GCC runtime libraries
  list(APPEND MSYS2_PACKAGES
    "mingw-w64-x86_64-gcc-libs"
    "15.2.0-13"
    "5ccdbb9cbc1356f92cb77c22429c67b8d7ab0976e84293c8297548284826fb72"
  )
  
  # WinPthread
  list(APPEND MSYS2_PACKAGES
    "mingw-w64-x86_64-libwinpthread"
    "13.0.0.r560.g3197fc7d6-1"
    "b12bd011c79bcac1c742a94a8222d782913dcdebd0298aebb08ae90556405c0f"
  )
  
  # Download and extract all packages
  set(ALL_FETCHED_PACKAGES "")
  list(LENGTH MSYS2_PACKAGES NUM_PKG_ENTRIES)
  math(EXPR NUM_PKGS "${NUM_PKG_ENTRIES} / 3")
  
  set(PKG_INDEX 0)
  while(PKG_INDEX LESS NUM_PKGS)
    math(EXPR NAME_INDEX "${PKG_INDEX} * 3")
    math(EXPR VERSION_INDEX "${NAME_INDEX} + 1")
    math(EXPR SHA256_INDEX "${NAME_INDEX} + 2")
    
    list(GET MSYS2_PACKAGES ${NAME_INDEX} PKG_NAME)
    list(GET MSYS2_PACKAGES ${VERSION_INDEX} PKG_VERSION)
    list(GET MSYS2_PACKAGES ${SHA256_INDEX} PKG_SHA256)
    
    # Skip placeholder SHA256 values (all hex digits but length might be wrong)
    string(LENGTH "${PKG_SHA256}" SHA256_LEN)
    if(NOT SHA256_LEN EQUAL 64)
      message(WARNING "Skipping ${PKG_NAME}: invalid SHA256 length (${SHA256_LEN}, expected 64)")
      message(WARNING "  Get correct SHA256 from https://packages.msys2.org/package/${PKG_NAME}")
      math(EXPR PKG_INDEX "${PKG_INDEX} + 1")
      continue()
    endif()
    
    fetch_msys2_package(
      PACKAGE_NAME "${PKG_NAME}"
      VERSION "${PKG_VERSION}"
      SHA256 "${PKG_SHA256}"
      DESTINATION "${CMAKE_BINARY_DIR}/third_party_msys2"
    )
    
    list(APPEND ALL_FETCHED_PACKAGES "${PKG_NAME}")
    
    math(EXPR PKG_INDEX "${PKG_INDEX} + 1")
  endwhile()
  
  # Set up FLINT_ROOT from downloaded package if available
  if(ALL_FETCHED_PACKAGES AND "mingw-w64-x86_64-flint" IN_LIST ALL_FETCHED_PACKAGES)
    if(NOT FLINT_ROOT)
      set(FLINT_ROOT "${mingw-w64-x86_64-flint_SOURCE_DIR}/mingw64" CACHE PATH "FLINT root from MSYS2" FORCE)
      message(STATUS "Set FLINT_ROOT to: ${FLINT_ROOT}")
    endif()
  endif()
  
  # Collect all DLLs and copy to dll/ directory
  if(ALL_FETCHED_PACKAGES)
    add_custom_target(copy_msys2_dlls ALL
      COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/dll"
      COMMENT "Preparing DLL directory..."
    )
    
    foreach(PKG_NAME IN LISTS ALL_FETCHED_PACKAGES)
      if(DEFINED ${PKG_NAME}_DLLS)
        add_custom_command(TARGET copy_msys2_dlls POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E copy_if_different
                  ${${PKG_NAME}_DLLS}
                  "${CMAKE_BINARY_DIR}/dll/"
          COMMENT "Copying DLLs from ${PKG_NAME}..."
        )
      endif()
    endforeach()
  endif()
  

endif()

# ============================================================
# C standard
# ============================================================
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)     # GNU extensions used by FLINT headers

# ============================================================
# LTO
# ============================================================
if(DRSOLVE_ENABLE_LTO)
  check_ipo_supported(RESULT _ipo_ok OUTPUT _ipo_msg)
  if(_ipo_ok)
    # Applied per-target below
    message(STATUS "LTO: supported")
  else()
    set(DRSOLVE_ENABLE_LTO OFF)
    message(STATUS "LTO: not supported (${_ipo_msg})")
  endif()
endif()

# ============================================================
# Arch / tuning flags
# ============================================================
if(NOT CMAKE_CROSSCOMPILING)
  check_c_compiler_flag(-march=native _march_native)
  if(_march_native)
    add_compile_options($<$<CONFIG:Release>:-march=native>)
  endif()
else()
  # MinGW cross-build: mirror the flags from Makefile-windows
  add_compile_options(-msse4.1 -mpclmul)
endif()

# ============================================================
# AddressSanitizer
# ============================================================
if(DRSOLVE_ENABLE_ASAN)
  check_c_compiler_flag(-fsanitize=address _asan_ok)
  if(_asan_ok)
    add_compile_options(-fsanitize=address)
    add_link_options(-fsanitize=address)
    message(STATUS "ASan: enabled")
  else()
    message(WARNING "ASan requested but -fsanitize=address is not supported; disabling.")
    set(DRSOLVE_ENABLE_ASAN OFF)
  endif()
endif()

# ============================================================
# Find FLINT  (required)
# ============================================================
set(_flint_include_hints
  "${FLINT_ROOT}/include"
)
set(_flint_lib_hints
  "${FLINT_ROOT}/lib"
)
if(DRSOLVE_WINDOWS)
  list(PREPEND _flint_include_hints
    "${CMAKE_SOURCE_DIR}/mingw/include"
  )
  list(PREPEND _flint_lib_hints
    "${CMAKE_SOURCE_DIR}/mingw/lib"
  )
endif()
list(APPEND _flint_include_hints
  /usr/local/include
  /opt/homebrew/include
  /opt/local/include
)
list(APPEND _flint_lib_hints
  /usr/local/lib
  /opt/homebrew/lib
  /opt/local/lib
)
# Pull in lib dirs from all MSYS2-fetched packages so that GMP/MPFR/etc. are found
foreach(_msys2_dep
    "mingw-w64-x86_64-gmp"
    "mingw-w64-x86_64-mpfr"
    "mingw-w64-x86_64-openblas"
    "mingw-w64-x86_64-gcc-libs"
    "mingw-w64-x86_64-libwinpthread")
  if(DEFINED ${_msys2_dep}_SOURCE_DIR)
    list(APPEND _flint_lib_hints     "${${_msys2_dep}_SOURCE_DIR}/mingw64/lib")
    list(APPEND _flint_include_hints "${${_msys2_dep}_SOURCE_DIR}/mingw64/include")
  endif()
endforeach()

# 1. Try the official CMake config installed by FLINT >= 3.x
find_package(FLINT CONFIG QUIET
  HINTS "${FLINT_ROOT}" "${FLINT_ROOT}/lib/cmake/flint"
)

if(FLINT_FOUND AND TARGET FLINT::flint)
  message(STATUS "FLINT: found via CMake config (${FLINT_VERSION})")
  # Alias to our internal name
  if(NOT TARGET FLINT::FLINT)
    add_library(FLINT::FLINT ALIAS FLINT::flint)
  endif()
else()
  # 2. Manual detection
  if(DRSOLVE_WINDOWS AND FLINT_INCLUDE_DIR)
    if(NOT EXISTS "${FLINT_INCLUDE_DIR}/flint/flint.h" OR NOT EXISTS "${FLINT_INCLUDE_DIR}/gmp.h")
      unset(FLINT_INCLUDE_DIR CACHE)
    endif()
  endif()

  find_path(FLINT_INCLUDE_DIR
    NAMES flint/flint.h
    HINTS ${_flint_include_hints}
  )

  if(DRSOLVE_STATIC)
    find_library(FLINT_LIBRARY
      NAMES libflint.a
      HINTS ${_flint_lib_hints}
    )
  else()
    find_library(FLINT_LIBRARY
      NAMES flint
      HINTS ${_flint_lib_hints}
    )
  endif()

  if(NOT FLINT_INCLUDE_DIR)
    message(FATAL_ERROR
      "FLINT header (flint/flint.h) not found.\n"
      "Install FLINT or set -DFLINT_ROOT=/path/to/flint")
  endif()
  if(NOT FLINT_LIBRARY)
    message(FATAL_ERROR
      "FLINT library not found.\n"
      "Install FLINT or set -DFLINT_ROOT=/path/to/flint")
  endif()

  message(STATUS "FLINT: ${FLINT_LIBRARY}")

  add_library(FLINT::FLINT UNKNOWN IMPORTED GLOBAL)
  set_target_properties(FLINT::FLINT PROPERTIES
    IMPORTED_LOCATION             "${FLINT_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${FLINT_INCLUDE_DIR}"
  )
endif()

# FLINT always needs HAVE_FLINT defined
target_compile_definitions(FLINT::FLINT INTERFACE HAVE_FLINT)

# ============================================================
# Find GMP and MPFR  (FLINT runtime deps)
# ============================================================
find_library(GMP_LIBRARY   NAMES gmp
  HINTS ${_flint_lib_hints})
find_library(MPFR_LIBRARY  NAMES mpfr
  HINTS ${_flint_lib_hints})

if(NOT GMP_LIBRARY)
  message(FATAL_ERROR "libgmp not found (required by FLINT).")
endif()
if(NOT MPFR_LIBRARY)
  message(FATAL_ERROR "libmpfr not found (required by FLINT).")
endif()

# ============================================================
# Find OpenBLAS  (Windows only, optional on others)
# ============================================================
if(DRSOLVE_WINDOWS)
  set(_openblas_hints "")
  if(DEFINED mingw-w64-x86_64-openblas_SOURCE_DIR)
    list(APPEND _openblas_hints "${mingw-w64-x86_64-openblas_SOURCE_DIR}/mingw64/lib")
  endif()
  list(APPEND _openblas_hints "${FLINT_ROOT}/lib")
  
  find_library(OPENBLAS_LIBRARY NAMES openblas
    HINTS ${_openblas_hints})
  if(NOT OPENBLAS_LIBRARY)
    message(WARNING "OpenBLAS not found; linking may fail on Windows.")
  else()
    message(STATUS "OpenBLAS found: ${OPENBLAS_LIBRARY}")
  endif()
endif()

# ============================================================
# Find PML or build from source  (optional)
# ============================================================
set(PML_AVAILABLE FALSE)
set(PML_STATUS "disabled")

if(DRSOLVE_ENABLE_PML)
  set(PML_SOURCE_DIR "${CMAKE_SOURCE_DIR}/pml_det")

  if(EXISTS "${PML_SOURCE_DIR}/src/pml.h")
    message(STATUS "PML: using bundled determinant subset from ${PML_SOURCE_DIR}")

    set(PML_FILTERED_SOURCES
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/nmod_poly_mat_det.c"
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/nmod_poly_mat_weak_popov.c"
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/nmod_poly_mat_utils.c"
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/degrees_pivots_leadingmatrix.c"
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/kernel.c"
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/approximant_basis.c"
      "${PML_SOURCE_DIR}/src/nmod_poly_mat_extra/middle_product.c"
      "${PML_SOURCE_DIR}/src/nmod_mat_poly_extra/nmod_mat_poly_mbasis.c"
      "${PML_SOURCE_DIR}/src/nmod_mat_poly_extra/nmod_mat_poly_mem.c"
      "${PML_SOURCE_DIR}/src/nmod_mat_poly_extra/nmod_mat_poly_arith.c"
      "${PML_SOURCE_DIR}/src/nmod_mat_poly_extra/nmod_mat_poly_set_from.c"
      "${PML_SOURCE_DIR}/src/nmod_mat_poly_extra/nmod_mat_poly_shift.c"
      "${PML_SOURCE_DIR}/src/nmod_mat_extra/left_nullspace.c"
      "${PML_SOURCE_DIR}/src/nmod_extra/nmod_find_root.c"
      "${PML_SOURCE_DIR}/src/nmod_vec_extra/nmod_vec_dot_product.c"
    )

    # Create PML library target
    add_library(pml STATIC ${PML_FILTERED_SOURCES})

    # PML include directories
    target_include_directories(pml PUBLIC
      "${PML_SOURCE_DIR}"
      "${PML_SOURCE_DIR}/src"
    )

    # Link PML with FLINT, GMP, MPFR, and pthread
    target_link_libraries(pml PUBLIC FLINT::FLINT ${GMP_LIBRARY} ${MPFR_LIBRARY})
    if(NOT DRSOLVE_WINDOWS)
      target_link_libraries(pml PUBLIC pthread)
    endif()

    # PML compile definitions - avoid conflicting with FLINT's macros
    target_compile_definitions(pml PUBLIC
      HAVE_PML
      PACKAGE_NAME="PML"
      PACKAGE_VERSION="0.0.5-dev"
      PACKAGE_STRING="PML 0.0.5-dev"
      PACKAGE_BUGREPORT="https://github.com/vneiger/pml/issues/"
      STDC_HEADERS=1
      FLINT_UNROLL_LOOPS=1
      PML_USES_PTHREAD=1
      PML_USES_TLS=1
      PML_WANT_GMP_INTERNALS=1
      PML_WANT_PRETTY_TESTS=1
      PML_HAVE_MACHINE_VECTORS=0
    )

    # Platform-specific definitions
    if(DRSOLVE_WINDOWS)
      target_compile_definitions(pml PUBLIC HAVE_WINDOWS_H=1)
    else()
      target_compile_definitions(pml PUBLIC
        HAVE_ALLOCA_H=1
        HAVE_DLFCN_H=1
        HAVE_MALLOC_H=1
        HAVE_SYS_PARAM_H=1
        HAVE_UNISTD_H=1
        HAVE_ALIGNED_ALLOC=1
      )
    endif()

    # Create alias target
    add_library(PML::PML ALIAS pml)
    set(PML_AVAILABLE TRUE)
    set(PML_STATUS "bundled subset: ${PML_SOURCE_DIR}")
  else()
    message(STATUS "PML: bundled source not found under ${PML_SOURCE_DIR}; trying system installation...")

    find_path(PML_INCLUDE_DIR
      NAMES pml.h pml/pml.h
      HINTS "${PML_ROOT}/include"
            /usr/local/include
            /opt/homebrew/include
    )

    if(DRSOLVE_STATIC)
      find_library(PML_LIBRARY NAMES libpml.a
        HINTS "${PML_ROOT}/lib" /usr/local/lib /opt/homebrew/lib)
    else()
      find_library(PML_LIBRARY NAMES pml
        HINTS "${PML_ROOT}/lib" /usr/local/lib /opt/homebrew/lib)
    endif()

    if(PML_INCLUDE_DIR AND PML_LIBRARY)
      set(_pml_extra_inc "")
      if(EXISTS "${PML_INCLUDE_DIR}/nmod_poly_mat_utils.h")
        set(_pml_extra_inc "")
      elseif(EXISTS "${PML_INCLUDE_DIR}/pml/nmod_poly_mat_utils.h")
        set(_pml_extra_inc "${PML_INCLUDE_DIR}/pml")
      else()
        message(STATUS "PML: nmod_poly_mat_utils.h not found under ${PML_INCLUDE_DIR}")
        set(PML_LIBRARY "")
      endif()

      if(PML_LIBRARY)
        if(NOT EXISTS "${PML_INCLUDE_DIR}/config.h" AND
           NOT EXISTS "${PML_INCLUDE_DIR}/pml/config.h")
          message(STATUS "PML: config.h not found under ${PML_INCLUDE_DIR}")
          set(PML_LIBRARY "")
        endif()
      endif()
    endif()

    if(PML_INCLUDE_DIR AND PML_LIBRARY)
      set(PML_AVAILABLE TRUE)
      set(PML_STATUS "${PML_LIBRARY}")
      message(STATUS "PML: using installed library ${PML_LIBRARY}")

      add_library(PML::PML UNKNOWN IMPORTED GLOBAL)
      set_target_properties(PML::PML PROPERTIES
        IMPORTED_LOCATION             "${PML_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PML_INCLUDE_DIR}"
        INTERFACE_COMPILE_DEFINITIONS "HAVE_PML"
      )
      if(_pml_extra_inc)
        set_property(TARGET PML::PML APPEND PROPERTY
          INTERFACE_INCLUDE_DIRECTORIES "${_pml_extra_inc}")
      endif()
    else()
      message(STATUS "PML: not found (optional) -- building without PML support")
    endif()
  endif()
endif()

# ============================================================
# Find OpenMP
# ============================================================
if(DRSOLVE_ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS C)
  if(OpenMP_C_FOUND)
    message(STATUS "OpenMP: ${OpenMP_C_FLAGS}")
  else()
    message(STATUS "OpenMP: not found")
  endif()
endif()

# ============================================================
# System libraries list  (per platform)
# ============================================================
set(_sys_libs "")
if(APPLE)
  list(APPEND _sys_libs ${MPFR_LIBRARY} ${GMP_LIBRARY} m)
elseif(DRSOLVE_WINDOWS)
  list(APPEND _sys_libs ${MPFR_LIBRARY} ${GMP_LIBRARY} m pthread)
  if(OPENBLAS_LIBRARY)
    list(APPEND _sys_libs ${OPENBLAS_LIBRARY})
  endif()
else()  # Linux
  list(APPEND _sys_libs ${MPFR_LIBRARY} ${GMP_LIBRARY} m pthread stdc++)
endif()

# ============================================================
# Common link-library list  (shared between all targets)
# ============================================================
set(_drsolve_deps)
if(PML_AVAILABLE)
  list(APPEND _drsolve_deps PML::PML)
endif()
list(APPEND _drsolve_deps FLINT::FLINT)
if(OpenMP_C_FOUND)
  list(APPEND _drsolve_deps OpenMP::OpenMP_C)
endif()
list(APPEND _drsolve_deps ${_sys_libs})

# ============================================================
# Math source list  (shared between library and exe targets)
# ============================================================
set(MATH_SOURCES
  src/dixon_complexity.c
  src/dixon_recursive.c
  src/dixon_flint.c
  src/dixon_interface_flint.c
  src/dixon_test.c
  src/dixon_with_ideal_reduction.c
  src/fq_mat_det.c
  src/macaulay_flint.c
  src/fq_mpoly_mat_det.c
  src/fq_multivariate_interpolation.c
  src/fq_mvpoly.c
  src/fq_nmod_roots.c
  src/fq_poly_mat_det.c
  src/fq_sparse_interpolation.c
  src/fq_unified_interface.c
  src/fmpq_acb_roots.c
  src/gf2n_mpoly.c
  src/gf2n_field.c
  src/gf2n_poly.c
  src/large_prime_system_solver.c
  src/polynomial_system_solver.c
  src/rational_system_solver.c
  src/resultant_with_ideal_reduction.c
  src/unified_mpoly_det.c
  src/unified_mpoly_interface.c
  src/unified_mpoly_resultant.c
)

# Convenience function: apply standard settings to a drsolve target
function(drsolve_target_common _tgt)
  target_include_directories(${_tgt} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/pml_det>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/drsolve>
  )
  target_link_libraries(${_tgt} PUBLIC ${_drsolve_deps})
  target_compile_options(${_tgt} PRIVATE
    $<$<CONFIG:Release>:-O3>
    $<$<CONFIG:Debug>:-O0 -g>
  )
  if(DRSOLVE_ENABLE_LTO)
    set_property(TARGET ${_tgt} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  endif()
endfunction()

# ============================================================
# Static library: libdrsolve.a  (libdrsolve-1.a on Windows)
# ============================================================
if(DRSOLVE_BUILD_STATIC_LIB OR DRSOLVE_STATIC)
  add_library(drsolve_static STATIC ${MATH_SOURCES})
  drsolve_target_common(drsolve_static)
  if(DRSOLVE_WINDOWS)
    set_target_properties(drsolve_static PROPERTIES OUTPUT_NAME "libdrsolve-1")
  else()
    set_target_properties(drsolve_static PROPERTIES OUTPUT_NAME "drsolve")
  endif()
  set_target_properties(drsolve_static PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
  )
endif()

# ============================================================
# Shared library: libdrsolve.so / .dylib / libdrsolve-1.dll
# ============================================================
if(DRSOLVE_BUILD_SHARED_LIB AND NOT DRSOLVE_STATIC)
  add_library(drsolve_shared SHARED ${MATH_SOURCES})
  drsolve_target_common(drsolve_shared)

  if(DRSOLVE_WINDOWS)
    # Windows DLL with import library
    set_target_properties(drsolve_shared PROPERTIES
      OUTPUT_NAME            "libdrsolve-1"
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/dll"
      ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"   # .dll.a import lib
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    )
  else()
    set_target_properties(drsolve_shared PROPERTIES
      OUTPUT_NAME            "drsolve"
      LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
      VERSION                "${PROJECT_VERSION}"
      SOVERSION              "1"
    )
  endif()

  # macOS: set install name so binaries find the dylib at runtime
  if(APPLE)
    set_target_properties(drsolve_shared PROPERTIES
      INSTALL_NAME_DIR "@rpath"
      BUILD_WITH_INSTALL_RPATH FALSE
    )
  endif()
endif()

# ============================================================
# Main executable: drsolve
# On Windows it lives in bin/ as drsolve_cli_real.exe;
# the top-level drsolve.exe is a thin launcher (built below).
# ============================================================
# The original Makefile compiles all sources together in one
# translation unit for maximum inlining ("unity" compilation).
# With CMake IPO/LTO enabled this is equivalent at link-time.
# We replicate the unity approach via a separate object library
# that is linked privately into the executable.

add_executable(drsolve_exe drsolve.c ${MATH_SOURCES})
drsolve_target_common(drsolve_exe)

if(DRSOLVE_WINDOWS)
  set_target_properties(drsolve_exe PROPERTIES
    OUTPUT_NAME              "drsolve_cli_real"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
  )
else()
  set_target_properties(drsolve_exe PROPERTIES
    OUTPUT_NAME "drsolve"
  )

  # RPATH: find shared libs relative to binary and from install prefix
  if(APPLE)
    set_target_properties(drsolve_exe PROPERTIES
      INSTALL_RPATH "@loader_path/../lib"
      BUILD_WITH_INSTALL_RPATH FALSE
      BUILD_RPATH "${CMAKE_BINARY_DIR}/lib"
    )
  else()
    # Linux: embed RPATH to bundled lib dir and FLINT/PML locations
    set(_rpath_dirs "${CMAKE_BINARY_DIR}/lib;\$ORIGIN/../lib")
    get_filename_component(_flint_lib_dir "${FLINT_LIBRARY}" DIRECTORY)
    list(APPEND _rpath_dirs "${_flint_lib_dir}")
    if(PML_AVAILABLE AND PML_LIBRARY)
      get_filename_component(_pml_lib_dir "${PML_LIBRARY}" DIRECTORY)
      list(APPEND _rpath_dirs "${_pml_lib_dir}")
    endif()
    string(REPLACE ";" ":" _rpath_str "${_rpath_dirs}")
    set_target_properties(drsolve_exe PROPERTIES
      INSTALL_RPATH           "\$ORIGIN/../lib"
      BUILD_RPATH             "${_rpath_str}"
      BUILD_WITH_INSTALL_RPATH FALSE
    )
  endif()
endif()

# Static-link mode: add linker flags for duplicate-symbol tolerance
# (FLINT and PML share some symbol names intentionally)
if(DRSOLVE_STATIC)
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    target_link_options(drsolve_exe PRIVATE
      -Wl,--allow-multiple-definition
      -no-pie
    )
  endif()
  # macOS: Apple ld resolves duplicates by first-definition-wins; no flag needed.
  # (In CI, nmedit is used to localise FLINT symbols before linking.)
endif()

# ============================================================
# Windows: launcher (drsolve.exe) and GUI (drsolve_win_gui.exe)
# ============================================================
if(DRSOLVE_WINDOWS AND DRSOLVE_BUILD_GUI)

  # Thin launcher: sets DLL path then spawns bin/drsolve_cli_real.exe
  add_executable(drsolve_launcher win_gui/dixon_launcher.c)
  target_compile_definitions(drsolve_launcher PRIVATE _WIN32_WINNT=0x0601)
  set_target_properties(drsolve_launcher PROPERTIES
    OUTPUT_NAME              "drsolve"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  )

  # Windows GUI frontend
  add_executable(drsolve_win_gui WIN32
    win_gui/dixon_gui.c
    win_gui/dixon_wrapper.c
  )
  target_compile_definitions(drsolve_win_gui PRIVATE _WIN32_WINNT=0x0601)
  target_link_libraries(drsolve_win_gui PRIVATE comctl32 comdlg32 shlwapi)
  set_target_properties(drsolve_win_gui PROPERTIES
    OUTPUT_NAME              "drsolve_win_gui"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  )

  # Copy runtime DLLs:
  # - For native MSYS2 builds: don't copy (deps come from pacman)
  # - For cross-compile with bundled deps: copy from runtime/
  # - For cross-compile without bundled deps: DLLs are copied via copy_msys2_dlls target
  if(DRSOLVE_USE_BUNDLED_DEPS)
    file(GLOB _runtime_dlls "${CMAKE_SOURCE_DIR}/runtime/*.dll")
    if(_runtime_dlls)
      add_custom_command(TARGET drsolve_exe POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/dll"
        COMMENT "Copying bundled runtime DLLs → dll/"
        VERBATIM
      )
      add_custom_command(TARGET drsolve_exe POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${_runtime_dlls} "${CMAKE_BINARY_DIR}/dll/"
        COMMENT "Copying runtime DLLs..."
        VERBATIM
      )
    endif()
  elseif(TARGET copy_msys2_dlls)
    # Make sure copy_msys2_dlls runs before building drsolve_exe
    add_dependencies(drsolve_exe copy_msys2_dlls)
    message(STATUS "Using MSYS2-fetched DLLs (via copy_msys2_dlls)")
  else()
    message(STATUS "Native/cross build: runtime DLLs will be provided by system/toolchain")
  endif()
endif()

# ============================================================
# Attack programs  (optional, ../Attack directory)
# ============================================================
if(DRSOLVE_BUILD_ATTACK)
  cmake_path(GET CMAKE_SOURCE_DIR PARENT_PATH _attack_parent)
  set(_attack_dir "${_attack_parent}/Attack")
  if(EXISTS "${_attack_dir}" AND IS_DIRECTORY "${_attack_dir}")
    file(GLOB_RECURSE _attack_srcs CONFIGURE_DEPENDS "${_attack_dir}/*.c")
    foreach(_src IN LISTS _attack_srcs)
      # Skip Jupyter checkpoint directories
      if(_src MATCHES "\\.ipynb_checkpoints/")
        continue()
      endif()
      cmake_path(GET _src STEM   _aname)
      cmake_path(GET _src PARENT_PATH _adir)
      set(_tgt "attack_${_aname}")
      add_executable(${_tgt} EXCLUDE_FROM_ALL ${_src} ${MATH_SOURCES})
      drsolve_target_common(${_tgt})
      set_target_properties(${_tgt} PROPERTIES
        OUTPUT_NAME              "${_aname}"
        RUNTIME_OUTPUT_DIRECTORY "${_adir}"
      )
    endforeach()
    message(STATUS "Attack programs: found at ${_attack_dir}")
  else()
    message(STATUS "Attack programs: ${_attack_dir} not found -- skipping")
  endif()
endif()

# ============================================================
# Install rules
# ============================================================

# Executable
install(TARGETS drsolve_exe
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Static library
if(DRSOLVE_BUILD_STATIC_LIB OR DRSOLVE_STATIC)
  install(TARGETS drsolve_static
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()

# Shared library
if(DRSOLVE_BUILD_SHARED_LIB AND NOT DRSOLVE_STATIC)
  install(TARGETS drsolve_shared
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}   # Windows DLL → bin/
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}   # Windows import lib → lib/
  )
  # Linux: run ldconfig after install
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    install(CODE "
      execute_process(
        COMMAND ldconfig \"${CMAKE_INSTALL_FULL_LIBDIR}\"
        ERROR_QUIET
      )")
  endif()
endif()

# Headers
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/drsolve
  FILES_MATCHING PATTERN "*.h"
)

# ============================================================
# CTest
# ============================================================
enable_testing()

# Helper to register a test and route output files to a tmp dir
set(_test_wd "${CMAKE_BINARY_DIR}")

add_test(NAME drsolve_help
  COMMAND drsolve_exe --help
  WORKING_DIRECTORY ${_test_wd}
)
set_tests_properties(drsolve_help PROPERTIES PASS_REGULAR_EXPRESSION ".")

add_test(NAME drsolve_basic
  COMMAND drsolve_exe "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_basic_rational
  COMMAND drsolve_exe "x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z" "x,y" 0
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_solve
  COMMAND drsolve_exe --solve "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_solve_linear
  COMMAND drsolve_exe --solve "x+y-3, x-y+1" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_comp
  COMMAND drsolve_exe --comp "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_comp_short
  COMMAND drsolve_exe -c "x^2+y^2+1, x*y+z, x+y+z^2" "x,y" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_comp_omega
  COMMAND drsolve_exe --comp --omega 2.81
    "x^4+y^4+z^4+w^4+1, x^3*y+z+1, x+y^3+z^2+w, x*y*z*w+1"
    "x,y,z" 65537
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_random
  COMMAND drsolve_exe --random "[3,3,2]" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_random_solve
  COMMAND drsolve_exe -r --solve "[2]*3" 257
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_random_comp
  COMMAND drsolve_exe -r --comp "[3,2]*2" 65537
  WORKING_DIRECTORY ${_test_wd}
)
add_test(NAME drsolve_extension_field
  COMMAND drsolve_exe --silent "x+y^2+t, x*y+t*y+1" "x" "2^8"
  WORKING_DIRECTORY ${_test_wd}
)

# ============================================================
# Configuration summary
# ============================================================
message(STATUS "")
message(STATUS "=== DRsolve ${PROJECT_VERSION} — Build Configuration ===")
message(STATUS "  Platform        : ${CMAKE_SYSTEM_NAME} (${CMAKE_SYSTEM_PROCESSOR})")
message(STATUS "  Build type      : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Compiler        : ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "  FLINT           : ${FLINT_LIBRARY}")
if(PML_AVAILABLE)
  message(STATUS "  PML             : ${PML_STATUS}")
else()
  message(STATUS "  PML             : disabled")
endif()
if(OpenMP_C_FOUND)
  message(STATUS "  OpenMP          : ${OpenMP_C_VERSION}")
else()
  message(STATUS "  OpenMP          : not found")
endif()
message(STATUS "  LTO             : ${DRSOLVE_ENABLE_LTO}")
message(STATUS "  ASan            : ${DRSOLVE_ENABLE_ASAN}")
message(STATUS "  Static linking  : ${DRSOLVE_STATIC}")
message(STATUS "  Install prefix  : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "==============================================")
