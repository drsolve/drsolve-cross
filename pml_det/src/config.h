/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/*
    Copyright (C) 2019 Seung Gyu Hyun, Vincent Neiger, Eric Schost
    Copyright (C) 2025 Vincent Neiger, Eric Schost, Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 3 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Define to locally unroll some loops */
#define FLINT_UNROLL_LOOPS 1

/* Define to 1 if you have the `aligned_alloc' function. */
#define HAVE_ALIGNED_ALLOC 1

/* Define to 1 if you have the <alloca.h> header file. */
#define HAVE_ALLOCA_H 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the <fenv.h> header file. */
#define HAVE_FENV_H 1

/* Define to 1 if you have the <float.h> header file. */
#define HAVE_FLOAT_H 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <pthread_np.h> header file. */
/* #undef HAVE_PTHREAD_NP_H */

/* Define to 1 if you have the <stdarg.h> header file. */
#define HAVE_STDARG_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/param.h> header file. */
#define HAVE_SYS_PARAM_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define to 1 if you have the `_aligned_malloc' function. */
/* #undef HAVE__ALIGNED_MALLOC */

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://github.com/vneiger/pml/issues/"

/* Define to the full name of this package. */
#define PACKAGE_NAME "PML"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "PML 0.0.5-dev"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "pml"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.0.5-dev"

/* Define if system is big endian. */
/* #undef PML_BIG_ENDIAN */

/* Define to enable coverage. */
/* #undef PML_COVERAGE */

/* Define to enable reentrant. */
/* #undef PML_REENTRANT */

/* Define to enable BLAS. */
/* #undef PML_USES_BLAS */

/* Define to enable the use of pthread. */
#define PML_USES_PTHREAD 1

/* Define to enable thread-local storage. */
#define PML_USES_TLS 1

/* Define to enable use of asserts. */
/* #undef PML_WANT_ASSERT */

/* Define to enable use of GMP internals. */
#define PML_WANT_GMP_INTERNALS 1

/* Define to enable pretty printing for tests. */
#define PML_WANT_PRETTY_TESTS 1

/* Define the following to what diagnostic pragmas your compiler allows.
   These are used to silence certain warnings. */
#ifndef DIAGNOSTIC_PUSH
#define DIAGNOSTIC_PUSH _Pragma("GCC diagnostic push")
#endif
#ifndef DIAGNOSTIC_POP
#define DIAGNOSTIC_POP _Pragma("GCC diagnostic pop")
#endif
#ifndef DIAGNOSTIC_IGNORE_INCOMPATIBLE_FUNCTION_POINTER_TYPES
#define DIAGNOSTIC_IGNORE_INCOMPATIBLE_FUNCTION_POINTER_TYPES
#endif
#ifndef DIAGNOSTIC_IGNORE_DISCARDED_QUALIFIERS
#define DIAGNOSTIC_IGNORE_DISCARDED_QUALIFIERS _Pragma("GCC diagnostic ignored \"-Wdiscarded-qualifiers\"")
#endif
#ifndef DIAGNOSTIC_IGNORE_FORMAT
#define DIAGNOSTIC_IGNORE_FORMAT _Pragma("GCC diagnostic ignored \"-Wformat\"")
#endif
#ifndef DIAGNOSTIC_IGNORE_DANGLING_POINTER
#define DIAGNOSTIC_IGNORE_DANGLING_POINTER
#endif
#ifndef DIAGNOSTIC_IGNORE_CAST_FUNCTION_TYPE
#define DIAGNOSTIC_IGNORE_CAST_FUNCTION_TYPE _Pragma("GCC diagnostic ignored \"-Wcast-function-type\"")
#endif
#ifndef DIAGNOSTIC_IGNORE_OVERLENGTH_STRINGS
#define DIAGNOSTIC_IGNORE_OVERLENGTH_STRINGS _Pragma("GCC diagnostic ignored \"-Woverlength-strings\"")
#endif
#ifndef DIAGNOSTIC_IGNORE_UNUSED_VARIABLE
#define DIAGNOSTIC_IGNORE_UNUSED_VARIABLE _Pragma("GCC diagnostic ignored \"-Wunused-variable\"")
#endif
#ifndef DIAGNOSTIC_IGNORE_MAYBE_UNINITIALIZED
#define DIAGNOSTIC_IGNORE_MAYBE_UNINITIALIZED _Pragma("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
#endif

/* Define the following to what optimization pragmas your compiler allows. */
#ifndef PUSH_OPTIONS
#define PUSH_OPTIONS _Pragma("GCC push_options")
#endif
#ifndef POP_OPTIONS
#define POP_OPTIONS
#endif
#ifndef OPTIMIZE_O2
#define OPTIMIZE_O2 _Pragma("GCC optimize (\"O2\")")
#endif
#ifndef OPTIMIZE_OSIZE
#define OPTIMIZE_OSIZE _Pragma("GCC optimize (\"Os\")")
#endif
#ifndef OPTIMIZE_UNROLL_LOOPS
#define OPTIMIZE_UNROLL_LOOPS _Pragma("GCC optimize (\"unroll-loops\")")
#endif

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif
