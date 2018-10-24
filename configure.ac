dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(poing2, 0.1)
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects])

AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_CC_C99
AC_PROG_RANLIB
AC_OPENMP
AC_GNU_SOURCE

AM_MISSING_PROG([GPERF], [gperf])
AC_SEARCH_LIBS([cos], [m], [], [AC_MSG_ERROR([Could not find cosine function])])
AC_SEARCH_LIBS([clock_gettime], [rt],
               [
                AC_DEFINE([HAVE_CLOCK_GETTIME], [1], [Used for profiler])
                AM_CONDITIONAL([HAVE_CLOCK_GETTIME_AM], [true])
               ],
               [
                AC_MSG_WARN([Could not find clock_gettime function])
                AM_CONDITIONAL([HAVE_CLOCK_GETTIME_AM], [false])
               ])

dnl Set the HAVE_OPENMP flag if using openmp
AS_IF([test -n "$OPENMP_CFLAGS"],
      [AC_DEFINE([HAVE_OPENMP], [1], [Check for OpenMP])],
      [])

AC_CONFIG_FILES(Makefile)
AC_OUTPUT

