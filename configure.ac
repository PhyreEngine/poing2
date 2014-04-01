dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(poing2, 0.1)
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_SYSTEM
AM_INIT_AUTOMAKE([-Wall -Werror subdir-objects])

AC_PROG_CC
AC_PROG_CC_C99
AC_PROG_RANLIB
AC_OPENMP
AC_GNU_SOURCE

AM_MISSING_PROG([GPERF], [gperf])
AC_SEARCH_LIBS([cos], [m], [], [AC_MSG_ERROR([Could not find cosine function])])

AC_CONFIG_FILES(Makefile)
AC_OUTPUT

