# generated automatically by aclocal 1.15 -*- Autoconf -*-

# Copyright (C) 1996-2014 Free Software Foundation, Inc.

# This file is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.

m4_ifndef([AC_CONFIG_MACRO_DIRS], [m4_defun([_AM_CONFIG_MACRO_DIRS], [])m4_defun([AC_CONFIG_MACRO_DIRS], [_AM_CONFIG_MACRO_DIRS($@)])])
# ===========================================================================
#    http://www.gnu.org/software/autoconf-archive/ax_require_defined.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_REQUIRE_DEFINED(MACRO)
#
# DESCRIPTION
#
#   AX_REQUIRE_DEFINED is a simple helper for making sure other macros have
#   been defined and thus are available for use.  This avoids random issues
#   where a macro isn't expanded.  Instead the configure script emits a
#   non-fatal:
#
#     ./configure: line 1673: AX_CFLAGS_WARN_ALL: command not found
#
#   It's like AC_REQUIRE except it doesn't expand the required macro.
#
#   Here's an example:
#
#     AX_REQUIRE_DEFINED([AX_CHECK_LINK_FLAG])
#
# LICENSE
#
#   Copyright (c) 2014 Mike Frysinger <vapier@gentoo.org>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved. This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_REQUIRE_DEFINED], [dnl
  m4_ifndef([$1], [m4_fatal([macro ]$1[ is not defined; is a m4 file missing?])])
])dnl AX_REQUIRE_DEFINED

m4_include([m4/m4_ax_cxx_compile_stdcxx.m4])
#m4_include([m4/m4_ax_cxx_compile_stdcxx_11.m4])
