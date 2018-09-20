AC_PREREQ([2.69])
AC_INIT([analyzeUnstrMesh], [1.0.0], [cfog@mech.ubc.ca])
AC_CONFIG_SRCDIR([analyzeVolMesh.cxx])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])

# Checks for programs.
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11

AC_ARG_ENABLE([debug],
[AS_HELP_STRING([--enable-debug],[Turn on warnings and debugging support])],
[CXXFLAGS+=" -g"],
[CXXFLAGS+=" -O3"])


# Checks for libraries.
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [sqrt])

AC_SUBST(INT_SIZE)
AC_ARG_ENABLE(64bit-ints,
[AS_HELP_STRING([--enable-64bits-ints],[Use 64 bit integers instead of 32 bit])],[CXXFLAGS+=" -DGMGW_INT64"],[CXXFLAGS+=" -DGMGW_INT32"])
echo $CXXFLAGS

AC_ARG_WITH(geode,
	[AS_HELP_STRING([--with-geode],[Provide a directory where the Pointwise GEODE geometry kernel can be found])],,with_geode="no")
AC_ARG_WITH(geode-arch,
	[AS_HELP_STRING([--with-geode-arch],[Provide an architecture identifier for GEODE])],,with_geode_arch="no")

AC_SUBST(GEODE_DEFINES)
AC_SUBST(GEODE_LFLAGS)
AC_SUBST(GEODE_INCLUDE)
AC_SUBST(GEODE_LIBS)

if test "x$with_geode" == xno; then
   true
else
  if test "x$with_geode_arch" == xno; then
    AC_MSG_ERROR("You must also specify --with-geode-arch!")
  else
    GEODE_DEFINES="-DHAVE_GEODE -DGE_SHARED_LINKAGE -DGE_EXPORT_SYMBOLS"
    GEODE_LFLAGS="-Wl,-O1 -static"
    GEODE_INCLUDE="-I$with_geode/include"
    GEODE_LIBS="-L$with_geode/lib/$with_geode_arch -lnmb -lgeom -lvmath"
  fi
fi

# FIXME: Replace `main' with a function in `-lpthread':
AC_CHECK_LIB([pthread], [pthread_detach])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt strstr sysinfo])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
