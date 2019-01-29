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

AC_ARG_ENABLE(64bit_ints,
[AS_HELP_STRING([--enable-64bit-ints],[Use 64 bit integers instead of 32 bit])])
if test "x$enable_64bit_ints" == "xyes" ; then
    AC_TYPE_INT64_T
	CXXFLAGS+=" -DGMGW_INT64"
else
    AC_TYPE_INT32_T
	CXXFLAGS+=" -DGMGW_INT32"
fi

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

AC_ARG_WITH(tau,
	[AS_HELP_STRING([--with-tau],[Provide a directory where the TAU installation can be found])],,with_tau="no")
AC_ARG_WITH(netcdf,
	[AS_HELP_STRING([--with-netcdf],[Provide a directory where the NETCDF libraies can be found])],,with_netcdf="no")
AC_ARG_WITH(mpi,
	[AS_HELP_STRING([--with-mpi],[Provide a directory where the MPI io libraies can be found])],,with_mpi="no")

AC_SUBST(TAU_DEFINES)
AC_SUBST(TAU_INCLUDE)
AC_SUBST(TAU_LIBS)
AC_SUBST(TAU_OBJ)
AC_SUBST(NETCDF_LIB)
AC_SUBST(MPI_LIB)

if test "x$with_tau" == xno; then
   true
else
  TAU_DEFINES="-DHAVE_TAU"
  TAU_INCLUDE="-I$with_tau/include"
  TAU_LIBS="-L$with_tau/lib -lprimgrid_io -lgrid_util -lstream_interface -ltau_structs -lnetcdf_interface -lqueue -lcutil"
  TAU_OBJ="TAUFileWrapper.o"
  if test "x$with_netcdf" == xno; then
    NETCDF_LIB="-lnetcdf"
  else
    NETCDF_LIB="-L$with_netcdf -lnetcdf"
  fi
  if test "x$with_mpi" == xno; then
    MPI_LIB="-lmpi"
  else
    MPI_LIB="-L$with_mpi -lmpi"
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
