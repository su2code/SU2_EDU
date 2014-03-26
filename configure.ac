# SU2_EDU autoconf file
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated:         02/04/14
#
# Created by: 		M. Colonno
# Updated by: 		T. Economon
#
##########################

AC_INIT([SU2_EDU], [1.1.0], [susquared-dev@lists.stanford.edu],
        [SU2_EDU], [http://su2.stanford.edu])
AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([subdir-objects])  # [SU^2], [2.0]
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES(yes)])

# AC_CONFIG_HEADERS([config.h])
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_RANLIB
AC_PROG_MKDIR_P

# --------------------------------------------------------------
# Check for important type sizes
# --------------------------------------------------------------
AC_CHECK_SIZEOF(short int)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(unsigned int)
AC_CHECK_SIZEOF(long int)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(void *)

# --------------------------------------------------------------
# Build up contributed packages that we will build from source
# --------------------------------------------------------------
su2_externals_INCLUDES=""
su2_externals_LIBS=""

# Metis
CONFIGURE_METIS
if (test $enablemetis = yes); then
  su2_externals_INCLUDES="$METIS_INCLUDE $su2_externals_INCLUDES"
  su2_externals_LIBS="$METIS_LIB $su2_externals_LIBS"
else
  su2_externals_INCLUDES="$su2_externals_INCLUDES"
fi
AM_CONDITIONAL(BUILD_METIS, test x$enablemetis = xyes)
AC_CONFIG_FILES([externals/metis/Makefile])

AC_SUBST([su2_externals_INCLUDES])
AC_SUBST([su2_externals_LIBS])

AC_ARG_WITH(Metis-include,
        AS_HELP_STRING([--with-Metis-include[=ARG]], [Metis include directory, ARG = path to metis.h]),
        [with_Metis_include=$withval], [with_Metis_include="no"])

###########################
# Force the NO_MPI directive and -O3 optimization

CXXFLAGS=""$CXXFLAGS" -DNO_MPI -O3"

###########################

AC_PROG_CXX
AC_PROG_RANLIB
AC_OUTPUT(Makefile \
        externals/Makefile)

abs_srcdir=`(cd $srcdir && pwd)`

##########################

abs_prefix=`(cd $prefix && pwd)`

#########################


AC_MSG_RESULT([

   _____. __    __   ___       _______  _______   __    __  
  /     ||  |  |  | |__ \     |   ____||       \ |  |  |  | 
 |   (--`|  |  |  |    ) |    |  |__   |  .--.  ||  |  |  | 
  \   \  |  |  |  |   / /     |   __|  |  |  |  ||  |  |  | 
.--)   | |  `--'  |  / /_     |  |____ |  '--'  ||  `--'  | 
|_____/   \______/  |____|    |_______||_______/  \______/  

Release 1.1.0 (su2.stanford.edu)

-------------------------------------------------------------------------
| SU2, Copyright (C) 2012-2014 Aerospace Design Laboratory (ADL).       |
| SU2 is distributed in the hope that it will be useful,                |
| but WITHOUT ANY WARRANTY; without even the implied warranty of        |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |
| Lesser General Public License (version 2.1) for more details.         |
-------------------------------------------------------------------------

Build Configuration Summary:

    Source code location: ${abs_srcdir}
    Version:              ${VERSION}
    Compiler:             ${CXX}
    Compiler flags:       ${CXXFLAGS}
    Metis support:        $enablemetis

    External includes:    $su2_externals_INCLUDES
    External libs:        $su2_externals_LIBS


The build has been configured. Run the 'make' command to compile SU2_EDU.

After building, the SU2_EDU executable can be found in the bin/ directory.

Detailed instructions for running SU2_EDU can be found in the README file.

])
