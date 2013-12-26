# SU2_EDU autoconf file
# Stanford University Aerospace Design Lab
#
# Created: 		10/26/12
# Last updated:         10/11/13
#
# Created by: 		M. Colonno  10/26/12
# Updated by: 		T. Economon  10/30/13
#
##########################

AC_INIT([SU2_EDU], [1.0.0], [susquared-dev@lists.stanford.edu],
        [SU2_EDU], [http://su2.stanford.edu])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([subdir-objects])

# optional bindings
AC_ARG_WITH(MPI,
    	AS_HELP_STRING([--with-MPI[=ARG]], [Build parallel version with MPI tools, ARG = MPI C++ compiler]), 
	[with_MPI=$withval], [with_MPI="no"])

###########################

have_MPI="no"
if test "$with_MPI" != "no"
then
	CXX="$with_MPI"	
	have_MPI="yes"
else
	CXXFLAGS=""$CXXFLAGS" -DNO_MPI"
fi

###########################

AC_PROG_CXX
AC_PROG_RANLIB
AC_OUTPUT(Makefile)

abs_srcdir=`(cd $srcdir && pwd)`
abs_prefix=`(cd $prefix && pwd)`

AC_MSG_RESULT([

   _____. __    __   ___       _______  _______   __    __  
  /     ||  |  |  | |__ \     |   ____||       \ |  |  |  | 
 |   (--`|  |  |  |    ) |    |  |__   |  .--.  ||  |  |  | 
  \   \  |  |  |  |   / /     |   __|  |  |  |  ||  |  |  | 
.--)   | |  `--'  |  / /_     |  |____ |  '--'  ||  `--'  | 
|_____/   \______/  |____|    |_______||_______/  \______/  

Release 1.0.0 (su2.stanford.edu)

-------------------------------------------------------------------------
| Stanford University Unstructured (SU2).                               | 
| Copyright (C) 2012-2013 Aerospace Design Laboratory (ADL).            |
| SU2 is distributed in the hope that it will be useful,                |
| but WITHOUT ANY WARRANTY; without even the implied warranty of        |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      |
| Lesser General Public License (version 2.1) for more details.         |
-------------------------------------------------------------------------

Build Configuration Summary:

    Source code location: ${abs_srcdir}
    Install location: 	  ${abs_prefix}
    Version:              ${VERSION}
    Compiler:             ${CXX}
    Preprocessor flags:   ${CPPFLAGS}
    Compiler flags:       ${CXXFLAGS}
    Linker flags:         ${LDFLAGS}
    MPI support:          $have_MPI

Please be sure to add the [\$]SU2_EDU_RUN environment variables, and update your [\$]PATH with [\$]SU2_EDU_RUN

Based on the input to this configuration, add these lines to your .bashrc file -

export SU2_EDU_RUN="${abs_prefix}/bin"
export PATH=[\$]PATH:[\$]SU2_EDU_RUN

])
