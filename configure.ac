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
AM_INIT_AUTOMAKE([subdir-objects])

###########################
# Force the NO_MPI directive and -O3 optimization

CXXFLAGS=""$CXXFLAGS" -DNO_MPI -O3"

###########################

AC_PROG_CXX
AC_PROG_RANLIB
AC_OUTPUT(Makefile)

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


The build has been configured. Run the 'make' command to compile SU2_EDU.

After building, the SU2_EDU executable can be found in the bin/ directory.

Detailed instructions for running SU2_EDU can be found in the README file.

])
