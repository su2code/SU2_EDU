# AC_C_OPENMP
# -----------
# Check which options need to be passed to the C compiler to support OpenMP.
# Set the OPENMP_CFLAGS variable to these options.
# The options are necessary at compile time (so the #pragmas are understood)
# and at link time (so the appropriate library is linked with).
# This macro takes care to not produce redundant options if $CC $CFLAGS already
# supports OpenMP. It also is careful to not pass options to compilers that
# misinterpret them; for example, most compilers accept "-openmp" and create
# an output file called 'penmp' rather than activating OpenMP support.

# Adapted from (2007-05-17) Bruno Haible, http://www.spinics.net/lists/ac/msg05891.html
 AC_DEFUN([AC_C_OPENMP],
 [
   AC_ARG_ENABLE(openmp,
     [  --disable-openmp        do not use OpenMP],
     [enableopenmp="$enableval"],
     [enableopenmp=yes])
   AC_MSG_RESULT(<<< Configuring for OpenMP Support >>>)
   OPENMP_CFLAGS=
   if test "$enableopenmp" = yes; then
     AC_MSG_CHECKING([for $CC option to support OpenMP])
     AC_CACHE_VAL([ac_cv_prog_cc_openmp], [
       ac_cv_prog_cc_openmp=no
       AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
         ])], [ac_cv_prog_cc_openmp=none])
       if test "$ac_cv_prog_cc_openmp" = no; then
         dnl Try these flags:
         dnl   GCC >= 4.2           -fopenmp
         dnl   SunPRO C             -xopenmp
         dnl   Intel C              -openmp
         dnl   SGI C, PGI C         -mp
         dnl   Tru64 Compaq C       -omp
         dnl   AIX IBM C            -qsmp=omp
         if test "$GCC" = yes; then
           dnl --- Test for GCC.
           gt_save_CFLAGS="$CFLAGS"
           CFLAGS="$CFLAGS -fopenmp"
           AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
             ])], [ac_cv_prog_cc_openmp="-fopenmp"])
           CFLAGS="$gt_save_CFLAGS"
         else
           dnl --- Test for SunPRO C.
           AC_EGREP_CPP([Brand], [
 #if defined __SUNPRO_C || defined __SUNPRO_CC
  Brand
 #endif
             ], ac_openmp_result=yes, ac_openmp_result=no)
           if test $ac_openmp_result = yes; then
             gt_save_CFLAGS="$CFLAGS"
             CFLAGS="$CFLAGS -xopenmp"
             AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
               ])], [ac_cv_prog_cc_openmp="-xopenmp"])
             CFLAGS="$gt_save_CFLAGS"
           else
             dnl --- Test for Intel C.
             AC_EGREP_CPP([Brand], [
 #if defined __INTEL_COMPILER
  Brand
 #endif
               ], ac_openmp_result=yes, ac_openmp_result=no)
             if test $ac_openmp_result = yes; then
               gt_save_CFLAGS="$CFLAGS"
               CFLAGS="$CFLAGS -openmp"
               AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
                 ])], [ac_cv_prog_cc_openmp="-openmp"])
               CFLAGS="$gt_save_CFLAGS"
             else
               dnl --- Test for SGI C, PGI C.
               AC_EGREP_CPP([Brand], [
 #if defined __sgi || defined __PGI || defined __PGIC__
  Brand
 #endif
                 ], ac_openmp_result=yes, ac_openmp_result=no)
               if test $ac_openmp_result = yes; then
                 gt_save_CFLAGS="$CFLAGS"
                 CFLAGS="$CFLAGS -mp"
                 AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
                   ])], [ac_cv_prog_cc_openmp="-mp"])
                 CFLAGS="$gt_save_CFLAGS"
               else
                 dnl --- Test for Compaq C.
                 AC_EGREP_CPP([Brand], [
 #if defined __DECC || defined __DECCXX
  Brand
 #endif
                   ], ac_openmp_result=yes, ac_openmp_result=no)
                 if test $ac_openmp_result = yes; then
                   gt_save_CFLAGS="$CFLAGS"
                   CFLAGS="$CFLAGS -omp"
                   AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
                     ])], [ac_cv_prog_cc_openmp="-omp"])
                   CFLAGS="$gt_save_CFLAGS"
                 else
                   dnl --- Test for AIX IBM C.
                   AC_EGREP_CPP([Brand], [
 #if defined _AIX
  Brand
 #endif
                     ], ac_openmp_result=yes, ac_openmp_result=no)
                   if test $ac_openmp_result = yes; then
                     gt_save_CFLAGS="$CFLAGS"
                     CFLAGS="$CFLAGS -qsmp=omp"
                     AC_COMPILE_IFELSE([AC_LANG_SOURCE([
 #ifndef _OPENMP
  Unlucky
 #endif
                       ])], [ac_cv_prog_cc_openmp="-qsmp=omp"])
                     CFLAGS="$gt_save_CFLAGS"
                   else
                     :
                   fi
                 fi
               fi
             fi
           fi
         fi
       fi
       ])
     case $ac_cv_prog_cc_openmp in
       none)
         AC_MSG_RESULT([none needed]) ;;
       no)
         AC_MSG_RESULT([unsupported]) ;;
       *)
         AC_MSG_RESULT([$ac_cv_prog_cc_openmp]) ;;
     esac
     case $ac_cv_prog_cc_openmp in
       none | no)
         OPENMP_CFLAGS= ;;
       *)
         OPENMP_CFLAGS=$ac_cv_prog_cc_openmp ;;
     esac
   fi
   AC_SUBST([OPENMP_CFLAGS])
 ])
