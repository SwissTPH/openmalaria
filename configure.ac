#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.61)
AC_INIT([malariaModel],[0.9])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
#AC_CONFIG_HEADER([config.h])


# Checks for programs.

#check for CC and C++
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL

AC_PROG_RANLIB
AC_ARG_ENABLE(deploy,
 [  --enable-deploy  Use this for static linking against xerces, gsl and BOINC ],
 [ deploy=yes
 ])


# Checks for Darwin architecture (MacOS)
AM_CONDITIONAL(DARWIN, [test `uname -s | tr '[A-Z]' '[a-z]' | grep darwin`])

# needed paths for MacOS
if test `uname -s | tr '[A-Z]' '[a-z]' | grep darwin` ; then
  export LDFLAGS="-L/opt/local/lib -L/sw/lib -L../boinc/mac_build/build/Deployment -L../../boinc/mac_build/build/Deployment -isysroot /Developer/SDKs/MacOSX10.4u.sdk -Wl,-syslibroot,/Developer/SDKs/MacOSX10.4u.sdk -arch i386"
  export CPPFLAGS="-I/usr/local/include -I/opt/local/include -I/sw/include -I../boinc/lib -I../../boinc/lib -I../xerces-c-src_2_7_0/include/xercesc -isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386"
  export CFLAGS="-I/usr/local/include -I/opt/local/include -I/sw/include -I../boinc/lib -I../../boinc/lib -I../xerces-c-src_2_7_0/include/xercesc -isysroot /Developer/SDKs/MacOSX10.4u.sdk -arch i386"
  export SDKROOT="/Developer/SDKs/MacOSX10.4u.sdk"
  export MACOSX_DEPLOYMENT_TARGET=10.4
  export PATH=/opt/local/bin:$PATH
fi

#check libraries
if test "${deploy}" = yes ; then
    build_state='Deploy'
    AC_CHECK_FILE(lib/libxerces-c.a,[AC_SUBST(LIB_XERC, "../lib/libxerces-c.a")],[AC_MSG_FAILURE([No static xerces lib present])])
    AC_CHECK_FILE(lib/libgsl.a,[AC_SUBST(LIB_GSL, "../lib/libgsl.a")],[AC_MSG_FAILURE([No static gsl lib present])])
    AC_CHECK_FILE(lib/libgslcblas.a,[AC_SUBST(LIB_GSLCBLAS, "../lib/libgslcblas.a")],[AC_MSG_FAILURE([No static gslcblas lib present])])
    AC_CHECK_FILE(../boinc/api/libboinc_api.a,[AC_SUBST(LIB_BOINC_API, "../../boinc/api/libboinc_api.a")],[AC_MSG_FAILURE([No static boinc api library present])]) 
    AC_CHECK_FILE(../boinc/lib/libboinc.a,[AC_SUBST(LIB_BOINC, "../../boinc/api/libboinc.a")],[AC_MSG_FAILURE([No static boinc library present])]) 
    #check for boinc libs
    CHECK_BOINC
else
    build_state='local'
    AC_SUBST(LIB_XERC, "")
    AC_SUBST(LIB_GSL, "")
    AC_SUBST(LIB_GSLCBLAS,"")
    AC_SUBST(LIB_BOINC, "")
    AC_SUBST(LIB_BOINC_API, "")
     AC_DEFINE(WITHOUT_BOINC)
    AC_CHECK_LIB(gslcblas,main,,[AC_MSG_FAILURE([gslcblas not present])])
    AC_CHECK_LIB(gsl,main,, [AC_MSG_FAILURE([gsl not present])])
    if test `uname -s | grep Linux` ; then
        PKG_CHECK_MODULES(gsl, gsl >= 1.11,,[AC_MSG_FAILURE([Old gsl version. You must update your current version. Please follow the instructions in the wiki])])
    fi
    AC_CHECK_HEADERS([gsl/gsl_cdf.h gsl/gsl_multimin.h gsl/gsl_rng.h gsl/gsl_randist.h],, [AC_MSG_FAILURE([Some gsl headers are missing])])
    AC_CHECK_LIB(xerces-c,main,,[AC_MSG_FAILURE([xerces-c not present])])
fi



#check posix threads
AC_CHECK_LIB(pthread,main,, [AC_MSG_FAILURE([pthread not present])])

dnl Check for glut In the future we will check for the good glut
AX_CHECK_GL="-lGL"
AX_CHECK_GLUL="-lGLU"
AX_CHECK_GLUTL="-lglut"
if test "$no_gl" = yes -o "$no_glu" = yes -o "$no_glut" = yes; then
   have_glut=no
   AC_MSG_WARN([
================================================================================WARNING: openGL/GLU/glut not found.
The GL, GLU and glut libraries are required in order to build the graphical parts
of the BOINC application API library.
================================================================================]) 
else
   have_glut=yes
   echo "DEBUG: GLUT_CFLAGS = $GLUT_CFLAGS" >&5
   echo "DEBUG: GLUT_LIBS = $GLUT_LIBS" >&5

   AC_CHECK_HEADERS([gl.h glu.h glut.h glaux.h GL/gl.h GL/glu.h GL/glut.h GL/glaux.h OpenGL/gl.h OpenGL/glu.h OpenGL/glut.h OpenGL/glaux.h GLUT/glut.h MesaGL/gl.h MesaGL/glu.h MesaGL/glut.h MesaGL/glaux.h])
fi

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h],,[AC_MSG_FAILURE([Some headers are missing])])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T


# Checks for library functions.
AC_FUNC_ERROR_AT_LINE

AC_CONFIG_FILES([Makefile
                 model/Makefile
                 xsdcxx/Makefile
                 graphics/Makefile])

AC_OUTPUT
