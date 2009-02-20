# Modified from seti_boinc: http://setiathome.berkeley.edu/sah_porting.php

AC_DEFUN([CHECK_BOINC],[
  AC_ARG_VAR([BOINCDIR],[boinc directory])
  if test -z "$HEAD"
  then
    AC_PATH_PROG(HEAD,head)
  fi
  if test -z "$FIND"
  then
    AC_PATH_PROG(FIND,find)
  fi
  thisdir=`pwd`
  AC_MSG_CHECKING([for BOINC])
  boinc_search_path="$BOINCDIR boinc ../boinc $HOME/boinc "
  for boinc_dir in $boinc_search_path
  do
    if test -d $boinc_dir 
    then
      if test -f $boinc_dir/include/std_fixes.h -o -f $boinc_dir/lib/std_fixes.h
      then
        cd $boinc_dir
        BOINCDIR=`pwd`
	cd $thisdir
	break
      else
        if $FIND $boinc_dir -name "std_fixes.h" >& /dev/null
	then
	  BOINCDIR=`$FIND $boinc_dir -name "std_fixes.h" -print | $HEAD -1 | sed 's/\/std_fixes.h//'`         
          cd $BOINCDIR/..
          BOINCDIR=`pwd`
	  cd $thisdir
	  break
	fi
      fi
    fi
  done
  if test -n "$BOINCDIR" 
  then
    AC_MSG_RESULT($BOINCDIR)
  else
    AC_MSG_WARN([ boinc not found.
============================================================================
boinc not found. You can get it here: http://boinc.ssl.berkeley.edu/
============================================================================
])
    exit 1
  fi
  AC_SUBST([BOINCDIR])
  BOINC_CFLAGS="-I$BOINCDIR -I$BOINCDIR/api -I$BOINCDIR/lib"
  AC_SUBST([BOINC_CFLAGS])
  BOINC_LDFLAGS="-L$BOINCDIR -L$BOINCDIR/api -L$BOINCDIR/lib"
  AC_SUBST([BOINC_LDFLAGS])
  LIB_BOINC="-lboinc_api -lboinc"
  AC_SUBST([LIB_BOINC])  
])
	

