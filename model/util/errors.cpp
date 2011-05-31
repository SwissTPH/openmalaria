/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "util/errors.h"
#include <ostream>
#include <iostream>

#ifdef __GNU_LIBRARY__
#include <execinfo.h>
#include <stdlib.h>
#include <cstring>
#include <cxxabi.h>
#endif

/** Standard exception classes for OpenMalaria. */
namespace OM { namespace util {

traced_exception::traced_exception(const string& msg)
 : runtime_error(msg) {
#ifdef __GNU_LIBRARY__
    // http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
// probably 100 is big enough
#define MAX_STACK_SIZE 100
    void *array[MAX_STACK_SIZE];
    
    length = backtrace (array, MAX_STACK_SIZE);
    trace = backtrace_symbols (array, length);
    
    //TODO: work out why OM catches this as an unknown exception
    //until then, dump the stack trace here
    cerr << "Exception: "<<msg<<'\n';
    cerr << *this << endl;
#endif
}
traced_exception::~traced_exception() throw(){
#ifdef __GNU_LIBRARY__
    free (trace);
#endif
}
ostream& operator<<(ostream& stream, const traced_exception& e){
#ifdef __GNU_LIBRARY__
    // demangle output: http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
    // spec: http://www.ib.cnea.gov.ar/~oop/biblio/libstdc++/namespaceabi.html
    
    // start from 1, since the first entry is the traced_exception ctor
    if( e.length >= 1 ){
        char *line = e.trace[0];
        char *lb = strchr(line, '(');
        if( lb != 0 ){
            *lb = '\0';
            stream << "in " << line << ":\n";
        }
    }
    for( size_t i=1; i<e.length; ++i ){
        char *line = e.trace[i];
        char *lb = strchr(line, '(');
        char *plus = 0;
        char *rb = 0;
        char *demangled = 0;
        if( lb != 0 ){
            *lb = '\0';
            plus = strchr(lb+1,'+');
            if( plus != 0 ){
                *plus = '\0';
                int status;
                demangled = abi::__cxa_demangle(lb+1, 0, 0, &status);
                /*if( status != 0 ){
                    cerr << "status: "<<status<<"\tstring:"<<lb+1<<endl;
                }*/
                rb = strchr(plus+1, ')');
            }else{
                rb = strchr(lb+1, ')');
            }
            if( rb != 0 ){
                *rb = '\0';
            }
        }
        //stream << line << ": ";
        if( plus != 0 ){
            stream << "+" << (plus+1);
        }
        if( demangled != 0 ){
            stream << '\t' << demangled;
            free(demangled);
        }else if( lb != 0 && rb > lb ){
            stream << '\t' << (lb+1) << "()";
        }
        stream << '\n';
    }
#else
    stream << "sorry, no trace\n";
#endif
    return stream;
}
 
xml_scenario_error::xml_scenario_error(const string&  msg)
  : runtime_error(msg) { }

checkpoint_error::checkpoint_error(const string&  msg)
  : traced_exception(msg) { }

cmd_exit::cmd_exit(const string& msg)
  : runtime_error(msg) { }
   
} }
