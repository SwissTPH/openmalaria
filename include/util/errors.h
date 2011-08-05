/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef OM_util_errors
#define OM_util_errors

#include <stdexcept>
#include <boost/static_assert.hpp>

using namespace std;

/* Macros to ease use of traced_exception.
 * msg should be obvious, and code is the program exit code.
 */
#define TRACED_EXCEPTION( msg, code ) OM::util::traced_exception( (msg), __FILE__, __LINE__, (code) );
#define TRACED_EXCEPTION_DEFAULT( msg ) OM::util::traced_exception( (msg), __FILE__, __LINE__ );

/** Standard exception classes for OpenMalaria. */
namespace OM { namespace util {
    
    /** Different exit codes which may be used.
     * 
     * These are intended to categorise errors, not diagnose them. Always use
     * error messages or stack traces to diagnose errors. (The only reason we
     * don't just return 1 is that BOINC uses the codes to correlate errors.)
     * 
     * We start exit codes at 64, as in /usr/include/sysexits.h; as far as I am
     * aware, linux is the only platform targetted which has such narrow ranges
     * of exit codes available (hence these codes should work on Windows and
     * Mac OS X too). */
    namespace Error { enum ErrorCodes {
        /// no error (exit code 0)
        None = 0,
        Default = 64,
        TracedDefault,
        /// any error from Code Synthesis's XSD
        XSD,
        /// any checkpointing error (these don't occur often enough to need segregating here)
        Checkpoint,
        /// invalid scenario file
        XmlScenario,
        GSL,
        FileExists,		// wanted to create a file but it already exists!
        FileIO,			// any other file read/write error
        EffectiveEIR,
        NumNewInfections,
        InitialKappa,
        VectorWarmup,
        Checksum,
        CommandLine,
        SumWeight,
        VectorFitting,
        InfLambda,
        Max
    }; }
    
    // As in the "Advanced Bash-Scripting Guide"; not directly relevant to C++
    // but gives some idea what codes make sense to use.
    BOOST_STATIC_ASSERT( Error::Max <= 113 );
    
    /** Base OpenMalaria exception class. */
    class base_exception : public runtime_error {
    public:
        // for common errors, specify a unique TRACED_DEFAULT to help segregate
        // errors in BOINC server summaries; otherwise leave the default
        explicit base_exception(const string& msg, int code=Error::Default);
        inline int getCode() const{
            return errCode;
        }
    private:
        int errCode;
    };
    
    /** Extension of runtime_error which tries to get a stack trace. */
    class traced_exception : public base_exception {
    public:
        /** Create a stack-trace and store in this object.
         * 
         * @param msg Error message
         * @param file File in which it occurred (may pass NULL)
         * @param line Line of error (may pass 0)
         * @param start Index of first stack frame of interest. 0 is this
         * constructor, 1 the code creating this exception (usually the location
         * of the problem, but sometimes an error handler), 2 the next frame etc.
         * @param code Exit code.
         * 
         * See TRACED_EXCEPTION and TRACED_EXCEPTION_DEFAULT macros which
         * reduce the number of arguments which need to be passed.
         */
        explicit traced_exception(const string& msg, const char *file, int line, int code=Error::TracedDefault, int start=1);
        virtual ~traced_exception() throw();
    private:
        const char *file;
        int line, start;
#ifdef __GNU_LIBRARY__
        size_t length;
        char** trace;
#endif
        friend ostream& operator<<(ostream& stream, const traced_exception& e);
    };
    /// Print trace to stream (a series of lines, each with \\n appended)
    ostream& operator<<(ostream& stream, const traced_exception& e);
    
    /** Thrown to indicate an error in the scenario.xml file.  */
    class xml_scenario_error : public base_exception {
    public:
        explicit xml_scenario_error(const string& msg);
    };
    
    /** Thrown to indicate an error while loading/saving a checkpoint.
    *
    * Prepends "Error reading checkpoint: " to the message. */
    class checkpoint_error : public traced_exception {
    public:
        explicit checkpoint_error(const string& msg);
    };
    
    /** Thrown due to command-line error or to halt, when a command-line
     * argument prompts an early exit. Not an error when errCode==0. */
    class cmd_exception : public base_exception {
    public:
        explicit cmd_exception(const string& msg, int code=Error::CommandLine);
    };
    
    /** Set up a GSL handler to throw a traced exception. */
    void set_gsl_handler();
    
} }
#endif
