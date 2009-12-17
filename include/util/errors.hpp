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
#ifndef OM_util_errors
#define OM_util_errors

#include <stdexcept>
using namespace std;

/** Standard exception classes for OpenMalaria. */
namespace OM { namespace util {
    
    /** Thrown to indicate an error in the scenario.xml file.  */
    class xml_scenario_error : public runtime_error
    {
	public:
	    explicit xml_scenario_error(const string&  __arg);
    };
    
    /** Thrown to indicate an error while loading/saving a checkpoint.
    *
    * Prepends "Error reading checkpoint: " to the message. */
    class checkpoint_error : public runtime_error
    {
	public:
	    explicit checkpoint_error(const string&  __arg);
    };
    
    /** Thrown to halt, when a command-line argument prompts an early exit.
    *
    * (Not an error; exits with status 0.) */
    class cmd_exit : public runtime_error
    {
	public:
	    explicit cmd_exit(const string& __arg);
    };
    
} }
#endif
