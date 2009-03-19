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

/**@brief A wrapper around BOINC. This header should not include any BOINC headers!
 * 
 * The purpose of this is to avoid most other source files from having to
 * include BOINC headers, since these have caused some wierd issues with the
 * MSVC++ compiler. */

#ifndef Hmod_boincWrapper
#define Hmod_boincWrapper

#include <string>

int boincWrapper_fraction_done(double progress);

int boincWrapper_time_to_checkpoint();
int boincWrapper_checkpoint_completed(void);

int boincWrapper_resolve_filename(const char* name, char* outName, int len);

int boincWrapper_resolve_filename_s(const char* name, std::string& outName);

int boincWrapper_finish(int status);

#endif
