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

// A wrapper around BOINC. This header should not include any BOINC headers!

#include "boincWrapper.h"
#include "boinc_bridge.h"

int boincWrapper_fraction_done(double progress) {
  return boinc_fraction_done (progress);
}

int boincWrapper_time_to_checkpoint() {
  return boinc_time_to_checkpoint();
}
int boincWrapper_checkpoint_completed(void) {
  return boinc_checkpoint_completed();
}

int boincWrapper_resolve_filename(const char* name, char* outName, int len) {
  return boinc_resolve_filename (name, outName, len);
}

int boincWrapper_resolve_filename_s(const char* name, std::string& outName) {
  return boinc_resolve_filename_s (name, outName);
}

int boincWrapper_finish(int status) {
  return boinc_finish (status);
}
