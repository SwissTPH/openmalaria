/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "Transmission/Vector/HostCategoryAnopheles.h"

// NOTE: should Mosq be an extension of NonHumanHosts type or NonHumanHosts a sub-object of Mosq? This is duplicate functionality.
void HostCategoryAnopheles::operator= (const scnXml::Mosq& mosq) {
  entoAvailability = mosq.getMosqEntoAvailability();
  probMosqBiting = mosq.getMosqProbBiting();
  probMosqFindRestSite = mosq.getMosqProbFindRestSite();
  probMosqSurvivalResting = mosq.getMosqProbResting();
}

void HostCategoryAnopheles::operator= (const scnXml::NonHumanHosts& nnh) {
  entoAvailability = nnh.getMosqEntoAvailability();
  probMosqBiting = nnh.getMosqProbBiting();
  probMosqFindRestSite = nnh.getMosqProbFindRestSite();
  probMosqSurvivalResting = nnh.getMosqProbResting();
}
