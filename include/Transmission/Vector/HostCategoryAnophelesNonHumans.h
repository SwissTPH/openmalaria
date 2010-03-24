/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef HOSTCATEGORYANOPHELESNONHUMANS_H_
#define HOSTCATEGORYANOPHELESNONHUMANS_H_

#include "HostCategoryAnopheles.h"
#include "util/WeibullDecayedValue.h"
#include "scenario.hxx"

namespace OM { namespace Transmission {

/** Stores vector model data applicable between a category of host and a
 * mosquito species.
 *
 * This is a subclass of HostCategoryAnopheles which contains parameters
 * that are only relevant for Non Human hosts.
 *
 */

class HostCategoryAnophelesNonHumans : public HostCategoryAnopheles {

public:
	HostCategoryAnophelesNonHumans() :
	    relativeEntoAvailability(0.0)
	  {}

	/** non human host type name, we need this value to retrieve non human host population size */
	string nonHumanHostName;

	/** the relativeEntoAvailability is only useful if Non Human hosts types > 1, otherwise this value = 1 */
	double relativeEntoAvailability;

	/** this operator is used to set all the parameters for non human hosts */
	void operator= (const scnXml::NonHumanHosts& nnh);

	template<class S>
	  void operator& (S& stream) {
		nonHumanHostName & stream;
		relativeEntoAvailability & stream;
	}


};
typedef vector<HostCategoryAnophelesNonHumans> NonHumanHostsType;
}}
#endif /* HOSTCATEGORYANOPHELESNONHUMANS_H_ */
