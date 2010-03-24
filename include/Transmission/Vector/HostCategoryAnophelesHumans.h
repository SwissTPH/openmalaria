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

#ifndef HOSTCATEGORYANOPHELESHUMANS_H_
#define HOSTCATEGORYANOPHELESHUMANS_H_

#include "HostCategoryAnopheles.h"
#include "util/WeibullDecayedValue.h"
#include "scenario.hxx"


/** Stores vector model data applicable between a category of host and a
 * mosquito species.
 *
 * This is a subclass of HostCategoryAnopheles which contains parameters
 * that are only relevant for Human hosts.
 *
 */

namespace OM { namespace Transmission {

class HostCategoryAnophelesHumans : public HostCategoryAnopheles {

public:
	HostCategoryAnophelesHumans() :
		    humanBloodIndex(0.0), probMosqOvipositing(0.0)
		  {}
	double humanBloodIndex;
	double probMosqOvipositing;
	/**
	 * this operator is used to set all the parameters
	 * for the human hosts.
	 *
	 */
	void operator= (const scnXml::Mosq& mosq);

	template<class S>
	  void operator& (S& stream) {
		humanBloodIndex & stream;
		probMosqOvipositing & stream;
	}

};
#endif /* HOSTCATEGORYANOPHELESHUMANS_H_ */
}}
