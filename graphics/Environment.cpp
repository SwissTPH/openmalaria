/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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


#include "Environment.h"

Environment Environment :: interpolate(
									const Environment & a, const Environment & b, 
									float dt)
{
	float at = 1.0f - dt;
	return Environment ( 
					at * a.sunlight + dt * b.sunlight,
					at * a.sky + dt * b.sky,
					at * a.sun + dt * b.sun,
					at * a.ambient + dt * b.ambient,
					at * a.shadow + dt * b.shadow,
					at * a.haze + dt * b.haze );
}

Environment :: Environment(	const Color & sunlight,
														const Color & sky,
														const Color & sun,
														const Color & ambient,
														const Color & shadow,
														const Color & haze )
:		sunlight(sunlight), sky(sky), sun(sun), 
		ambient(ambient), shadow(shadow), haze(haze)
{}

Environment :: Environment(	const Environment& e )
:		sunlight(e.sunlight), sky(e.sky), sun(e.sun), 
		ambient(e.ambient), shadow(e.shadow), haze(e.haze)
{}

