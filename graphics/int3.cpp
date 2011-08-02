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


#include "int3.h"
#include <cmath>

int3::int3(int x, int y, int z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

int3::int3(double x, double y, double z)
:	x( (int) (x+0.5)),
	y( (int) (y+0.5)),
	z( (int) (z+0.5))
{
}

int3::int3(const int3& v)
{
	x = v.x;
	y = v.y;
	z = v.z;
}

int int3::length()
{
	return
		(int)
		sqrt(

			(double)
				(
					x*x
				+	y*y
				+	z*z
				)			
			);
}

int int3::lengthSquared()
{
	return		x*x
			+	y*y
			+	z*z ;
}
