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


#ifndef FLOATHEXAGON_H
#define FLOATHEXAGON_H

#include "float3.h"
#include "float4.h"
#include "plane.h"
#include <iostream>

class float6
{
	public:

		float x0,y0,z0,x1,y1,z1;

		float6()		
		{}

		float6(	float x0, float y0, float z0, 
				float x1, float y1, float z1 )
		:	x0(x0), y0(y0), z0(z0),
			x1(x1), y1(y1), z1(z1)
		{}

		float6(	double x0, double y0, double z0, 
		double x1, double y1, double z1 )
		:	x0((float)x0), y0((float)y0), z0((float)z0),
			x1((float)x1), y1((float)y1), z1((float)z1)
		{}

		float6(float3 begin, float3 end)
		:	x0(begin.x), y0(begin.y), z0(begin.z),
			x1(end.x), y1(end.y), z1(end.z)
		{}

		float6(const float6& q)
		:	x0(q.x0), y0(q.y0), z0(q.z0), 
			x1(q.x1), y1(q.y1), z1(q.z1)
		{}

		void set(const float6& q)
		{
			x0 = q.x0; y0 = q.y0; z0 = q.z0; 
			x1 = q.x1; y1 = q.y1; z1 = q.z1;
		}

		inline float3 intersectPlane(Plane plane)
		{	 
			float3 normal = plane.normal;
			float directionsWorthTowardsPlane 
				= -(x1*normal.x + y1*normal.y + z1*normal.z);
			float originsDistanceFromPlane
				= x0*normal.x + y0*normal.y + z0*normal.z
					- plane.offset;
			return float3(x0,y0,z0)
				+(originsDistanceFromPlane/directionsWorthTowardsPlane)
					*float3(x1,y1,z1);
		}

		inline float6 operator = (const float6& v)
		{
			return float6(	x0 = v.x0,
							y0 = v.y0,
							z0 = v.z0,
							x1 = v.x1,
							y1 = v.y1,
							z1 = v.z1 );
		};
};

inline std::ostream& operator << (std::ostream& s, float6& c)
{
	s	<< "( " << c.x0 << "," << c.y0 << "," << c.z0 << ", " 
		<< c.x1 << "," << c.y1 << "," << c.z1 <<" )";
	return s;
};



#endif
