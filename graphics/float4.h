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


#ifndef FLOATSQUARE_H
#define FLOATSQUARE_H

#include "float2.h"
#include "float3.h"
#include <iostream>

class float4
{
	public: 

		float x0,y0,x1,y1;

		float4()		
		{}

		float4(float x0, float y0, float x1, float y1)
		:	x0(x0), y0(y0), x1(x1), y1(y1)
		{}

		float4(float2 begin, float2 end)
		:	x0(begin.x), y0(begin.y), x1(end.x), y1(end.y)
		{}

		float4(const float4& q)
		:	x0(q.x0), y0(q.y0), 
			x1(q.x1), y1(q.y1)
		{}

		void set(const float4& q)
		{
			x0 = q.x0; y0 = q.y0; 
			x1 = q.x1; y1 = q.y1;
		}

		inline bool planeContainsPoint(float3 point)
		{	return (point.x*x0 + point.y*y0 + point.z*x1 == y1); }

		inline bool planeContainsVector(float3 point)
		{	return (point.x*x0 + point.y*y0 + point.z*x1 == 0); }

		inline float4 operator = (const float4& v)
		{
			return float4(	x0 = v.x0,
							y0 = v.y0,
							x1 = v.x1,
							y1 = v.y1 );
		};
};

inline std::ostream& operator << (std::ostream& s, float4& c)
{
	s << "( " << c.x0 << ", " << c.y0 << ", " << c.x1 << ", " << c.y1 <<" )";
	return s;
};



#endif

