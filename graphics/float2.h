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


#ifndef FLOAT_VECTOR_2D_H
#define FLOAT_VECTOR_2D_H

#include <cmath>
#include "int2.h"

class float2
{
	public:

		float2(	const float x, 
				const float y)
				: x(x), y(y) {}

		float2(	const double x, 
				const double y)
				: x((float)x), y((float)y) {}

		float2(const float2& v)
			: x(v.x), y(v.y) {}
								
		float2(){}

		float x, y;		

		inline float length()
		{
			return (float) sqrt(x*x + y*y);
		}

		float lengthSquared()
		{
			return x*x + y*y;
		}

		inline float2 direction()
		{
			float l = length();
			return float2( x/l, y/l );
		}	

		inline float2 scaledBy(const float2 & c)
		{
			return float2(x*c.x, y*c.y);
		}

		inline void scaleBy(const float2 & c)
		{
			x *= c.x;		y *= c.y;
		}

		inline void operator = (const float2& v)
		{
			x = v.x;
			y = v.y;
		};

		inline void set(const float2& v)
		{
			x = v.x;
			y = v.y;
		};
};

inline bool operator == (const float2& u, const float2& v)
{
	return (u.x == v.x)&&(u.y == v.y);
};

inline void operator += (float2& u, const float2& v)
{
	u.x += v.x;
	u.y += v.y;
};

inline void operator -= (float2& u, const float2& v)
{
	u.x -= v.x;
	u.y -= v.y;	
};

inline float2 operator + (const float2& u, const float2& v)
{
	return float2
		(	u.x + v.x,	
			u.y + v.y	);
};

inline float2 operator - (const float2& u, const float2& v)
{
	return float2
		(	u.x - v.x,	
			u.y - v.y	);
};

inline float2 operator * (const float lambda, const float2& v)
{
	return float2
		(	lambda * v.x,	
			lambda * v.y	);
};

inline float2 operator * (const float2& v, const float lambda)
{
	return float2
		(	lambda * v.x,	
			lambda * v.y	);
};

inline float2 operator / (const float2& v, const float lambda)
{
	return float2
		(	v.x / lambda,	
			v.y / lambda	);
};

inline float2 operator / ( float2& u, float2& v)
{
	if (v.x*v.y == 0.0f) return float2(0.0f,0.0f);
	return float2
		(	u.x / v.x,	
			u.y / v.y	);
};

inline void operator *= ( float2& v, const float lambda)
{
	v.x *= lambda;
	v.y *= lambda;
};

inline void operator /= ( float2& u, float2& v)
{
	if (v.x*v.y == 0.0f) return;
	u.x /= v.x;
	u.y /= v.y;
};

inline void operator /= ( float2& v, const float lambda)
{
	v.x /= lambda;
	v.y /= lambda;
};

inline float2 operator - (const float2& v)
{
	return float2
		(	-v.x,	
			-v.y	);
};

inline int2 operator * (const float2& lambda, const int2& c)
{
	return int2
		(	(int) (c.x * lambda.x),
			(int) (c.y * lambda.y)	);
};

inline int2 operator * (const int2& c, const float2& lambda)
{
	return int2
		(	(int) (c.x * lambda.x),
			(int) (c.y * lambda.y)	);
};

inline std::ostream& operator << (std::ostream& s, float2& c)
{
	s << "(" << c.x << "," << c.y << ")";
	return s;
};

#endif

