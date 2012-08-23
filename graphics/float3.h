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


#ifndef FLOAT_VECTOR_H
#define FLOAT_VECTOR_H

#include <iostream>
#include <cmath>
#include "double3.h"

class float4;

class float3
{
	public:

		float3(	const float x, 
				const float y, 
				const float z);

		float3(const float3& v);

		float3(const double3& v);

		float3(){}

		float x, y, z;		

		inline float length()
		{
			return
				(float) sqrt(
								x*x
							+	y*y
							+	z*z
									);
		}

		inline float lengthSquared()
		{
			return x*x + y*y + z*z;
		}

		inline double lengthSquaredDouble()
		{
			double x = (double)this->x;
			double y = (double)this->y;
			double z = (double)this->z;

			return x*x + y*y + z*z;
		}

		inline float3 direction()
		{
			float l = (float) 
				sqrt(	x*x + y*y + z*z );
			return float3( x/l, y/l, z/l );
		}

		inline float3 
			cross( const float3& v)
		{
			return float3
				(	y*v.z - z*v.y,
					z*v.x - x*v.z,
					x*v.y - y*v.x	);
		}

		float3& operator = (const float3 & b)
		{
			x = b.x;
			y = b.y;
			z = b.z;
			return (*this);
		}

		inline void writeInto(float* array)
		{
			array[0] = x; array[1] = y;
			array[2] = z;
		}
};

inline void operator += (float3& u, const float3& v)
{
	u.x += v.x;
	u.y += v.y;
	u.z += v.z;
};

inline void operator -= (float3& u, const float3& v)
{
	u.x -= v.x;
	u.y -= v.y;
	u.z -= v.z;
};

inline float3 operator + (const float3& u, const float3& v)
{
	return float3
		(	u.x + v.x,	
			u.y + v.y,
			u.z + v.z);
};

inline float3 operator - (const float3& u, const float3& v)
{
	return float3
		(	u.x - v.x,	
			u.y - v.y,
			u.z - v.z);
};

inline float3 operator * (const float lambda, const float3& v)
{
	return float3
		(	lambda * v.x,	
			lambda * v.y,
			lambda * v.z);
};

inline float3 operator * (const float3& v, const float lambda)
{
	return float3
		(	lambda * v.x,	
			lambda * v.y,
			lambda * v.z);
};

inline float3 operator / (const float3& v, const float lambda)
{
	return float3
		(	v.x / lambda,	
			v.y / lambda,
			v.z / lambda);
};

inline void operator *= ( float3& v, const float lambda)
{
	v.x *= lambda;
	v.y *= lambda;
	v.z *= lambda;
};

inline void operator /= ( float3& v, const float lambda)
{
	v.x /= lambda;
	v.y /= lambda;
	v.z /= lambda;
};

inline float3 operator - (const float3& v)
{
	return float3
		(	-v.x,	
			-v.y,
			-v.z);
};

inline float operator * (const float3& v, const float3& w)
{
	return v.x*w.x + v.y*w.y + v.z*w.z;
};

inline float3 operator || (const float3& v, const float3& w)
{
	float length_w = (float) sqrt( w.x*w.x	+ w.y*w.y + w.z*w.z );
	float length_v = (float) sqrt( v.x*v.x	+ v.y*v.y + v.z*v.z );

	if (length_v*length_w == 0.0f) return float3(0.0f,0.0f,0.0f);

	return		v *	(v.x*w.x + v.y*w.y + v.z*w.z) 
					/	(length_v * length_w);			
}

inline std::ostream& operator << (std::ostream& s, float3& f)
{
	s << "(" << f.x << "," << f.y << "," << f.z << ")";
	return s;
};

#endif

