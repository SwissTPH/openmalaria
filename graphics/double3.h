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


#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>

class float3;

class double3
{
	public:
	
		double x, y, z;	
		
		double3(	const double x, 
					const double y, 
					const double z);
				
		double3(const double3& v);
		double3(const float3& v);
		double3(){}	

		inline double length()
		{
			return
				sqrt(
						x*x
					+	y*y
					+	z*z
								);
		}

		inline double lengthSquared()
		{
			return	x*x
				+	y*y
				+	z*z ;
		}

		inline double3 direction()
		{
			double l = length();
			return 
				double3( x/l, y/l, z/l );
		}

		inline double3 
			cross( const double3& v)
		{
			return double3
				(	y*v.z - z*v.y,
					z*v.x - x*v.z,
					x*v.y - y*v.x	);
		}
};

inline void operator += (double3& u, const double3& v)
{
	u.x += v.x;
	u.y += v.y;
	u.z += v.z;
};

inline void operator -= (double3& u, const double3& v)
{
	u.x -= v.x;
	u.y -= v.y;
	u.z -= v.z;
};

inline double3 operator + (const double3& u, const double3& v)
{
	return double3
		(	u.x + v.x,	
			u.y + v.y,
			u.z + v.z);
};

inline double3 operator - (const double3& u, const double3& v)
{
	return double3
		(	u.x - v.x,	
			u.y - v.y,
			u.z - v.z);
};

inline double3 operator * (const double lambda, const double3& v)
{
	return double3
		(	lambda * v.x,	
			lambda * v.y,
			lambda * v.z);
};

inline double3 operator * (const double3& v, const double lambda)
{
	return double3
		(	lambda * v.x,	
			lambda * v.y,
			lambda * v.z);
};

inline double3 operator / (const double3& v, const double lambda)
{
	return double3
		(	v.x / lambda,	
			v.y / lambda,
			v.z / lambda);
};

inline void operator *= ( double3& v, const double lambda)
{
	v.x *= lambda;
	v.y *= lambda;
	v.z *= lambda;
};

inline void operator /= ( double3& v, const double lambda)
{
	v.x /= lambda;
	v.y /= lambda;
	v.z /= lambda;
};

inline double3 operator - (const double3& v)
{
	return double3
		(	-v.x,	
			-v.y,
			-v.z);
}; 

inline double operator * (const double3& v, const double3& w)
{
	return v.x*w.x + v.y*w.y + v.z*w.z;
};

inline double3 operator || (const double3& v, const double3& w)
{
	double length_w = sqrt( w.x*w.x	+ w.y*w.y + w.z*w.z );
	double length_v = sqrt( v.x*v.x	+ v.y*v.y + v.z*v.z );

	if (length_v*length_w == 0.0) return double3(0.0f,0.0f,0.0f);

	return		v *	(v.x*w.x + v.y*w.y + v.z*w.z) 
					/	(length_v * length_w);			
};

inline std::ostream& operator << (std::ostream& s, double3& f)
{
	s << "(" << f.x << "," << f.y << "," << f.z << ")";
	return s;
};

#endif

