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


#ifndef COORDINATES3D_H
#define COORDINATES3D_H

#include <iostream>

class int3
{
	public:

		int3(int x, int y, int z = 0);
		int3(double x, double y, double z);
		int3(const int3& v);
		int3(){}

		int x, y, z;

		int length();
		int lengthSquared();
};

inline void operator += (int3& u, int3& v)
{
	u.x += v.x;
	u.y += v.y;
	u.z += v.z;
};

inline void operator -= (int3& u, int3& v)
{
	u.x -= v.x;
	u.y -= v.y;
	u.z -= v.z;
};

inline int3 operator + (int3& u, int3& v)
{
	return int3
		(	u.x + v.x,	
			u.y + v.y,
			u.z + v.z);
};

inline int3 operator - (int3& u, int3& v)
{
	return int3
		(	u.x - v.x,	
			u.y - v.y,
			u.z - v.z);
};

inline int3 operator * (double lambda, int3& v)
{
	return int3
		(	lambda * v.x,	
			lambda * v.y,
			lambda * v.z);
};

inline int3 operator * (int3& v, double lambda)
{
	return int3
		(	lambda * v.x,	
			lambda * v.y,
			lambda * v.z);
};

inline int3 operator / (int3& v, double lambda)
{
	return int3
		(	v.x / lambda,	
			v.y / lambda,
			v.z / lambda);
};

inline void operator *= (int3& v, double lambda)
{
	v.x = (int) (v.x * lambda + 0.5);
	v.y = (int) (v.y * lambda + 0.5);
	v.z = (int) (v.z * lambda + 0.5);	
};

inline void operator /= (int3& v, double lambda)
{
	v.x = (int) (v.x / lambda + 0.5);
	v.y = (int) (v.y / lambda + 0.5);
	v.z = (int) (v.z / lambda + 0.5);
};

inline int3 operator - (int3& v)
{
	return int3
		(	-v.x,	
			-v.y,
			-v.z);
};

inline std::ostream& operator << (std::ostream& s, int3& c)
{
	s << "(" << c.x << "," << c.y << ")";
	return s;
};

#endif
