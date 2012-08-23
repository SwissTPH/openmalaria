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


#ifndef COORDINATES_H
#define COORDINATES_H

#include <iostream>

#pragma	warning( disable : 4244 )

class int2
{
	public:

		int2(int x, int y);
		int2(const int2& v);
		int2(){}

		int x, y;

		int length();
		int lengthSquared();
};

inline bool operator == (const int2& u, const int2& v)
{
	return (u.x == v.x)&&(u.y == v.y);
};

inline void operator += (int2& u, int2& v)
{
	u.x += v.x;
	u.y += v.y;	
};

inline void operator -= (int2& u, int2& v)
{
	u.x -= v.x;
	u.y -= v.y;
};

inline int2 operator + (int2& u, int2& v)
{
	return int2
		(	u.x + v.x,	
			u.y + v.y);
};

inline int2 operator - (int2& u, int2& v)
{
	return int2
		(	u.x - v.x,	
			u.y - v.y);
};

inline int2 operator * (double lambda, int2& v)
{
	return int2
		(	(int)lambda * v.x,	
			(int)lambda * v.y);
};

inline int2 operator * (int2& v, double lambda)
{
	return int2
		(	(int)lambda * v.x,	
			(int)lambda * v.y);
};

inline int2 operator / (int2& v, double lambda)
{
	return int2
		(	(int)(v.x / lambda),	
			(int)(v.y / lambda));
};

inline void operator *= (int2& v, double lambda)
{
	v.x = (int) (v.x * lambda + 0.5);
	v.y = (int) (v.y * lambda + 0.5);	
};

inline void operator /= (int2& v, double lambda)
{
	v.x = (int) (v.x / lambda + 0.5);
	v.y = (int) (v.y / lambda + 0.5);
};

inline int2 operator - (int2& v)
{
	return int2
		(	-v.x,	
			-v.y);
};

inline std::ostream& operator << (std::ostream& s, int2& c)
{
	s << "(" << c.x << "," << c.y << ")";
	return s;
};

#endif
