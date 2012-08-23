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


#ifndef COLOR_H
#define COLOR_H

#include <iostream>
#include "gl_headers.h"

class Color
{
	public:

		//gcc doesn't like these inline defs inside /* ... */
		inline /*Color ::*/ Color()
		{}

		inline /*Color ::*/ Color(const Color& c)
		:		r(c.r), g(c.g), b(c.b), a(c.a)
		{}

		inline /*Color ::*/ Color(float r, float g, float b)
		:		r(r), g(g), b(b), a(1.0f)
		{}

		inline /*Color ::*/ Color(float r, float g, float b, float a)
		:		r(r), g(g), b(b), a(a)
		{}

		inline void set()
		{
			glColor4f(r,g,b,a);
		}

		inline void setOpaque()
		{
			glColor4f(r,g,b, 1.0f);
		}

		inline void setTransparent()
		{
			glColor4f(r,g,b, 0.0f);
		}

		inline void set(const Color& c)
		{
			r = c.r; g = c.g; b = c.b; a = c.a;
		}

		inline void set(float red, float green, float blue, float alpha)
		{
			r = red; g = green; b = blue; a = alpha;
		}

		inline void writeTo(float* array)
		{
			array[0] = r;
			array[1] = g;
			array[2] = b;
			array[3] = a;
		}

		float r, g, b, a;
};

inline void operator += (Color& u, const Color& v)
{
	u.r += v.r;
	u.g += v.g;
	u.b += v.b;
	u.a += v.a;
};

inline void operator -= (Color& u, const Color& v)
{
	u.r -= v.r;
	u.g -= v.g;
	u.b -= v.b;
	u.a -= v.a;
};

inline Color operator + (const Color& u, const Color& v)
{
	return Color
		(	u.r + v.r,	
			u.g + v.g,
			u.b + v.b,
			u.a + v.a);
};

inline Color operator - (const Color& u, const Color& v)
{
	return Color
		(	u.r - v.r,	
			u.g - v.g,
			u.b - v.b,
			u.a - v.a);
};

inline Color operator * (const float lambda, const Color& v)
{
	return Color
		(	lambda * v.r,	
			lambda * v.g,
			lambda * v.b,
			lambda * v.a);
};

inline Color operator * (const Color& v, const float lambda)
{
	return Color
		(	lambda * v.r,	
			lambda * v.g,
			lambda * v.b,
			lambda * v.a);
};

inline Color operator / (const Color& v, const float lambda)
{
	return Color
		(	v.r / lambda,	
			v.g / lambda,
			v.b / lambda,
			v.a / lambda);
};

inline void operator *= ( Color& v, const float lambda)
{
	v.r *= lambda;
	v.g *= lambda;
	v.b *= lambda;
	v.a *= lambda;
};

inline void operator /= ( Color& v, const float lambda)
{
	v.r /= lambda;
	v.g /= lambda;
	v.b /= lambda;
	v.a /= lambda;
}; 

inline Color operator * (const Color& v, const Color& w)
{
	return Color(v.r*w.r, v.g*w.g, v.b*w.b, v.a*w.a);
};

inline std::ostream& operator << (std::ostream& s, const Color& f)
{
	s << "(" << f.r << "," << f.g << "," << f.b << "," << f.a << ")";
	return s;
};

#endif
