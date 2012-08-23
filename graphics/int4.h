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


#ifndef INTSQUARE_H
#define INTSQUARE_H

#include "int2.h"
#include "int3.h"
#include <iostream>

class int4
{
	public: 

		int x0,y0,x1,y1;

		int4()		
		{}

		int4(int x0, int y0, int x1, int y1)
		:	x0(x0), y0(y0), x1(x1), y1(y1)
		{}

		int4(int2 begin, int2 end)
		:	x0(begin.x), y0(begin.y), x1(end.x), y1(end.y)
		{}

		int4(const int4& q)
		:	x0(q.x0), y0(q.y0), 
			x1(q.x1), y1(q.y1)
		{}

		void set(const int4& q)
		{
			x0 = q.x0; y0 = q.y0; 
			x1 = q.x1; y1 = q.y1;
		}

		inline int4 operator = (const int4& v)
		{
			return int4(	x0 = v.x0,
							y0 = v.y0,
							x1 = v.x1,
							y1 = v.y1 );
		};
};

inline std::ostream& operator << (std::ostream& s, int4& c)
{
	s << "( " << c.x0 << ", " << c.y0 << ", " << c.x1 << ", " << c.y1 <<" )";
	return s;
};

#endif
