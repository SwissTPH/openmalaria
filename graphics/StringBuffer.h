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


#ifndef STRING_BUFFER_H
#define STRING_BUFFER_H

#include <list>
#include "StringWrapper.h"
#include "macros.h"

typedef std::list<char> CharList;

class StringBuffer
{
	public:

		StringBuffer(){}
		StringBuffer(String s)
		{
			ITERATE(i, s.length())
			{ list.push_back(s[i]); }
		}

		std::list<char> list;

		void append(char c)
		{
			list.push_back(c);
		}

		void undo()
		{
			list.pop_back();
		}

		String getString()
		{
			unsigned int l = 0;
			char* c = new char[list.size() + 1];
			ITERATE_LIST(char, list, i, t, 
				c[l] = *i; l++;)
			c[list.size()] = 0;
			String out = String(c);
			delete c;
			return out;
		}
};

#endif
