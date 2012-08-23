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


#ifndef KEY_H
#define KEY_H

typedef enum 
{	
	L_SHIFT, L_CTRL, L_ALT, R_SHIFT, R_CTRL, R_ALT, 
	BACKSPACE, DELETE_KEY, RETURN_KEY, SPACE,
	L_CURSOR, R_CURSOR, U_CURSOR, D_CURSOR,
	INSERT_KEY, HOME, END, PG_UP, PG_DN, ESCAPE,
	F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12
} 
	SpecialKey;

class Key
{
	public:

		Key() {}

		Key(const Key& c) 
		{
			isCharacter = c.isCharacter;
			character		= c.character;
			specialKey	= c.specialKey;
		}
		
		bool					isCharacter;
		int						character;
		SpecialKey		specialKey;
};

inline bool operator < (const Key & a, const Key & b)
{
	return a.character < b.character;
};

#endif

