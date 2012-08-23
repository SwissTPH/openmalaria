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


#ifndef FONT_H
#define FONT_H

#include "StringWrapper.h"
#include "int2.h"

class FontMM
{
	public:

		FontMM(	String filename, int2 charSize, int2 tileSize, int2 tileOffset);
		~FontMM();

		unsigned char**		data;
		int2							charSize, tileSize, tileOffset;

	private:

		void	extractCharacters(unsigned char* image, int w, int h);
};

#endif

