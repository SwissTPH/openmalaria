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


#ifndef TEXT_SURFACE_H
#define TEXT_SURFACE_H

#include <list>

#include "Line.h"

class SurfaceProvider;

class Surface
{
	friend class Line;

	public:

		Surface(	int size,		unsigned int texture, 
							FontMM* font, SurfaceProvider* provider	);

		~Surface();
    		
		unsigned int		charsPerLine, lineCount, linesInUse;		

		LineList&	getLines()
		{
			return lines;
		}

		unsigned int getTexture()
		{
			return texture;
		}

		void unload();


	private:

		SurfaceProvider*	provider;		
		FontMM*							font;
		int								size;
		unsigned int			texture;
		LineList					lines;
};

typedef std::list<Surface*> SurfaceList;

#endif

