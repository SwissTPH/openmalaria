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


#ifndef TEXT_LINE_H
#define TEXT_LINE_H

#include "StringWrapper.h"
#include "Color.h"
#include "math_headers.h"
#include <list>

class Surface;
class SurfaceProvider;
class FontMM;

class Line
{
	public:

		Line(Surface* parent, int h);
		~Line();

		void clear();
		void print(String s);
		void print(int i);
		void setUsed(bool whether);
		void render(float2 charSize, float2 alignment);
		void render(float2 charSize, float2 alignment, Color top, Color bottom);
		void changeFont(FontMM* font);

	private:

		FontMM*							font;
		SurfaceProvider*	provider;
		Surface*					parent;
		int								textureYOffset, textureSize, cursor, maxCharWidth;
		float							widthInChars, texCoordWidth, texCoordHeight, texCoordOffsetY;
		unsigned int			texture;
};

typedef std::list<Line*> LineList;

#endif

