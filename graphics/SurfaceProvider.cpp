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


#include "SurfaceProvider.h"
#include "Surface.h"
#include "Font.h"
#include "gl_headers.h"
#include "macros.h"

SurfaceProvider*  SurfaceProvider :: instance = 0;

//	c'tor
SurfaceProvider :: SurfaceProvider(int baseSize, FontMM* font)
:	size(baseSize),
	defaultFont(font)
{}

//	getInstance
SurfaceProvider* SurfaceProvider :: getInstance()
{
	return instance;
}

//		init
void SurfaceProvider :: init(int baseSize, FontMM* font)
{
	instance = new SurfaceProvider(baseSize, font);
}

//		getLine
Line*		SurfaceProvider :: getLine()
{ 	
	if (freeLines.empty()) 
	{ 	
		Surface* s = createSurface();
		allSurfaces.insert(allSurfaces.end(),s);		
		LineList l = s->getLines();
		freeLines.insert(freeLines.end(), l.begin(), l.end());
	}
	Line* l = *freeLines.begin();
	freeLines.remove(l);
	l->setUsed(true);

	return l;
}

//		createSurface
Surface*	SurfaceProvider :: createSurface(FontMM* font)
{
	if (!font) font = defaultFont;

	GLuint index;
	glGenTextures(1, &index);
	glBindTexture(GL_TEXTURE_2D, index);
	loadBlank();
	textures.insert(textures.end(),index);

	Surface* s = new Surface(size, index, font, this);
	return s;
}

//		loadBlank
void SurfaceProvider :: loadBlank()
{
	int s = size*size*4;
	unsigned char* blank = new unsigned char[s];

	for (int i = 0; i < s; i++)
		blank[i] = 0;
  	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	glTexImage2D(	GL_TEXTURE_2D, 0, GL_RGBA, size, size, 0, 
								GL_RGBA, GL_UNSIGNED_BYTE, blank);

	delete blank;
}

//		remove
void SurfaceProvider :: remove(Surface* s, unsigned int texture)
{
	allSurfaces.remove(s);
	LineList& l = s->getLines();
	ITERATE_LIST(Line*,l,i,t, freeLines.remove(*i);)
	glDeleteTextures(1, &texture);
}
