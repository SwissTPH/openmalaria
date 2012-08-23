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


#include "Font.h"
#include "IL/il.h"
#include "IL/ilu.h"
#include "macros.h"
#include "gl_headers.h"
#include "Color.h"
#include <iostream>
#include <cmath>

FontMM :: FontMM(String filename, int2 charSize, int2 tileSize, int2 tileOffset)
:	data(new unsigned char*[256]),
	charSize(charSize),
	tileSize(tileSize),
	tileOffset(tileOffset)
{
	LOG("loading '" << filename << "'")
	ILuint image;

  ilGenImages(1, &image);
	ilBindImage(image);

 	if (!ilLoadImage(ILstring(filename.c_str())))	
	{ CRASH( "unable to load " << filename << "!" )	} 

	extractCharacters (
		ilGetData(), 
		ilGetInteger(IL_IMAGE_WIDTH), 
		ilGetInteger(IL_IMAGE_HEIGHT)
	);

	ilDeleteImages(1, &image);
}

FontMM :: ~FontMM()
{
	for (int i = 0; i < 256; i++)
		delete data[i];
	delete data;
}

void FontMM :: extractCharacters(unsigned char* image, int w, int h)
{ 
	int		x, y, x0, y0, imgIndex, charIndex, 
				s = charSize.x*charSize.y, lpr = w / tileSize.x;

	unsigned char* character;

	for (int i = 0; i < 256; i++)
	{
		character = new unsigned char[4*s];
		x0 = tileSize.x*(i%lpr) + tileOffset.x;
		y0 = tileSize.y*(i/lpr) + tileOffset.y;
		for (int r = 0; r < charSize.y; r++)
		for (int c = 0; c < charSize.x; c++)
		{
			x = x0 + c; y = y0 + r;
			imgIndex = 4*(y*w + x);
			charIndex = 4*(r*charSize.x + c);
			if (imgIndex < 0 || imgIndex >= 4*w*h)
			{	
				for (int j = 0; j < 4; j++) character[charIndex + j] = 0;
			}
			else
			{
				for (int j = 0; j < 4; j++) 
					character[charIndex + j] = image[imgIndex + j];				
			}
		}
		data[i] = character;
	}
}

