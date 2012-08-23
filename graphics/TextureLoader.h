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


#ifndef TEXTURE_LOADER_H
#define TEXTURE_LOADER_H

#include "IL/il.h"
#include "IL/ilu.h"
#include "macros.h"
#include "gl_headers.h"
#include "Color.h"
#include <iostream>
#include <cmath>
#include <string>

typedef enum 
{ RGBA_TEXTURE, GRAYSCALE_TEXTURE, DESATURATED_TEXTURE, RGB_TEXTURE }
TextureType;

typedef enum
{ EMPTY, FULL, DOME, ALL_EQUAL }
Emptiness;

typedef enum
{	Z_VALUE_MAP, EMPTY_RGB_MAP }
ProceduralMapType;

class TextureLoader
{
	public:

		TextureLoader(){}

		GLuint	loadTexture2D(std::string filename, TextureType type, 
															GLint edge = GL_CLAMP_TO_EDGE); 

		GLuint	loadCubeMap(std::string prefix, TextureType type, 
															Emptiness emptiness = EMPTY);

		GLuint	generateCubeMap(ProceduralMapType type, unsigned int size);


	private:

		GLuint	generateAndBind(GLenum target);
		void		loadImage(std::string filename, GLenum target, TextureType type);
		void		generateImage(ProceduralMapType type, GLenum target, unsigned int size);
};

#endif
