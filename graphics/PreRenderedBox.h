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


#ifndef PRE_RENDERED_BOX_H
#define PRE_RENDERED_BOX_H

#include "SkyBox.h"
#include "macros.h"
#include "gl_headers.h"
#include "TextureLoader.h"
#include <map>

class PreRenderedBox
{
	public:

		typedef std::map< SkyBox::Side, GLenum > Sidemap;
		typedef std::map< SkyBox::Side, SkyBox::Side > SideToSideMap;

		PreRenderedBox(SkyBox* skyBox, int size)
			:
			skyBox(skyBox), size(size),
			currentSide(SkyBox::EAST),
			deltaT(0.0f),
			initialized(false)
		{
			sidemap[SkyBox::EAST] = GL_TEXTURE_CUBE_MAP_POSITIVE_X;
			sidemap[SkyBox::WEST] = GL_TEXTURE_CUBE_MAP_NEGATIVE_X;
			sidemap[SkyBox::NORTH] = GL_TEXTURE_CUBE_MAP_POSITIVE_Z;
			sidemap[SkyBox::SOUTH] = GL_TEXTURE_CUBE_MAP_NEGATIVE_Z;
			sidemap[SkyBox::TOP] = GL_TEXTURE_CUBE_MAP_NEGATIVE_Y;
			sidemap[SkyBox::BOTTOM] = GL_TEXTURE_CUBE_MAP_POSITIVE_Y;

			
			sideOrder[SkyBox::EAST] = SkyBox::WEST;
			sideOrder[SkyBox::WEST] = SkyBox::NORTH;
			sideOrder[SkyBox::NORTH] = SkyBox::SOUTH;
			sideOrder[SkyBox::SOUTH] = SkyBox::TOP;
			sideOrder[SkyBox::TOP] = SkyBox::BOTTOM;
			sideOrder[SkyBox::BOTTOM] = SkyBox::EAST;

			sides = new unsigned short*[6];
			ITERATE(i, 6)
			{
				sides[i] = new unsigned short[size*size*3];
			}

			texCubeFront = textureLoader.generateCubeMap(EMPTY_RGB_MAP, size);
			texCubeBack = textureLoader.generateCubeMap(EMPTY_RGB_MAP, size);			
		}

		void readPixels(SkyBox::Side side);
		void update();

		SkyBox* skyBox;	 
		Sidemap sidemap; 
		SideToSideMap sideOrder;
		TextureLoader textureLoader;
		bool initialized;
		int size;
		float deltaT;
		SkyBox::Side currentSide;
		unsigned int texCubeFront, texCubeBack;
		unsigned short** sides;
};

#endif
