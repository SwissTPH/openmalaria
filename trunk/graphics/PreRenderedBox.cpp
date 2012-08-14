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


#include "PreRenderedBox.h"

void PreRenderedBox :: readPixels(SkyBox::Side side)
{
	glReadPixels( 0, 0, size, size, GL_RGB, GL_UNSIGNED_BYTE, sides[side]);
	glEnable(GL_TEXTURE_CUBE_MAP);
	glDisable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_CUBE_MAP, texCubeBack);
	glTexSubImage2D(sidemap[side], 0, 0, 0, size, size, GL_RGB, GL_UNSIGNED_BYTE, sides[currentSide]);
}

void PreRenderedBox :: update()
{ 		
	skyBox->render(currentSide);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	currentSide = sideOrder[currentSide]; 

	if (currentSide == SkyBox::EAST)
	{ 		
		unsigned int swap = texCubeFront;
		texCubeFront = texCubeBack;
		texCubeBack = swap;
		skyBox->activeMode = true;
		skyBox->update(deltaT/6);
		skyBox->activeMode = false;
		deltaT = 0.0f;
	}
}
