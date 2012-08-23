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


#include "Display.h"
#include "SkyBox.h"
#include "gl_headers.h"
#include "macros.h"
#include "AlphaConfiguration.h"

DisplayMM :: DisplayMM(SkyBox* skyBox, Scene* scene)
	:	skyBox(skyBox),
		scene(scene)
{
	myLines = new Line*[20];
	data = new AlphaConfiguration(this);
}

void DisplayMM :: render()
{ 
	float fogColor[4];
	skyBox->hazeColor.writeTo(fogColor);
	glFogi (GL_FOG_MODE, GL_EXP2);
  glFogfv (GL_FOG_COLOR, fogColor);
  glFogf (GL_FOG_DENSITY, 0.07f + 0.05f*fogColor[3]);
  glHint (GL_FOG_HINT, GL_NICEST);
	
	data->render();
	
}

void DisplayMM :: update(float deltaT)
{
	data->update(deltaT);
}

