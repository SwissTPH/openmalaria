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


#include "GL_Window.h"
#include "Bridge.h"
#include "gl_headers.h"
#include "SystemTimer.h"
#include "macros.h"

GL_Window :: GL_Window(Bridge* bridge)
:	bridge(bridge), frameCount(0)
{
	overlay.scene = &scene;
	init();
}

void GL_Window :: init()
{
	glMatrixMode(GL_MODELVIEW);
	glShadeModel(GL_SMOOTH);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glMaterialf(GL_FRONT, GL_SHININESS, 128); 

	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glEnable(GL_NORMALIZE);

	glEnable(GL_TEXTURE_2D);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);	
	glCullFace(GL_BACK);

	glDisable(GL_LIGHT0);

	glInitNames();
	glPushName(-1);
	
	glEnable(GL_LINE_SMOOTH);
}

void GL_Window :: resize(int w, int h)
{
	glViewport(0,0,w,h);
}

void GL_Window :: render()
{
	scene.render();
	overlay.render();

	glFlush();

	frameCount++;
	if (frameCount == 20)
	{
		frameCount = 0;
		int msecs = SystemTimer::getMsecs();
		scene.fps = 20000.0f/(float)msecs;
	}
}
