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


#include "ProgressBar.h"
#include "macros.h"
#include "GraphicsBridge.h"
#include "gl_headers.h"
#include "config.h"

#define W_HALF 3.3f
#define H_HALF 0.2f

#if (defined(_GRAPHICS_6))
#include "ShmStruct.h"
extern UC_SHMEM* shmem;
#else
extern float fdone;
#endif

ProgressBar :: ProgressBar(unsigned int insideTexture, unsigned int outsideTexture)
:	inside(insideTexture), outside(outsideTexture),
	a(-W_HALF,  H_HALF, 0.0f),
	b( W_HALF,  H_HALF, 0.0f),
	c( W_HALF, -H_HALF, 0.0f),
	d(-W_HALF, -H_HALF, 0.0f),
	value(0.3f)
{
	GraphicsBridge::progressBar = this;
}

void ProgressBar :: render(Color in, Color out)
{
	value += 0.0002f;
	if (value > 1.0f) value = 0.0f;
#if (defined(_GRAPHICS_6))
		value = shmem->fraction_done;
#else
    value = fdone;
#endif
	float clamp = 0.15f;
	float correctedValue = (value + clamp)/(1.0f + 2*clamp);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	glBindTexture(GL_TEXTURE_2D, outside);	
	out.set();
	RENDER_IMAGE_QUAD(a,b,c,d);
	
	glBindTexture(GL_TEXTURE_2D, inside);
	in.set();
	glBegin(GL_QUADS);
		glTexCoord2f(0,0);
		glVertex3f(a.x, a.y, a.z);
		glTexCoord2f(correctedValue, 0);
		glVertex3f(a.x + (b.x - a.x)*correctedValue, b.y , b.z);
		glTexCoord2f(correctedValue, 1);
		glVertex3f(d.x + (c.x - d.x)*correctedValue, c.y, c.z);
		glTexCoord2f(0,1);
		glVertex3f(d.x, d.y, d.z);
	glEnd();
}
