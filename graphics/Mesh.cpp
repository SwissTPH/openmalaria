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


#include "Mesh.h"
#include "macros.h"
#include "boinc_api.h"

Mesh :: Mesh(String textureDirectory)
:	textureDirectory(textureDirectory)
{}

void 
Mesh :: addSegment(Segment* s)
{   const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	
	segments.push_back(s);
	if (s->texture[0] == '*') s->glTexID = 0;
	else
	{
		retval = boinc_resolve_filename((s->texture).c_str(),imagefile,MAX_LENGTH);
		s->glTexID = texLoader.loadTexture2D(imagefile, RGBA_TEXTURE, GL_CLAMP);
		s->buildVertexArray();
	}
}

void 
Mesh :: render()
{
	int tex;
	ITERATE_LIST(Segment*, segments, s, ts,
		
		tex = (*s)->glTexID;
		if (tex)
		{
			glBindTexture(GL_TEXTURE_2D, (*s)->glTexID);
			glEnable(GL_TEXTURE_2D);
		}
		else glDisable(GL_TEXTURE_2D);

		ITERATE_LIST(Triangle*, (*s)->triangles, t, tt,
			(*t)->render();
	)	)
}
