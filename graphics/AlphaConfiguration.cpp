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


#include "AlphaConfiguration.h"
#include "macros.h"
#include "Scene.h"
#include "GraphicsBridge.h"

#include <cmath>  

//		c'tor
AlphaConfiguration :: AlphaConfiguration(DisplayMM* display)
:	DataConfiguration(display),
	updateTime(0), sampleSize(GraphicsBridge::sampleSize), chart(0)
{		
	chart = new FieldDisplay(display, 5, sampleSize, float3(5.0f, 5.0f, 5.0f));
	GraphicsBridge :: display = chart;	  	
}

//		d'tor
AlphaConfiguration :: ~AlphaConfiguration()
{ 
	if (chart) delete chart;
}

//		render
void AlphaConfiguration :: render()
{ 
	glPushMatrix();
	glColor4f(1,0,0,1);

	chart->render();		
  	
	glLoadIdentity();
	glTranslatef(0,1.6f,-1.6f);	
	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);     

	glPopMatrix();
}

void AlphaConfiguration :: update(float deltaT)
{
	chart->update(deltaT);	
}

