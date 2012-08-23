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


#include "DemoConfiguration.h"
#include "macros.h"
#include "Scene.h"
#include <cmath>

//		c'tor
DemoConfiguration :: DemoConfiguration(DisplayMM* display)
:	DataConfiguration(display),
	title(0),
	chart1(display, Color(0.0f, 0.4f, 0.0f, 1)),
	chart2(display, Color(0.8f, 0.8f, 0.0f, 1)),
	chart3(display, Color(0.8f, 0.0f, 0.0f, 1))
{
	unsigned int samples = 25;
	float* data = new float[samples];
	ITERATE(i, samples) data[i] = RANDOM();
	chart1.setData(data, samples);
	delete data;

	samples = 13;
	data = new float[samples];
	ITERATE(i, samples) data[i] = RANDOM();
	chart2.setData(data, samples);
	delete data;

	samples = 33;
	data = new float[samples];
	ITERATE(i, samples) data[i] = RANDOM();
	chart3.setData(data, samples);
	delete data;

	title = SurfaceProvider :: getInstance() -> getLine();
	title->print("demo 04");
}

//		d'tor
DemoConfiguration :: ~DemoConfiguration()
{ 
	if (title) delete title;
}

//		render
void DemoConfiguration :: render()
{
	glPushMatrix();
	glTranslatef(-0.7f,-0.333f,0);
	chart1.render();
	glTranslatef(0.7f,0,0);
	chart2.render();
	glTranslatef(0.7f,0,0);
	chart3.render();
	glLoadIdentity();
	glTranslatef(0,1.3f,-1.6f);
	glColor4f(1,0.9f,0.2f,1); 
	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	int fps = (int)display->scene->fps;
	String s = "fps:     ";	 
	int i = 8;
	while (fps > 0){
		s[i] = '0' + fps%10;
		fps /= 10;
		i--;
	}
	title->print(s);
	title->render(float2(0.17f, 0.23f), float2(0.5f,0.0f));
	glPopMatrix();
}

