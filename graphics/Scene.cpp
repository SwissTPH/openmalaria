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


#include "Scene.h"
#include "SkyBox.h"
#include "GraphicsBridge.h"
#include "gl_headers.h"
#include "GL/glu.h"
#include "StringWrapper.h"

//		c'tor
Scene :: Scene()
: phi(-30), theta(27), r(5.6f), rDot(0.0f),
	fov(98.0),
	frames(0), fps(0),
	deltaS(0,0,0),
	deltaSDot(0,0,0),
	controller(this),
	screenshotIndex(0),
	skyBox(new SkyBox(GraphicsBridge::imagePath + "savanna")),
	anophelesCount(40),
	overlayPresence(1.0f),
	overlayOn(true)
{
	dataDisplay = new DisplayMM(skyBox, this);	

	ITERATE(i, anophelesCount)
	{
		anopheles.push_back(
			new Anopheles(float3(SYMM_RAND(3.0f),SYMM_RAND(3.0f),SYMM_RAND(3.0f))));
	}

	ITERATE(i, 3) { switches[i] = true; }

	viewController = new ViewController(anopheles[0], this);
}

//		render
void Scene :: render()
{  	
	float fTime = 1.0f/(float)fps;
	if (fTime > 0.3f) fTime = 0.3f;
	else if (fTime < 0.002f) fTime = 0.002f;
	deltaS += deltaSDot;
	r += rDot;
	rDot *= 0.9f;
	viewController->update(fTime);

	float lambda = 1.0f - exp(-fTime*30.0f);
	if (overlayOn)
	{
		overlayPresence = 1.0f - (1.0f - overlayPresence)*lambda;
	}
	else
	{
		overlayPresence = (overlayPresence)*lambda;
	}
	
  int vPort[4];
	glGetIntegerv(GL_VIEWPORT, vPort);

	int y0 = (int)(overlayPresence*(vPort[3]/5));
	int y1 = vPort[3] - y0;

	glViewport(vPort[0], y0, vPort[2], y1 - y0);

	fov = 96 + 9*overlayPresence;
	setPerspectiveMatrix((float)(vPort[2] - vPort[0])/(float)(y1 - y0));
	glClearColor(0,0,0,1);
	glClearDepth(0.0f);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	
	glLoadIdentity();		
	glPushMatrix();  	
		viewController->setFarView();
		if (switches[0]) 
		{
			skyBox->update(fTime);
			skyBox->render();
		}
	glPopMatrix();

	glClearDepth(0.0f);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthRange(0.5f, 1.0f);

	float illumAngle = skyBox->night ? skyBox->moonAngle : skyBox->sunAngle;
	light = float3(	cos(illumAngle),
									cos(skyBox->inclination)*sin(illumAngle),
									sin(skyBox->inclination)*sin(illumAngle));
	occlusion = skyBox->night ? skyBox->nightiness*skyBox->glowOcclusion : skyBox->glowOcclusion;

	viewController->setNearView();
	glPushMatrix();
  
		dataDisplay->update(fTime);
		if (switches[1]) dataDisplay->render();

	glPopMatrix();

	float lite[4];
	
	Color specular = skyBox->sunColor*occlusion,
				diffuse = skyBox->sunlightColor*occlusion, 
				ambient = skyBox->ambientColor;
	specular.writeTo(lite); lite[3] = 1.0f;
	glLightfv(GL_LIGHT0, GL_SPECULAR, lite);
	diffuse.writeTo(lite); lite[3] = 1.0f;
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lite);
	ambient.writeTo(lite); lite[3] = 1.0f;
	glLightfv(GL_LIGHT0, GL_AMBIENT, lite);
  	
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LESS);

	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_FOG);  
	
	glColor4f(1,1,1,1);	
	lite[0] = light.z;
	lite[1] = light.x;
	lite[2] = -light.y;
	lite[3] = 0.0f;	
	glLightfv(GL_LIGHT0, GL_POSITION, lite);

	/*glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_TEXTURE_COORD_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );		*/
  
	glEnable(GL_TEXTURE_2D);

	if (switches[2])
	{
		ITERATE(i, anophelesCount)
		{
			anopheles[i]->update(fTime);
			anopheles[i]->render();
		}
	}

	glDisable(GL_LIGHTING);
	glDisable(GL_FOG);

	//viewController->renderDebugInfo();

	/*glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_TEXTURE_COORD_ARRAY );
	glDisableClientState( GL_NORMAL_ARRAY ); */

	glViewport(vPort[0], vPort[1], vPort[2], vPort[3]);
}

//		setPerspectiveMatrix
void Scene :: setPerspectiveMatrix(float aspect)
{ 	
	glMatrixMode (GL_PROJECTION);

	glLoadIdentity();
	gluPerspective( fov, aspect, 0.01f, 200.0f );

	glMatrixMode (GL_MODELVIEW);
}

//	 saveScreenshot
void Scene :: saveScreenshot()
{
	unsigned short* data = new unsigned short[1024*768*3];
	glReadPixels( 0, 0, 1024, 768, GL_RGB, GL_UNSIGNED_BYTE, data);

	ILuint image;
  ilGenImages(1, &image);
	ilBindImage(image);

	ilTexImage( 1024, 768, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, data);

	String filename = "screenshotxxxx.png";
	filename[13] = '0' + screenshotIndex%10;
	filename[12] = '0' + (screenshotIndex/10)%10;
	filename[11] = '0' + (screenshotIndex/100)%10;
	filename[10] = '0' + (screenshotIndex/1000)%10;
	ilSaveImage(ILstring(filename.c_str()));
	screenshotIndex++;

	ilDeleteImages(1, &image);
	delete data;
}

