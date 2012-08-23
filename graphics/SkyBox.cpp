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


#include "SkyBox.h"
#include "PreRenderedBox.h"
#include "GraphicsBridge.h"
#include "GL/glu.h"
#include "IL/il.h"
#include "IL/ilu.h"
#include "macros.h"
#include <iostream>
#include <cmath>
#include "boinc_api.h" 

SkyBox :: SkyBox(std::string directory)
:	diffuse(0), normal(0),
	sunrise(0), sunset(0),
	sun(0), glowOpacity(1.0f),	
	sunsetTime(0.535f), sunriseTime(0.98f),
	softShadows(true),
	preRendering(false),
	activeMode(true),
	cubeRenderer(this),
	controller(this),
	sunAngle(0.3f),
	moonPhase(1.2f),
	inclination((float)(PI*10.0f/180.0f)),
	timeDot(0.001f),
	time(5.4f),
	humidityIst(0.1f),
	humiditySoll(0.1f),
	skyColor(0.2f,0.4f,0.8f,1.0f),
	sunlightColor(0.93f,0.9f,0.8f,1.0f),
	ambientColor(0.33f,0.33f,0.33f,1.0f),
	shadowColor(0.03f,0.03f,0.03f,1.0f),
	hazeColor(0.2f,0.4f,1,0.3f),
	sunColor(0.93f,0.9f,0.8f,1.0f),
	currentColor(0.0f,0.0f,0.0f,0.0f),
	timeFract(0.0f),
	updateTotal(16), updateCurrent(0)
{
	loadTextures(directory);

	for (int i = 0; i < 8; i++)
		pipelineSwitches[i] = true;

	sunset0 = 2.0f*PI*sunsetTime - 0.08f;
	sunset1 = 2.0f*PI*sunsetTime;
	sunrise0 = 2.0f*PI*sunriseTime;
	sunrise1 = 2.0f*PI;

	if (GraphicsBridge :: preRenderedBoxResolution > 4)
		preRenderedBox = new PreRenderedBox(this, GraphicsBridge :: preRenderedBoxResolution);
	else
		preRenderedBox = new PreRenderedBox(this, 512);
}

void SkyBox :: loadTextures(std::string directory)
{	
	const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;
	glPixelTransferf(GL_RED_SCALE, 1.5f);
	glPixelTransferf(GL_GREEN_SCALE, 1.4f);
	glPixelTransferf(GL_BLUE_SCALE, 1.6f);
	diffuse = textureLoader.loadCubeMap("texture", RGBA_TEXTURE);
	glPixelTransferf(GL_RED_SCALE, 1.0f);
	glPixelTransferf(GL_GREEN_SCALE, 1.0f);
	glPixelTransferf(GL_BLUE_SCALE, 1.0f);

	gray = textureLoader.loadCubeMap("texture", DESATURATED_TEXTURE);
	normal = textureLoader.loadCubeMap("normal", RGBA_TEXTURE);
	sunrise = textureLoader.loadCubeMap("sunrise", GRAYSCALE_TEXTURE);
	sunset = textureLoader.loadCubeMap("sunset", GRAYSCALE_TEXTURE);
	haze = textureLoader.loadCubeMap("haze", GRAYSCALE_TEXTURE, DOME);
	sky = textureLoader.loadCubeMap("sky", GRAYSCALE_TEXTURE, FULL);

	retval = boinc_resolve_filename("sun.png",imagefile,MAX_LENGTH);
	sun = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE);

	retval = boinc_resolve_filename("moon.png",imagefile,MAX_LENGTH);
	moon = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE);
	
	retval = boinc_resolve_filename("moon_mask.png",imagefile,MAX_LENGTH);
	moonMask = textureLoader.loadTexture2D(imagefile, GRAYSCALE_TEXTURE);
	
	retval = boinc_resolve_filename("moon_shadow.png",imagefile,MAX_LENGTH);
	moonShadow = textureLoader.loadTexture2D(imagefile, GRAYSCALE_TEXTURE);
	
	retval = boinc_resolve_filename("glow.png",imagefile,MAX_LENGTH);
	glow = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE);
	
	retval = boinc_resolve_filename("starfield.png",imagefile,MAX_LENGTH);
	stars = textureLoader.loadCubeMap(imagefile, GRAYSCALE_TEXTURE, ALL_EQUAL);
	
	retval = boinc_resolve_filename("blank.png",imagefile,MAX_LENGTH);
	blank = textureLoader.loadTexture2D(imagefile, GRAYSCALE_TEXTURE);
	
	afterglowMask = textureLoader.generateCubeMap(Z_VALUE_MAP, 512);
}

void SkyBox :: render(Side side)
{
	int vPort[4];
	glGetIntegerv(GL_VIEWPORT, vPort);
	glViewport(0,0,preRenderedBox->size,preRenderedBox->size);

	glMatrixMode (GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluPerspective( 90, 1, 0.01f, 30.0f );

	glMatrixMode (GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	switch(side)
	{
		case NORTH: glRotatef(0,0,1,0); break;
		case WEST: glRotatef(270,0,1,0); break;
		case SOUTH: glRotatef(180,0,1,0); break;
		case EAST: glRotatef(90,0,1,0); break;
		case TOP: glRotatef(-90,1,0,0); break;
		case BOTTOM: glRotatef(90,1,0,0); break;
	}

	render();

	glPopMatrix();
	glMatrixMode (GL_PROJECTION);
	glPopMatrix();
	glMatrixMode (GL_MODELVIEW);

	preRenderedBox->readPixels(side);
	glViewport(vPort[0],vPort[1],vPort[2],vPort[3]);
}

void SkyBox :: update(float dt)
{ 
	if (!activeMode)
	{
		preRenderedBox->deltaT += dt;
		
		return;
	}

	humidityIst = humiditySoll;

	time += 50.0f*timeDot*dt;
	timeFract = time - (float)(int)(time);
	
	sunAngle = (float)(2.0f*timeFract*PI);
	moonPhase = (float)(time/29.5);
	moonPhase = moonPhase - (float)(int)(moonPhase);
	moonPhase = 2.0f*PI_F*moonPhase;

	moonFulness = (float)(1.0f - cos(moonPhase))/1.3f;
	moonlightOffset = (float)sin((PI + moonPhase)/2.0f);

	moonAngle = sunAngle - moonPhase;
	if (moonAngle < 0.0f) moonAngle += 2.0f*PI_F;
}

void SkyBox :: render()
{	 	
	if (!activeMode)
	{
		activeMode = true;
		preRendering = true;
		updateCurrent++;
		if (!preRenderedBox->initialized)
		{
			for (int i = 0; i < 6; i++) preRenderedBox->update();
			preRenderedBox->initialized = true;
		}
		if (updateCurrent == updateTotal)
		{
			preRenderedBox->update();
			updateCurrent = 0;
		}
		glBindTexture(GL_TEXTURE_CUBE_MAP, preRenderedBox->texCubeFront);
		cubeRenderer.render(PRE_RENDERED, 1.0f);
		box = preRenderedBox->texCubeFront;
		activeMode = false;
		preRendering = false;

		return;
	}
  
	box = 0;
	float epsilon = 0.05f;
	night = timeFract > sunsetTime && timeFract < sunriseTime;
	if (night)
	{ 		
		if (timeFract < sunsetTime + epsilon) nightiness = 1.0f - (sunsetTime + epsilon - timeFract)/epsilon;
		else if (timeFract > sunriseTime - epsilon) nightiness = (sunriseTime - timeFract)/epsilon;
		else nightiness = 1.0f;
	}
	else nightiness = 0.0f;

	float  rise = (sunrise0 + sunrise1)/2,
					set = (sunset0 + sunset1)/2;

	illuminationAngle = night ? moonAngle : sunAngle;
	illuminationAngle = PI*(illuminationAngle + 2*PI - rise)/(set + 2*PI_F - rise);

	saveSunLocation();
	
	controller.setHumidity(humidityIst);

	if (timeFract > (sunsetTime + sunriseTime) / 2.0f)
		controller.setEnvironment((timeFract - sunriseTime)/(sunsetTime + 1.0f - sunriseTime));	
	else controller.setEnvironment((timeFract + 1.0f - sunriseTime)/(sunsetTime + 1.0f - sunriseTime));
  
	glDepthRange(1.0f, 1.0f);	
	glDepthFunc(GL_ALWAYS);	
	glEnable(GL_DEPTH_TEST);

	cubeRenderer.render(SKY, 1.0f);

	glDepthRange(0.95f, 0.95f);
	glDepthFunc(GL_LEQUAL);

	cubeRenderer.render(HAZE_ON_SKY, 1.0f);

	glPushMatrix();

		setModelViewMatrix(sunAngle);
		glDepthRange(0.9f, 0.9f);

		cubeRenderer.render(AFTERGLOW, 1.0f);

	glPopMatrix();
    
	glDepthFunc(GL_GEQUAL);
	glDepthRange(0.1f, 0.1f);

	if (softShadows)
	{
		float	a0 = 1.0f/3.0f, a, 
					b0 = 2.0f/3.0f, b,
					c0 = 1.0f, c;
		if (night)
		{
			a = a0 * nightiness*moonFulness;
			b = b0 * nightiness*moonFulness;
			c = c0 * nightiness*moonFulness;
		}
		else 
		{
			a = a0; b = b0; c = c0;
		}

		angleBias = 1.0f/256.0f;
		currentColor.set((1.0f - a)*shadowColor + a*sunlightColor);
		cubeRenderer.render(SUNRISE, 1.0f);
		
		angleBias = 0.0f;
		currentColor.set((1.0f - b)*shadowColor + b*sunlightColor);
		cubeRenderer.render(SUNRISE, 1.0f);
		
		angleBias = -1.0f/256.0f;
		currentColor.set((1.0f - c)*shadowColor + c*sunlightColor);
		cubeRenderer.render(SUNRISE, 1.0f);
	
		angleBias = 1.0f/256.0f;
		currentColor.set((1.0f - b)*shadowColor + b*sunlightColor);
		cubeRenderer.render(SUNSET, 1.0f);	
		
		angleBias = 0.0f;
		currentColor.set((1.0f - a)*shadowColor+a*sunlightColor);
		cubeRenderer.render(SUNSET, 1.0f);	

		angleBias = -1.0f/256.0f;
		currentColor.set(shadowColor);
		cubeRenderer.render(SUNSET, 1.0f);	
	}
	else
	{
		angleBias = 0.0f;
		currentColor.set(night?
			(1.0f - nightiness)*shadowColor + nightiness*moonFulness*sunlightColor
			:	sunlightColor);
		cubeRenderer.render(SUNRISE, 1.0f);
		currentColor.set(shadowColor);
		cubeRenderer.render(SUNSET, 1.0f);
	}	

	if (pipelineSwitches[0]) 
	{
		if (night)
		glColor3f(	0.5f + nightiness*0.5f*cos(illuminationAngle),
								0.5f + nightiness*0.5f*sin(illuminationAngle),
								0.5f);
		else
		glColor3f(	0.5f + 0.5f*cos(illuminationAngle),
								0.5f + 0.5f*sin(illuminationAngle),
								0.5f);
		
		cubeRenderer.render(NORMALS, 1.0f);	
	}
	if (pipelineSwitches[1]) cubeRenderer.render(DIFFUSE, 1.0f);
	if (pipelineSwitches[2]) cubeRenderer.render(AMBIENT, 1.0f);	

	glDisable(GL_DEPTH_TEST);
	if (pipelineSwitches[3]) cubeRenderer.render(HAZE, 1.0f);	

	float newMoonFade = 0.45f, newMoonOpacity;
	if (moonPhase < newMoonFade) newMoonOpacity = moonPhase/newMoonFade;
	else if (moonPhase < 2.0f*PI - newMoonFade) newMoonOpacity = 1.0f;
	else newMoonOpacity = (2.0f*PI - moonPhase)/newMoonFade;

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	if (night) 
	{ 	 
			glPushMatrix();
			setModelViewMatrix(moonAngle);

			glDepthRange(0.8f, 0.8f);
			cubeRenderer.render(MOON_SHADOW, 1.0f);				
			glDepthRange(0.7f, 0.7f);
			currentColor.set(nightiness*sunColor + 0.3f*sunColor*(1.0f - nightiness));
			currentColor.a *= newMoonOpacity;
			cubeRenderer.render(MOON, 1.0f);

			glDisable(GL_DEPTH_TEST);

			float projectionMatrix[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, projectionMatrix);
			glowOpacity = preRendering ? 0.6f : -projectionMatrix[10];
			CLAMP(glowOpacity,0.0f,1.0f)
			glowOcclusion = getOcclusion(moonAngle);
			currentColor.set(0.7f*(sunColor + hazeColor));				
			currentColor *= glowOpacity*nightiness*glowOcclusion*moonFulness;
			glTranslatef(0,moonlightOffset/10,0);
			cubeRenderer.render(GLOW, 1.0f);

			glPopMatrix();
			setModelViewMatrix(sunAngle);	
			glEnable(GL_DEPTH_TEST);
			glDepthRange(0.5f, 0.5f);
			currentColor.set(sunColor);
			currentColor *= nightiness;
			cubeRenderer.render(STARS, 1.0f);	
	}
	else 
	{	
		glPushMatrix();
		setModelViewMatrix(sunAngle);
		glDepthRange(0.85f, 0.85f);
		cubeRenderer.render(SUN, 1.0f);

		glDisable(GL_DEPTH_TEST);   
		float projectionMatrix[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, projectionMatrix);
		glowOpacity = preRendering ? 0.6f : -projectionMatrix[10];
		glowOcclusion = getOcclusion(sunAngle);
		CLAMP(glowOpacity,0.0f,1.0f)
		currentColor.set(glowOcclusion*(1.2f*sunColor)*glowOpacity);
		cubeRenderer.render(GLOW, 1.0f); 

		glPopMatrix();
		
		glEnable(GL_DEPTH_TEST);

		glPushMatrix();	
		{
			setModelViewMatrix(moonAngle);	

			glDepthRange(0.75f, 0.75f);
			cubeRenderer.render(MOON_SHADOW, 1.0f);				
			glDepthRange(0.7f, 0.7f);
			currentColor.set(sunColor*0.3f);
			currentColor.a *= newMoonOpacity;
			cubeRenderer.render(MOON, 1.0f);

			glDisable(GL_DEPTH_TEST);
		}
		glPopMatrix();
	}
	
	glDepthFunc(GL_LEQUAL);
	cubeRenderer.reset();
}

void SkyBox :: assumeEnvironment(Environment t)
{ 	
	skyColor.set(t.sky);
	sunlightColor.set(t.sunlight);	
	ambientColor.set(t.ambient);
	shadowColor.set(t.shadow);
	hazeColor.set(t.haze);
	sunColor.set(t.sun);
}

void SkyBox :: setModelViewMatrix(float angle)
{
	glRotatef(90, 0,1.0f,0);	
	glRotatef(-180*angle/PI,(float)cos(inclination),sin(inclination),0);
}

void SkyBox :: saveSunLocation()
{
	glPushMatrix();
	glLoadIdentity();

	glRotatef(90, 0,1.0f,0);	
	glRotatef(-180*(night ? moonAngle : sunAngle)/PI,(float)cos(inclination),sin(inclination),0);

	float mvm[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, mvm);
  
	sunX.x = mvm[0]; sunX.y = mvm[1]; sunX.z = mvm[2];
	sunY.x = mvm[4]; sunY.y = mvm[5]; sunY.z = mvm[6];
	sunPosition.x = mvm[8]; sunPosition.y = mvm[9]; sunPosition.z = mvm[10];

	glPopMatrix();
}

float SkyBox :: getOcclusion(float angle)
{
	float out;
	if (angle < sunset0) out = 1.0f;
	else if (angle < sunset1) out = (sunset1 - angle)/(sunset1 - sunset0);
	else if (angle < sunrise0) out = 0.0f;
	else if (angle < sunrise1) out = (angle - sunrise0)/(sunrise1 - sunrise0);
	else out = 1.0f;
	return out;
}
