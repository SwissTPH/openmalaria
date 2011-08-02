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


#include "ViewController.h"
#include "Anopheles.h"
#include "Scene.h"
#include "SkyBox.h"
#include "gl_headers.h"

ViewController :: ViewController(Anopheles* anopheles, Scene* scene)
:	anopheles(anopheles),
	scene(scene),
	roll(0), bearing(0),
	mosquitoCam(false),
	camera(anopheles->position),
	time(PI)
{
	bearingLine = SurfaceProvider :: getInstance() -> getLine();
	bearingTgtLine = SurfaceProvider :: getInstance() -> getLine();
	deltaBearingLine = SurfaceProvider :: getInstance() -> getLine();
}

void ViewController :: update(float deltaT)
{
	bearing = anopheles->bearing;
	roll = anopheles->bearingDot/2;
	
	time += deltaT/40.0f;
	if (time > 2.0f*PI) time -= 2.0f*PI;
	float h = cos(time);
	if (h < 0) h = 0;
	h *= h;
	scene->skyBox->humiditySoll = h;
}

void ViewController :: setNearView()
{
	if (mosquitoCam) 
	{
		glTranslatef(0, 0, -0.4f);
		glRotatef(-roll, 0, 0, 1 );		
		glRotatef(180-bearing, 0, 1, 0 );
		glTranslatef(-anopheles->position.x,-anopheles->position.y,-anopheles->position.z);
	}
	else
	{
		glTranslatef(0, 0, -scene->r - scene->overlayPresence*2.0f);
		glRotatef(scene->theta, 1, 0, 0 );
		glRotatef(scene->phi, 0, 1, 0 );
	}		
}

void ViewController :: unrotate(float additional)
{
	if (mosquitoCam) 
	{ 
		glRotatef(-additional, 0, 1, 0);
		glRotatef(bearing - 180, 0, 1, 0 );
		glRotatef(roll, 0, 0, 1 );			
	}
	else
	{ 		
		glRotatef(-additional, 0, 1, 0);
		glRotatef(-scene->phi, 0, 1, 0 );
		glRotatef(-scene->theta, 1, 0, 0 );		
	}		
}

void ViewController :: setFarView()
{
	if (mosquitoCam) glRotatef(-roll, 0, 0, 1 );	
	else glRotatef(scene->theta, 1, 0, 0 );
	glRotatef(mosquitoCam?180-bearing:scene->phi, 0, 1, 0 );
}

void ViewController :: renderDebugInfo()
{
	glDisable( GL_TEXTURE_2D );
	glBegin(GL_LINES);	
	float3 tgt(anopheles->target);
	float3 vel(anopheles->position + anopheles->velocity);
	float3 acc(anopheles->position + anopheles->acceleration);
	float3 pos(anopheles->position);
	glColor3f(1,0,0);
	glVertex3f(pos.x, pos.y, pos.z);
	glVertex3f(tgt.x, tgt.y, tgt.z);
	glColor3f(0,0.8f,0);
	glVertex3f(pos.x, pos.y, pos.z);
	glVertex3f(vel.x, vel.y, vel.z);
	glColor3f(1,0.9f,0);
	glVertex3f(pos.x, pos.y, pos.z);
	glVertex3f(acc.x, acc.y, acc.z);
	glEnd();

	glEnable(GL_TEXTURE_2D);
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPushMatrix();
	glLoadIdentity();

	bearingLine->print((int)anopheles->bearing);
	bearingTgtLine->print((int)anopheles->bearingTgt);
	deltaBearingLine->print((int)anopheles->bearingDot);

	glColor3f(1,0,0);
	glTranslatef(0,0.4f,-5);
	bearingLine->render(float2(0.17f, 0.23f), float2(0.5f,0.0f));

	glColor3f(0,0.8f,0);
	glTranslatef(0,0.3f,0);
	bearingTgtLine->render(float2(0.17f, 0.23f), float2(0.5f,0.0f));

	glColor3f(1,0.9f,0);
	glTranslatef(0,0.3f,0);
	deltaBearingLine->render(float2(0.17f, 0.23f), float2(0.5f,0.0f));

	glPopMatrix();
}
