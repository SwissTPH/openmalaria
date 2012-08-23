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


#include "SceneController.h"
#include "Scene.h"
#include "Debug.h"
#include "macros.h"
#include "SkyBox.h"
#include "AlphaConfiguration.h"
#include "FieldDisplay.h"

SceneController :: SceneController(Scene* scene)
:	scene(scene)
{}

void SceneController :: keyPressed(Key k)
{
	float dc = 0.05f;
	char c = (char) k.character + 32;
	if (k.isCharacter)
	{
		switch (c)
		{ 
			case 'q':
				scene->skyBox->timeDot = -0.004f;
				break;

			case 'w':
				scene->skyBox->timeDot = 0.004f;
				break;

			case 'a':
				scene->skyBox->timeDot = 0.0f;
				break;

			case 't':
				scene->skyBox->softShadows ^= true;
				break;

			case 's':
				scene->saveScreenshot();
				break;

			case 'f':
				scene->dataDisplay->data->chart->soft ^= true;
				break;

			case 'm':
				scene->viewController->mosquitoCam ^= true;
				break;
		}
	}
	else
	{
		double speed = 0.05;
		switch (k.specialKey)
		{
			case F1:
				scene->skyBox->pipelineSwitches[0] ^= true;
				break;
			case F2:
				scene->skyBox->pipelineSwitches[1] ^= true;
				break;
			case F3:
				scene->skyBox->pipelineSwitches[2] ^= true;
				break;
			case F4:
				scene->skyBox->pipelineSwitches[3] ^= true;
				break;

			case F5:
				scene->switches[0] ^= true;
				break;
			case F6:
				scene->switches[1] ^= true;
				break;
			case F7:
				scene->skyBox->activeMode ^= true;
				break;

			case END:
				LOG("time: " << scene->skyBox->time << "\nsun: " << scene->skyBox->sunAngle)
				break;

			case SPACE:
				scene->overlayOn ^= true;
				break;

			case U_CURSOR:
				scene->deltaSDot.z = speed;
				break;
			case D_CURSOR:
				scene->deltaSDot.z = -speed;
				break;
			case L_CURSOR:
				scene->deltaSDot.x = speed;
				break;
			case R_CURSOR:
				scene->deltaSDot.x = -speed;
				break;

			case PG_UP:
				scene->fov -= 1.0;
				break;
			case PG_DN:
				scene->fov += 1.0;
				break;

			case INSERT_KEY:
				Debug :: doubles[0] += 10.0;
				break;

			case DELETE_KEY:
				Debug :: doubles[0] -= 10.0;
				break;
    }
	}
}

void SceneController :: keyReleased(Key k)
{
	char c = (char) k.character + 32;
	if (k.isCharacter)
	{
		switch (c)
		{			
			case 'w':
			case 'q':
				scene->skyBox->timeDot = 0.001f;
				break;
		}
	}
	else
	{
		switch (k.specialKey)
		{
			case U_CURSOR:
				scene->deltaSDot.z = 0.0;
				break;
			case D_CURSOR:
				scene->deltaSDot.z = 0.0;
				break;
			case L_CURSOR:
				scene->deltaSDot.x = 0.0;
				break;
			case R_CURSOR:
				scene->deltaSDot.x = 0.0;
				break;
		}
	}
}

void SceneController :: rotate(int2 relative)
{
	double phi = scene->phi, theta = scene->theta;
	phi += relative.x; theta += relative.y;
	CLAMP(theta, -90, 90)
	scene->phi = phi; scene->theta = theta;
}

void SceneController :: zoom(int2 relative)
{
	double r = scene->r;
	r -= 0.02f*relative.y;
	CLAMP(r, 0.3f, 8.0f)
	scene->r = r;
}
