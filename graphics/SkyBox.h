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


#ifndef SKYBOX_H
#define SKYBOX_H

#include "gl_headers.h"
#include "CubeRenderer.h"
#include "Color.h"
#include "TextureLoader.h"
#include "Environment.h"
#include "EnvironmentController.h"
#include "float3.h"
#include <string>

class PreRenderedBox;

class SkyBox
{ 

	public:

		typedef enum {NORTH, WEST, SOUTH, EAST, TOP, BOTTOM} Side;

		SkyBox(std::string directory);

		bool		softShadows, night, preRendering, activeMode;
		bool		pipelineSwitches[8];

		int			updateTotal, updateCurrent;

		GLuint	diffuse, gray, normal, sunrise, sunset, haze, 
						sky, sun, moon, moonMask, moonShadow, glow, stars,
						blank, afterglowMask, box;

		float		sunAngle, moonAngle, illuminationAngle, moonPhase, 
						moonFulness, moonlightOffset, angleBias, inclination, 
						humiditySoll, humidityIst, glowOcclusion, time, timeDot, glowOpacity, 
						nightiness, timeFract;

		float3	sunPosition, sunX, sunY;

		Color		skyColor, sunlightColor, ambientColor, 
						shadowColor, hazeColor, sunColor,
						currentColor, afterglowColor;

		void update(float dt);
		void render();
		void render(Side side);
		void loadTextures(std::string directory);
		void assumeEnvironment(Environment t);

	private:

		float sunset0, sunset1, sunrise0, sunrise1,
					sunsetTime, sunriseTime;

		CubeRenderer cubeRenderer;
		EnvironmentController controller;
		TextureLoader textureLoader;
		PreRenderedBox* preRenderedBox;

		void setModelViewMatrix(float angle);
		void saveSunLocation();
		float getOcclusion(float angle);
};

#endif
