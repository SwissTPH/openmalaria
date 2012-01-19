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


#ifndef LINE_CHART_H
#define LINE_CHART_H

#include "gl_headers.h"
#include "math_headers.h"
#include "Display.h"
#include "TextureLoader.h"

class SkyBox;
class Scene;

class LineChart
{
	public:

		LineChart(DisplayMM* display, Color color);
		~LineChart();

		Color	diffuse;

		void	render();
		void	setData(float* newData, unsigned int count);


	private:

		inline	float2	getReflection(const float3 & vertex, float* mvm, 
													const float3 & lightX, const float3 & lightY, 
													const float3 & normal, const float3 & light,
													float facing)
		{
			float specularity;
			float2 out;
			float3 reflectedEye = vertex - 2.0f*facing*normal;
			reflectedEye /= reflectedEye.length();

			specularity = reflectedEye*light; 
			if (specularity < 0.0001f) 
			{
				specularity = 0.0001f;
				out.x = 0.5f + 0.7f*lightX*reflectedEye/specularity;
				out.y = 0.5f + 0.7f*lightY*reflectedEye/specularity;
			}
			else
			{
				out.x = 0.5f + 0.7f*lightX*reflectedEye;
				out.y = 0.5f + 0.7f*lightY*reflectedEye;
			}

			return out;
		}

		float					distance;
		DisplayMM*			display;
		SkyBox*				skyBox;
		Scene*				scene;
		TextureLoader textureLoader;
		unsigned int	sampleCount, specularTex;
		float*				data;
};

#endif

