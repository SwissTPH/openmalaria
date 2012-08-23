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


#include "LineChart.h"
#include "macros.h"
#include "math_headers.h"
#include "SkyBox.h"
#include "Scene.h"
#include "boinc_api.h"



//	BEGIN MACRO MESS:

#define RENDER_SPECULAR(X, Y, Z)	\
	vertex = float3(X, Y, Z); \
	eyeVertex.x = mvm[0]*vertex.x + mvm[4]*vertex.y + mvm[8]*vertex.z + mvm[12];\
	eyeVertex.y = mvm[1]*vertex.x + mvm[5]*vertex.y + mvm[9]*vertex.z + mvm[13];\
	eyeVertex.z = mvm[2]*vertex.x + mvm[6]*vertex.y + mvm[10]*vertex.z + mvm[14];\
  facing = normal*eyeVertex;\
	texCoord = getReflection(eyeVertex, mvm, lightY, eyeLightAxis, normal, eyeLight, facing); \
	glTexCoord2f(texCoord.x, texCoord.y); \
	glVertex3f(vertex.x, vertex.y, vertex.z);\

/*#define RENDER_SPECULAR_CAUTIOUSLY(X, Y, Z)	\
	vertex = float3(X, Y, Z); \
	eyeVertex.x = mvm[0]*vertex.x + mvm[4]*vertex.y + mvm[8]*vertex.z;\
	eyeVertex.y = mvm[1]*vertex.x + mvm[5]*vertex.y + mvm[9]*vertex.z;\
	eyeVertex.z = mvm[2]*vertex.x + mvm[6]*vertex.y + mvm[10]*vertex.z - distance;\
  facing = normal*eyeVertex;\
	if (facing <= 0) { \
		texCoord = getReflection(eyeVertex, mvm, lightY, eyeLightAxis, normal, eyeLight, facing); \
		glTexCoord2f(texCoord.x, texCoord.y); \
		glVertex3f(vertex.x, vertex.y, vertex.z);}\*/

#define RENDER_DIFFUSE(X, Y, Z)	\
	glColor3f(diffuseFinal.r, diffuseFinal.g, diffuseFinal.b);\
	glVertex3f(X, Y, Z);\


//	END MACRO MESS:

//		c'tor
LineChart :: LineChart(DisplayMM* display, Color color)
:	display(display),
	data(0),
	sampleCount(0),
	skyBox(display->skyBox),
	diffuse(color),
	scene(display->scene)
{
	const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	retval = boinc_resolve_filename("specular.png",imagefile,MAX_LENGTH);
	specularTex = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE, GL_CLAMP_TO_EDGE);
}

//		d'tor
LineChart :: ~LineChart()
{
	if (data) delete data;
}

//		render
void	LineChart :: render()
{
	distance = scene->r;
	float w = 3.0f/(float)(sampleCount - 1);
	float d = 0.1f, facing;
	
	float3	normal, vertex, eyeVertex, light, lightAxis, 
					eyeLightAxis, eyeLight, lightY, worldNormal;
	float2	texCoord;
	float*	mvm = new float[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, mvm);

	light = scene->light;
	
	lightAxis = float3(0, (float)sin(skyBox->inclination), -(float)cos(skyBox->inclination));
          	
	eyeLightAxis.x = mvm[0]*lightAxis.x + mvm[4]*lightAxis.y + mvm[8]*lightAxis.z;
	eyeLightAxis.y = mvm[1]*lightAxis.x + mvm[5]*lightAxis.y + mvm[9]*lightAxis.z;
	eyeLightAxis.z = mvm[2]*lightAxis.x + mvm[6]*lightAxis.y + mvm[10]*lightAxis.z;

	eyeLight.x = mvm[0]*light.x + mvm[4]*light.y + mvm[8]*light.z;
	eyeLight.y = mvm[1]*light.x + mvm[5]*light.y + mvm[9]*light.z;
	eyeLight.z = mvm[2]*light.x + mvm[6]*light.y + mvm[10]*light.z;

	lightY = eyeLight.cross(eyeLightAxis);
  
	float occlusion = scene->occlusion;
	float maxDeflection = -1000.1f;

	Color specularColor = skyBox->sunColor*occlusion,
				diffuseFinal, sun = skyBox->sunlightColor, 
				ambient = skyBox->ambientColor;

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, specularTex);
	glColor3f(COLOR_AS_ARGUMENT_NO_ALPHA(specularColor));

	glDisable(GL_BLEND);
	glDisable(GL_ALPHA_TEST);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_FOG);

	glBegin(GL_QUADS);
		normal = float3(-mvm[8], -mvm[9], -mvm[10]);
		if (normal.z >= maxDeflection)
		{
			RENDER_SPECULAR( -d, -1.0f, -1.5f)
			RENDER_SPECULAR( d, -1.0f, -1.5f)
			RENDER_SPECULAR( d, 1.333f*data[0] - 1.0f, -1.5f)
			RENDER_SPECULAR( -d, 1.333f*data[0] - 1.0f, -1.5f)
		}
		normal = float3(mvm[8], mvm[9], mvm[10]);
		if (normal.z >= maxDeflection)
		{
			RENDER_SPECULAR( -d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_SPECULAR( d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_SPECULAR( d, 1.333f*data[sampleCount-1] - 1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_SPECULAR( -d, 1.333f*data[sampleCount-1] - 1.0f, -1.5f + (sampleCount-1)*w)
		}
		normal = float3(-mvm[4], -mvm[5], -mvm[6]);
		if (normal.z >= maxDeflection)
		{
			RENDER_SPECULAR( -d, -1.0f, -1.5f)
			RENDER_SPECULAR( -d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_SPECULAR( d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_SPECULAR( d, -1.0f, -1.5f)
		}
		ITERATE(i, sampleCount - 1)
		{
			normal = float3(-mvm[0], -mvm[1], -mvm[2]);
			if (normal.z >= maxDeflection)
			{
				RENDER_SPECULAR( -d, -1.0f, -1.5f + i*w)
				RENDER_SPECULAR( -d, -1.0f, -1.5f + (i+1)*w)
				RENDER_SPECULAR( -d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)
				RENDER_SPECULAR( -d, 1.333f*data[i] - 1.0f, -1.5f + i*w)
			}

			normal = float3(mvm[0], mvm[1], mvm[2]);
			if (normal.z >= maxDeflection)
			{
				RENDER_SPECULAR( d, -1.0f, -1.5f + i*w)
				RENDER_SPECULAR( d, 1.333f*data[i] - 1.0f, -1.5f + i*w)
				RENDER_SPECULAR( d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)			
				RENDER_SPECULAR( d, -1.0f, -1.5f + (i+1)*w)
			}

			worldNormal = float3(0, w, data[i] - data[i+1]);
			worldNormal /= worldNormal.length();
			normal.x = mvm[4]*worldNormal.y + mvm[8]*worldNormal.z;
			normal.y = mvm[5]*worldNormal.y + mvm[9]*worldNormal.z;
			normal.z = mvm[6]*worldNormal.y + mvm[10]*worldNormal.z; 
			
			RENDER_SPECULAR( d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)
			RENDER_SPECULAR( -d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)			
			RENDER_SPECULAR( -d, 1.333f*data[i] - 1.0f, -1.5f + i*w)			
			RENDER_SPECULAR( d, 1.333f*data[i] - 1.0f, -1.5f + i*w)
			
		}
	glEnd();

	glDisable(GL_TEXTURE_2D);	
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
	float lambertCoefficient;	
	glDisable(GL_FOG);

	glBegin(GL_QUADS);
		normal = float3(-mvm[8], -mvm[9], -mvm[10]);
		if (normal.z >= maxDeflection)
		{
			lambertCoefficient = -light.x*occlusion;
			if (lambertCoefficient < 0.0f) lambertCoefficient = 0.0f;
			diffuseFinal = ambient*diffuse + lambertCoefficient*diffuse*sun;
			RENDER_DIFFUSE( -d, -1.0f, -1.5f)
			RENDER_DIFFUSE( d, -1.0f, -1.5f)
			RENDER_DIFFUSE( d, 1.333f*data[0] - 1.0f, -1.5f)
			RENDER_DIFFUSE( -d, 1.333f*data[0] - 1.0f, -1.5f)
		}
		normal = float3(mvm[8], mvm[9], mvm[10]);
		if (normal.z >= maxDeflection)
		{
			lambertCoefficient = -light.x*occlusion;
			if (lambertCoefficient < 0.0f) lambertCoefficient = 0.0f;
			diffuseFinal = ambient*diffuse + lambertCoefficient*diffuse*sun;
			RENDER_DIFFUSE( -d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_DIFFUSE( d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_DIFFUSE( d, 1.333f*data[sampleCount-1] - 1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_DIFFUSE( -d, 1.333f*data[sampleCount-1] - 1.0f, -1.5f + (sampleCount-1)*w)
		}
		normal = float3(-mvm[4], -mvm[5], -mvm[6]);
		if (normal.z >= maxDeflection)
		{
			lambertCoefficient = -light.x*occlusion;
			if (lambertCoefficient < 0.0f) lambertCoefficient = 0.0f;
			diffuseFinal = ambient*diffuse + lambertCoefficient*diffuse*sun;
			RENDER_DIFFUSE( -d, -1.0f, -1.5f)
			RENDER_DIFFUSE( -d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_DIFFUSE( d, -1.0f, -1.5f + (sampleCount-1)*w)
			RENDER_DIFFUSE( d, -1.0f, -1.5f)
		}
		ITERATE(i, sampleCount - 1)
		{
			normal = float3(-mvm[0], -mvm[1], -mvm[2]);
			if (normal.z >= maxDeflection)
			{
				lambertCoefficient = -light.x*occlusion;
				if (lambertCoefficient < 0.0f) lambertCoefficient = 0.0f;
				diffuseFinal = ambient*diffuse + lambertCoefficient*diffuse*sun;

				RENDER_DIFFUSE( -d, -1.0f, -1.5f + i*w)
				RENDER_DIFFUSE( -d, -1.0f, -1.5f + (i+1)*w)
				RENDER_DIFFUSE( -d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)
				RENDER_DIFFUSE( -d, 1.333f*data[i] - 1.0f, -1.5f + i*w)
			}

			normal = float3(mvm[0], mvm[1], mvm[2]);
			if (normal.z >= maxDeflection)
			{
				lambertCoefficient = light.x*occlusion;
				if (lambertCoefficient < 0.0f) lambertCoefficient = 0.0f;
				diffuseFinal = ambient*diffuse + lambertCoefficient*diffuse*sun;

				RENDER_DIFFUSE( d, -1.0f, -1.5f + i*w)
				RENDER_DIFFUSE( d, 1.333f*data[i] - 1.0f, -1.5f + i*w)
				RENDER_DIFFUSE( d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)			
				RENDER_DIFFUSE( d, -1.0f, -1.5f + (i+1)*w)
			}

			worldNormal = float3(0, w, data[i] - data[i+1]);
			worldNormal /= worldNormal.length();
			normal.x = mvm[4]*worldNormal.y + mvm[8]*worldNormal.z;
			normal.y = mvm[5]*worldNormal.y + mvm[9]*worldNormal.z;
			normal.z = mvm[6]*worldNormal.y + mvm[10]*worldNormal.z;
			lambertCoefficient = light*worldNormal*occlusion;
			if (lambertCoefficient < 0.0f) lambertCoefficient = 0.0f;
			diffuseFinal = ambient*diffuse + lambertCoefficient*diffuse*sun;

			RENDER_DIFFUSE( d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)
			RENDER_DIFFUSE( -d, 1.333f*data[i+1] - 1.0f, -1.5f + (i+1)*w)			
			RENDER_DIFFUSE( -d, 1.333f*data[i] - 1.0f, -1.5f + i*w)			
			RENDER_DIFFUSE( d, 1.333f*data[i] - 1.0f, -1.5f + i*w)
		}
	glEnd();	

	delete mvm;
}

//		setData
void	LineChart :: setData(float* newData, unsigned int count)
{
	if (data) delete data;
	data = new float[count];
	ITERATE(i, count) data[i] =	newData[i];
	sampleCount = count;
}


