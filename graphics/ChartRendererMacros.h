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



#define PUSH_NORMAL(X1,X0,Y1,Y0) \
	normals.push_back(-((vertices[X1] - vertices[X0]).cross(vertices[Y1] - vertices[Y0])).direction())

#define RENDER_DEPTH(V) \
	dist = (V - origin).length(); \
  glColor4f(0,0,0, 0);	\
	glVertex3f(V.x, V.y, V.z);

#define RENDER_FRONT(V) \
	dist = (V - origin).length(); \
	if (dist > 0.99f) dist = 0.99f; \
	glColor4f(0,0,0, (dist - dMin)/dSpan);	\
	glVertex3f(V.x, V.y, V.z);

#define RENDER_BACK(V) \
	dist = (V - origin).length(); \
	if (dist < 0.01f) dist = 0.01f; \
	glColor4f(0,0,0, (dist - dMin)/dSpan);	\
	glVertex3f(V.x, V.y, V.z);

#define RENDER_DEPTH_DEBUG(V) \
	dist = (V - origin).length(); \
  glColor4f(dist*calib,dist*calib,dist*calib, dist*calib);	\
	glVertex3f(V.x, V.y, V.z);

#define RENDER_SPECULAR(V, N) \
	in = V - origin; \
	out = in - 2.0f*(N*in)*N; \
	out /= out.length(); \
	facing = out*sun;  \
	fresnel = 1.0f + in*N/in.length();	\
	fresnel = 0.2f + 0.8f * fresnel;	\
	fresnel = fresnel*fresnel;	\
	if (N*sun < 0.0f) fresnel = 0.0f; \
  if (facing < 0.0001f) \
  {	\
		fresnel = 0.0f; \
    facing = 0.6f/0.0001f;	\
		texCoord.x = 0.5f - facing*lightX*out;	\
		texCoord.y = 0.5f - facing*lightY*out;	\
	}	\
	else \
	{	\
		texCoord.x = 0.5f + 0.6f*lightX*out;	\
		texCoord.y = 0.5f + 0.6f*lightY*out;	\
	}	\
	glColor4f(fresnel*specularColor.r, fresnel*specularColor.g, fresnel*specularColor.b, 1.0f);	\
	glTexCoord2f(texCoord.x, texCoord.y); \
	glVertex3f(V.x, V.y, V.z);

#define RENDER_ENVIRONMENT(V, N)	\
	in = V - origin; \
	out = in - 2.0f*(N*in)*N; \
	fresnel = 1.0f + in*N/in.length();	\
	/*fresnel = 0.2f + 0.8f * fresnel;*/	\
	fresnel = 0.3f + 0.7f * fresnel;	\
	glColor4f(ambientColor.r, ambientColor.g, ambientColor.b, fresnel);	\
	/*glColor4f(0.0f, 0.25f, 0.6f, fresnel);*/	\
	glTexCoord3f(out.x, z_scale*out.y, out.z); \
	glVertex3f(V.x, V.y, V.z); 

#define RENDER_DIFFUSE(V, N)	\
	/*glTexCoord2f(0.5f, 0.333f - V.y/3.0f);*/ \
	glVertex3f(V.x, V.y, V.z);
