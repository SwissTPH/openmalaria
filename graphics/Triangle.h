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


#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "float3.h"
#include "float2.h"

#include "gl_headers.h"

typedef struct
{
	float3 vertex, normal;
	float2 texture;
}  
Triplet;

#define RENDER_VERTEX(V, N, T) \
	glTexCoord2f(T.x, T.y); \
	glNormal3f(N.x, N.y, N.z); \
	glVertex3f(V.x, V.y, V.z); \


class Triangle
{
	public:

		Triangle(const Triplet & a, const Triplet & b, const Triplet & c)
		:	vertex0(a.vertex), normal0(a.normal), texture0(a.texture),
			vertex1(b.vertex), normal1(b.normal), texture1(b.texture),
			vertex2(c.vertex), normal2(c.normal), texture2(c.texture)
		{}

		float3	vertex0, normal0,
						vertex1, normal1,
						vertex2, normal2;

		float2	texture0, texture1, texture2;

		inline void render()
		{
			glBegin(GL_TRIANGLES);
			RENDER_VERTEX(vertex0, normal0, texture0)
			RENDER_VERTEX(vertex1, normal1, texture1)
			RENDER_VERTEX(vertex2, normal2, texture2)
			glEnd();
		}
};

#endif
