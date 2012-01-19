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


#include "Segment.h"
#include "gl_headers.h"
#include "macros.h"

void Segment :: render()
{
	glBindTexture(GL_TEXTURE_2D, glTexID);
	
	/*glBegin(GL_TRIANGLES);
	for (unsigned int i = 0; i < vertexCount; i++)
	{
		glNormal3f(normalBuffer[3*i], normalBuffer[3*i + 1], normalBuffer[3*i + 2]);
		glTexCoord2f(texCoordBuffer[2*i], texCoordBuffer[2*i + 1]);
		glVertex3f(vertexBuffer[3*i], vertexBuffer[3*i + 1], vertexBuffer[3*i + 2]);
	}
	glEnd();*/
	
	glVertexPointer( 3, GL_FLOAT, 0, vertexBuffer );
	glNormalPointer( GL_FLOAT, 0, normalBuffer );
	glTexCoordPointer( 2, GL_FLOAT, 0, texCoordBuffer );

	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_TEXTURE_COORD_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );	

	glDrawArrays( GL_TRIANGLES, 0, 1 );

	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_TEXTURE_COORD_ARRAY );
	glDisableClientState( GL_NORMAL_ARRAY );
}

void Segment :: buildVertexArray()
{	
	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_TEXTURE_COORD_ARRAY );
	glEnableClientState( GL_NORMAL_ARRAY );	
		//Todo: resolve this (should also work on Unix)
		#ifdef _WIN32
		glGenBuffers(1, &vertexBufferID);
		glGenBuffers(1, &normalBufferID);
		glGenBuffers(1, &texCoordBufferID);
		#endif

		float* vertices = new float[triangles.size()*3*3];
		float* normals = new float[triangles.size()*3*3];
		float* texCoords = new float[triangles.size()*3*2];

		int i = 0;
		ITERATE_LIST(Triangle*, triangles, t, end,

			vertices[9*i + 0] = (*t)->vertex0.x;
			vertices[9*i + 1] = (*t)->vertex0.y;
			vertices[9*i + 2] = (*t)->vertex0.z;
			vertices[9*i + 3] = (*t)->vertex1.x;
			vertices[9*i + 4] = (*t)->vertex1.y;
			vertices[9*i + 5] = (*t)->vertex1.z;
			vertices[9*i + 6] = (*t)->vertex2.x;
			vertices[9*i + 7] = (*t)->vertex2.y;
			vertices[9*i + 8] = (*t)->vertex2.z;

			normals[9*i + 0] = (*t)->normal0.x;
			normals[9*i + 1] = (*t)->normal0.y;
			normals[9*i + 2] = (*t)->normal0.z;
			normals[9*i + 3] = (*t)->normal1.x;
			normals[9*i + 4] = (*t)->normal1.y;
			normals[9*i + 5] = (*t)->normal1.z;
			normals[9*i + 6] = (*t)->normal2.x;
			normals[9*i + 7] = (*t)->normal2.y;
			normals[9*i + 8] = (*t)->normal2.z;

			texCoords[6*i + 0] = (*t)->texture0.x;
			texCoords[6*i + 1] = (*t)->texture0.y;
			texCoords[6*i + 2] = (*t)->texture1.x;
			texCoords[6*i + 3] = (*t)->texture1.y;
			texCoords[6*i + 4] = (*t)->texture2.x;
			texCoords[6*i + 5] = (*t)->texture2.y;

			i++;
		)

		vertexCount = 3*i;

		/*glBindBuffer(GL_ARRAY_BUFFER, vertexBufferID);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, normalBufferID);
		glBufferData(GL_ARRAY_BUFFER, sizeof(normals), normals, GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, texCoordBufferID);
		glBufferData(GL_ARRAY_BUFFER, sizeof(texCoords), texCoords, GL_STATIC_DRAW);	 */

		vertexBuffer = vertices;
		normalBuffer = normals;
		texCoordBuffer = texCoords;

	glDisableClientState( GL_VERTEX_ARRAY );
	glDisableClientState( GL_TEXTURE_COORD_ARRAY );
	glDisableClientState( GL_NORMAL_ARRAY );

	/*	delete[] vertices; 
		delete[] normals;
		delete[] texCoords;	 */
}
