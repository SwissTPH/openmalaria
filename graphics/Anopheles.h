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


#ifndef ANOPHELES_H
#define ANOPHELES_H

#include "Mesh.h"


	#define RENDER_SEGMENT( S ) \
	glBindTexture(GL_TEXTURE_2D, S->glTexID); \
	i = S->triangles.begin(); \
	t = S->triangles.end(); \
	while (i != t) \
	{ (*i)->render(); i++;} \


/*
#define RENDER_SEGMENT( S ) \
	glBindTexture(GL_TEXTURE_2D, S->glTexID); \
	glBindBuffer( GL_ARRAY_BUFFER, S->vertexBufferID ); \
	glVertexPointer( 3, GL_FLOAT, 0, NULL ); \
	glBindBuffer( GL_ARRAY_BUFFER, S->normalBufferID ); \
	glNormalPointer( GL_FLOAT, 0, NULL ); \
	glBindBuffer( GL_ARRAY_BUFFER, S->texCoordBufferID ); \
	glTexCoordPointer( 2, GL_FLOAT, 0, NULL ); \
	glDrawArrays( GL_TRIANGLES, 0, S->vertexCount );
*/
/*
#define RENDER_SEGMENT( S ) \
	glBindTexture(GL_TEXTURE_2D, S->glTexID); \
	glVertexPointer( 3, GL_FLOAT, 0, S->vertexBuffer ); \
	glNormalPointer( GL_FLOAT, 0, S->normalBuffer ); \
	glTexCoordPointer( 2, GL_FLOAT, 0, S->texCoordBuffer ); \
	glDrawArrays( GL_TRIANGLES, 0, S->vertexCount );
*/

//#define RENDER_SEGMENT( S ) S->render();

class Anopheles
{
	public:

		static Mesh* mesh;
		static void init();

		Anopheles(float3 position)
		:	position(position),
			target(position),
			velocity(0.0f, 0.0f, 0.0f),
			acceleration(0.0f, 0.0f, 0.0f),
			bearing(PLUS_RAND(360)), 			
			bearingDot(0.0f), 
			bearingDotDot(0.0f), 
			bearingTgt(0.0f),
			elevation(0.0f),
			elevationDot(0.0f),
			viscosity(0.04f)
		{}

		float3 position, velocity, acceleration, target;
		float bearing, elevation, bearingDot, bearingDotDot, bearingTgt, elevationDot, viscosity;

		inline void update(float deltaT)
		{			
			float tau = exp(-viscosity*deltaT);
			float3 newVelocity = (velocity + deltaT*acceleration)*tau;
			position += 0.5f * deltaT * (velocity + newVelocity);
			float d = position.length();
			const float innerBarrier = 3.5f;
			bool bounced = false;
			if (d < innerBarrier && d > 0.0001f)
			{
				position = innerBarrier*position/d;
				bounced = true;
			}
			velocity = newVelocity;

			float3 dir = velocity;
			float len = sqrt(dir.x*dir.x + dir.z*dir.z);
			if (len > 0) 
			{
				dir /= len;
				bearingTgt = asin(dir.x);
				if (dir.z < 0) bearingTgt = 180.0 - 180.0f*bearingTgt/PI;
				else bearingTgt = 180.0f*bearingTgt/PI;
				if (bearingTgt < 0.0f) bearingTgt += 360.0f;
			}
      			
			float scale = 0.6f;
			float maxRot = 120.0f;
			float slide = bearingTgt - bearing;
			if (slide > 180.0f) slide = slide - 360.0f;
			else if (slide < -180.0f) slide = 360.0f + slide;
			bearingDotDot = slide/scale;     
			bearingDot += bearingDotDot*deltaT;
			CLAMP(bearingDot, -maxRot, maxRot)

			bearing += deltaT * bearingDot;
			if (bearing > 360.0f) bearing -= 360.0f;
			else if (bearing < 0.0f) bearing += 360.0f;			

			elevation += deltaT * elevationDot;
			float closeness = (target - position).length();

			if (closeness < 0.8 || RANDOM() > 0.999f || bounced) 
			{
				float ang = PI*bearing/180.0f;
				target = float3(SYMM_RAND_3(8.0f));
				d = target.length();
				if (d < innerBarrier) 
				{
					if ( d > 0.0001f) target = innerBarrier*target/d;
					else target = float3(innerBarrier, 0.0f, 0.0f);
				}
				float lenT = target.length(), lenP = 4.2f;
				if (lenT > lenP) target = lenP*target/lenT;
			}

			acceleration = target - position;
			acceleration /= 4.0f;
			acceleration += float3(SYMM_RAND_3(0.8f));		
		}

		inline void render()
		{   						
			glPushMatrix();

			glEnable(GL_CULL_FACE);

				glTranslatef(position.x, position.y, position.z);

				glRotatef(elevation,1,0,0);
				glRotatef(bearing,0,1,0);
				glRotatef(-bearingDot/2,0,0,1);

				TriangleList::iterator i, t;
	      
				RENDER_SEGMENT(torso)
				RENDER_SEGMENT(abdomen)
				RENDER_SEGMENT(head)
				
				glDisable(GL_CULL_FACE);

				float wingAngle = RANDOM()*120.0f - 60.0f;
				glRotatef(wingAngle, 0, 0, 1);
				RENDER_SEGMENT(leftWing)
				glRotatef(-2.0f*wingAngle, 0, 0, 1);
				RENDER_SEGMENT(rightWing)
      
			glPopMatrix();
		}

	private:

		static Segment *abdomen, *head, *torso, *leftWing, *rightWing;
		
};

#endif
