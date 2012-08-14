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


#include "CubeRenderer.h"
#include "SkyBox.h"
#include "gl_headers.h"
#include <cmath>

CubeRenderer :: CubeRenderer(SkyBox* skyBox)
:	skyBox(skyBox),
	texturizer(skyBox)
{
	vertices[0] = 1.0f;			vertices[1] = 1.0f;			vertices[2] = 1.0f;
	vertices[3] = 1.0f;			vertices[4] = -1.0f;		vertices[5] = 1.0f;	
	vertices[6] = -1.0f;		vertices[7] = -1.0f;		vertices[8] = 1.0f;
	vertices[9] = -1.0f;		vertices[10] = 1.0f;		vertices[11] = 1.0f;

	vertices[12] = 1.0f;		vertices[13] = 1.0f;		vertices[14] = 1.0f;
	vertices[15] = 1.0f;		vertices[16] = -1.0f;		vertices[17] = 1.0f;	
	vertices[18] = 1.0f;		vertices[19] = -1.0f;		vertices[20] = -1.0f;
	vertices[21] = 1.0f;		vertices[22] = 1.0f;		vertices[23] = -1.0f;

	vertices[24] = 1.0f;		vertices[25] = 1.0f;		vertices[26] = 1.0f;
	vertices[27] = 1.0f;		vertices[28] = 1.0f;		vertices[29] = -1.0f;	
	vertices[30] = -1.0f;		vertices[31] = 1.0f;		vertices[32] = -1.0f;
	vertices[33] = -1.0f;		vertices[34] = 1.0f;		vertices[35] = 1.0f;

	vertices[36] = 1.0f;		vertices[37] = 1.0f;		vertices[38] = -1.0f;
	vertices[39] = 1.0f;		vertices[40] = -1.0f;		vertices[41] = -1.0f;	
	vertices[42] = -1.0f;		vertices[43] = -1.0f;		vertices[44] = -1.0f;
	vertices[45] = -1.0f;		vertices[46] = 1.0f;		vertices[47] = -1.0f;

	vertices[48] = -1.0f;		vertices[49] = 1.0f;		vertices[50] = 1.0f;
	vertices[51] = -1.0f;		vertices[52] = -1.0f;		vertices[53] = 1.0f;	
	vertices[54] = -1.0f;		vertices[55] = -1.0f;		vertices[56] = -1.0f;
	vertices[57] = -1.0f;		vertices[58] = 1.0f;		vertices[59] = -1.0f;

	vertices[60] = 1.0f;		vertices[61] = -1.0f;		vertices[62] = 1.0f;
	vertices[63] = 1.0f;		vertices[64] = -1.0f;		vertices[65] = -1.0f;	
	vertices[66] = -1.0f;		vertices[67] = -1.0f;		vertices[68] = -1.0f;
	vertices[69] = -1.0f;		vertices[70] = -1.0f;		vertices[71] = 1.0f;
}

void CubeRenderer :: render(CubeMode mode, float w)
{
	texturizer.setUp(mode);

	switch (mode)
	{
		case NORMALS			:			
		case HAZE					:
		case HAZE_ON_SKY	:
		case SKY					:		
		case DIFFUSE			:
		case AMBIENT			:
		case STARS				:
		case SUNRISE			:
		case SUNSET				:		
		case AFTERGLOW		:
		{      			
			glBegin(GL_QUADS);
       
			for (int i = 0; i < 24; i++)
			{ 
				glTexCoord3f( vertices[3*i], vertices[3*i+1], vertices[3*i+2]);
				glVertex4f(	vertices[3*i], vertices[3*i+1], vertices[3*i+2], w);
			}

			glEnd();
		}
		break;

		case PRE_RENDERED :
		{      			
			glBegin(GL_QUADS);
       
			for (int i = 0; i < 24; i++)
			{ 
				glTexCoord3f( vertices[3*i], -vertices[3*i+1], vertices[3*i+2]);
				glVertex4f(	vertices[3*i], vertices[3*i+1], vertices[3*i+2], w);
			}

			glEnd();

			glDisable(GL_TEXTURE_CUBE_MAP);
			glEnable(GL_TEXTURE_2D);
		}
		break;

		/*case SUNRISE			:
		case SUNSET				:
		{      			
			glBegin(GL_QUADS);
       
			for (int i = 0; i < 24; i++)
			{ 
				glMultiTexCoord3f( GL_TEXTURE0, vertices[3*i], vertices[3*i+1], vertices[3*i+2]);
				glMultiTexCoord3f( GL_TEXTURE1, vertices[3*i], vertices[3*i+1], vertices[3*i+2]);
				glMultiTexCoord3f( GL_TEXTURE2, vertices[3*i], vertices[3*i+1], vertices[3*i+2]);
				glVertex4f(	vertices[3*i], vertices[3*i+1], vertices[3*i+2], w);
			}

			glEnd();
		}
		break;

		case SUNSET				:
			break;	*/

		case SUN:
		{
			float s = mode == MOON ? 0.10f : 0.23f;
			glBegin(GL_QUADS);
                                        
			glTexCoord2f(0,0);
			glVertex4f(-s,  s, 1.0f, w);
			glTexCoord2f(1,0);
			glVertex4f( s,  s, 1.0f, w);
			glTexCoord2f(1,1);
			glVertex4f( s, -s, 1.0f, w);
			glTexCoord2f(0,1);
			glVertex4f(-s, -s, 1.0f, w);

			glEnd();  
		}
		break;

		case MOON:
		{
			float s = mode == MOON ? 0.10f : 0.23f;
			glBegin(GL_QUADS);
                                        
			glTexCoord2f(0,0);
			glVertex4f(-s,  s, 1.0f, w);
			glTexCoord2f(1,0);
			glVertex4f( s,  s, 1.0f, w);
			glTexCoord2f(1,1);
			glVertex4f( s, -s, 1.0f, w);
			glTexCoord2f(0,1);
			glVertex4f(-s, -s, 1.0f, w);

			glEnd();

			glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_TRUE);

			glDepthRange(0.3, 0.3);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.02f);
			glEnable(GL_TEXTURE_2D);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glColor4f(0.0f,0.0f,0.0f,1.0f);

			glBegin(GL_QUADS);
      
				glTexCoord2f(0,0); glVertex4f(-s,  s, 1.0f, w);
				glTexCoord2f(1,0); glVertex4f( s,  s, 1.0f, w);
				glTexCoord2f(1,1); glVertex4f( s, -s, 1.0f, w);
				glTexCoord2f(0,1); glVertex4f(-s, -s, 1.0f, w);

			glEnd();

			glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
		}
		break;

		case MOON_SHADOW:
		{
			float phase = skyBox->moonPhase;
			float s = 0.10f, x0;			
      
			glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_TRUE);

			glDisable(GL_TEXTURE_2D);
			glDisable(GL_BLEND);			
			glColor4f(0.0f,0.0f,0.0f,0.0f);

			glBegin(GL_QUADS);
      
				glTexCoord2f(0,0); glVertex4f(-s,  s, 1.0f, w);
				glTexCoord2f(1,0); glVertex4f( s,  s, 1.0f, w);
				glTexCoord2f(1,1); glVertex4f( s, -s, 1.0f, w);
				glTexCoord2f(0,1); glVertex4f(-s, -s, 1.0f, w);

			glEnd();

			glEnable(GL_TEXTURE_2D);			
			glBlendFunc(GL_ONE, GL_ZERO);
			glColor4f(1.0f,1.0f,1.0f,1.0f);

			if (phase < PI/2.0)
			{
				x0 = s*cos(phase);
				glBegin(GL_QUADS);
	      
				glTexCoord2f(0, 1); glVertex4f(-s,  s, 1.0f, w);
				glTexCoord2f(1, 1); glVertex4f( s,  s, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s, 0, 1.0f, w);
				glTexCoord2f(0, 0.5f); glVertex4f(-s, 0, 1.0f, w);

				glEnd();

				glBindTexture(GL_TEXTURE_2D, skyBox->moonShadow);

				glBegin(GL_QUADS);
	      
				glTexCoord2f(0, 1); glVertex4f(-s,  x0, 1.0f, w);
				glTexCoord2f(1, 1); glVertex4f( s,  x0, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s, 0, 1.0f, w);
				glTexCoord2f(0, 0.5f); glVertex4f(-s, 0, 1.0f, w);

				glEnd();
			}
			else if (phase < PI)
			{
				x0 = s*cos(phase);

				glBegin(GL_QUADS);
	      
				glTexCoord2f(0, 1); glVertex4f(-s,  s, 1.0f, w);
				glTexCoord2f(1, 1); glVertex4f( s,  s, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s, 0, 1.0f, w);
				glTexCoord2f(0, 0.5f); glVertex4f(-s, 0, 1.0f, w);

				glTexCoord2f(0, 0.5f); glVertex4f(-s,  0, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s,  0, 1.0f, w);
				glTexCoord2f(1, 0); glVertex4f( s, x0, 1.0f, w);
				glTexCoord2f(0, 0); glVertex4f(-s, x0, 1.0f, w);

				glEnd();
			}
			else if (phase < 3.0*PI/2.0)
			{
				x0 = -s*cos(phase);

				glBegin(GL_QUADS);
	      
				glTexCoord2f(0, 1); glVertex4f(-s,  x0, 1.0f, w);
				glTexCoord2f(1, 1); glVertex4f( s,  x0, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s, 0, 1.0f, w);
				glTexCoord2f(0, 0.5f); glVertex4f(-s, 0, 1.0f, w);

				glTexCoord2f(0, 0.5f); glVertex4f(-s,  0, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s,  0, 1.0f, w);
				glTexCoord2f(1, 0); glVertex4f( s, -s, 1.0f, w);
				glTexCoord2f(0, 0); glVertex4f(-s, -s, 1.0f, w);

				glEnd();
			}
			else
			{
				x0 = -s*cos(phase);
				glBegin(GL_QUADS);
	      
				glTexCoord2f(0, 0.5f); glVertex4f(-s,  0, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s,  0, 1.0f, w);
				glTexCoord2f(1, 0); glVertex4f( s, -s, 1.0f, w);
				glTexCoord2f(0, 0); glVertex4f(-s, -s, 1.0f, w);

				glEnd();

				glBindTexture(GL_TEXTURE_2D, skyBox->moonShadow);

				glBegin(GL_QUADS);
	      
				glTexCoord2f(0, 0.5f); glVertex4f(-s,  0, 1.0f, w);
				glTexCoord2f(1, 0.5f); glVertex4f( s,  0, 1.0f, w);
				glTexCoord2f(1, 0); glVertex4f( s, x0, 1.0f, w);
				glTexCoord2f(0, 0); glVertex4f(-s, x0, 1.0f, w);

				glEnd();
			}
			glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
		}
		break;

		case GLOW:
		{
			float s = 0.6f + 0.4f*skyBox->glowOpacity;
			glBegin(GL_QUADS);
      
			glTexCoord2f(0,0);
			glVertex4f(-s,  s, 1.0f, w);
			glTexCoord2f(1,0);
			glVertex4f( s,  s, 1.0f, w);
			glTexCoord2f(1,1);
			glVertex4f( s, -s, 1.0f, w);
			glTexCoord2f(0,1);
			glVertex4f(-s, -s, 1.0f, w);

			glEnd();
		}
		break;
	}
}

void CubeRenderer :: reset()
{
	texturizer.reset();
}