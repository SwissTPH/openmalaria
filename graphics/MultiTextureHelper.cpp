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


#include "MultiTextureHelper.h"
#include "SkyBox.h"
#include "Color.h"
#include "macros.h"
#include <cmath>

MultiTextureHelper :: MultiTextureHelper(SkyBox* skyBox)
:	skyBox(skyBox),
	alphaClamp(0.4f) 
{}

void MultiTextureHelper :: reset()
{
	glDisable(GL_TEXTURE_CUBE_MAP);
	glEnable(GL_TEXTURE_2D);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glDisable(GL_ALPHA_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void MultiTextureHelper :: setUp(CubeMode mode)
{
	float adaptation = 0.16f;	

	switch (mode)
	{
		case PRE_RENDERED:
		{
			glEnable(GL_TEXTURE_CUBE_MAP);
			glDisable(GL_TEXTURE_2D);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4f(1,1,1,1);
			glDisable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);
			glDisable(GL_DEPTH_TEST);

			break;
		}

		case AMBIENT:
		{			
			float ambient[4] = COLOR_AS_ARRAY(skyBox->ambientColor);

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->gray);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(ambient);

			glDisable(GL_ALPHA_TEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);

			break;
		}

		case HAZE:
		case HAZE_ON_SKY:
		{			
			float shade[4] = {1,1,1,1};
			if (mode == HAZE_ON_SKY)
			{
				shade[3] = cos(skyBox->sunAngle);
				shade[3] *= shade[3];
				shade[3] *= shade[3];
				shade[3] *= shade[3];
			}
			else skyBox->hazeColor.writeTo(shade);

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->haze);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);
			
			if (mode == HAZE_ON_SKY)
			{
				glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_TRUE);
				glBlendFunc(GL_ONE, GL_ZERO);
			}
			else glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

			break;
		}

		case AFTERGLOW:
		{
			float shade[4] = COLOR_AS_ARRAY(skyBox->afterglowColor);
						
			glEnable(GL_TEXTURE_CUBE_MAP);
			glDisable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->afterglowMask);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);
			glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
      			
			glEnable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);			
			glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
		}
		break;

		case SKY:
		{			
			float sky[4] = COLOR_AS_ARRAY(skyBox->skyColor);

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->sky);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(sky);
                         			
			glEnable(GL_ALPHA_TEST);			
			glAlphaFunc(GL_GREATER, alphaClamp);			
			glEnable(GL_BLEND);
			glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

			break;
		}

		case SUNRISE:
		{    
			/*float sunlight[4] = COLOR_AS_ARRAY_NO_ALPHA(skyBox->currentColor);
      			
			glActiveTexture(GL_TEXTURE0);
			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->sunrise);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);			
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_RGB, GL_CONSTANT);
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_SUBTRACT);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA, GL_TEXTURE);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_ALPHA, GL_SRC_ALPHA);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE1_ALPHA, GL_CONSTANT);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND1_ALPHA, GL_SRC_ALPHA);
			glTexEnvi(GL_TEXTURE_ENV, GL_ALPHA_SCALE, 4);

			glActiveTexture(GL_TEXTURE1);
			glDisable(GL_TEXTURE_CUBE_MAP);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);			
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA, GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_RGB, GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_ALPHA_SCALE, 4);

			glActiveTexture(GL_TEXTURE2);
			glDisable(GL_TEXTURE_CUBE_MAP);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);			
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA, GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_RGB, GL_PREVIOUS);
			glTexEnvi(GL_TEXTURE_ENV, GL_ALPHA_SCALE, 4);		

			glColor4f(0.0f,0.0f,0.0f,skyBox->illuminationAngle/PI - skyBox->angleBias);

			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.0f);
			glBlendFunc(GL_ONE, GL_ZERO);	

			break;										*/
			float sunlight[4] = COLOR_AS_ARRAY_NO_ALPHA(skyBox->currentColor);

			glEnable(GL_TEXTURE_CUBE_MAP);

			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->sunrise);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(sunlight);

			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_LESS, skyBox->illuminationAngle/PI + skyBox->angleBias);
			glBlendFunc(GL_ONE, GL_ZERO);

			break;
		}

		case SUNSET:
		{ 
			float shade[4] = COLOR_AS_ARRAY_NO_ALPHA(skyBox->currentColor);				

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->sunset);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_LEQUAL, skyBox->illuminationAngle/PI + skyBox->angleBias);
			glBlendFunc(GL_ONE, GL_ZERO);
		}
		break;

		case NORMALS:
		{			
			/*glActiveTexture(GL_TEXTURE1);
			glDisable(GL_TEXTURE_CUBE_MAP);

			glActiveTexture(GL_TEXTURE2);
			glDisable(GL_TEXTURE_CUBE_MAP);	*/

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->normal);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_COMBINE);
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_RGB, GL_DOT3_RGB);
			glTexEnvi(GL_TEXTURE_ENV, GL_COMBINE_ALPHA, GL_REPLACE);
			glTexEnvi(GL_TEXTURE_ENV, GL_SOURCE0_ALPHA, GL_TEXTURE);
			glTexEnvi(GL_TEXTURE_ENV, GL_OPERAND0_ALPHA, GL_SRC_ALPHA); 			

			glEnable(GL_BLEND);		
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, alphaClamp);
			glBlendFunc(GL_DST_COLOR, GL_ZERO);
		}
		break;

		case DIFFUSE:
		{
			float shade[4] = {1,1,1,1};

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->diffuse);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);		
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, alphaClamp);
			glBlendFunc(GL_DST_COLOR, GL_SRC_COLOR);
		}
		break;

		case STARS				:
		{
			float shade[4] = COLOR_AS_ARRAY(skyBox->currentColor);

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->stars);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);		
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.0f);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		}
		break;

		case MOON:
		{ 
			float shade[4] = COLOR_AS_ARRAY(skyBox->currentColor);			

			glDisable(GL_TEXTURE_CUBE_MAP);
			glEnable(GL_TEXTURE_2D);

			glBindTexture(GL_TEXTURE_2D, skyBox->moon);
			
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.02f);
			glBlendFunc(GL_DST_ALPHA, GL_ONE);
		}
		break;

		case MOON_SHADOW:
		{ 
			glDisable(GL_TEXTURE_CUBE_MAP);
			glEnable(GL_TEXTURE_2D);
      			
			glBindTexture(GL_TEXTURE_2D, skyBox->moonMask);     
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

			glEnable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);			
		}
		break;

		case SUN:
		{ 
			float shade0[4] = {1,1,1,1};

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->diffuse);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade0);

			glEnable(GL_BLEND);		
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, alphaClamp);
			glBlendFunc(GL_DST_COLOR, GL_SRC_COLOR);
			
			float shade[4] = COLOR_AS_ARRAY(skyBox->sunColor);		
			
			glDisable(GL_TEXTURE_CUBE_MAP);
			glEnable(GL_TEXTURE_2D);   
			glBindTexture(GL_TEXTURE_2D, skyBox->sun);   
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, 0.02f);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		}
		break;

		case GLOW:
		{ 
			float shade0[4] = {1,1,1,1};

			glEnable(GL_TEXTURE_CUBE_MAP);
			glBindTexture(GL_TEXTURE_CUBE_MAP, skyBox->diffuse);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade0);

			glEnable(GL_BLEND);		
			glEnable(GL_ALPHA_TEST);
			glAlphaFunc(GL_GREATER, alphaClamp);
			glBlendFunc(GL_DST_COLOR, GL_SRC_COLOR);

			float shade[4] = COLOR_AS_ARRAY(skyBox->currentColor);			
      			
			glDisable(GL_TEXTURE_CUBE_MAP);
			glEnable(GL_TEXTURE_2D);
			glBindTexture(GL_TEXTURE_2D, skyBox->glow);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glColor4fv(shade);

			glEnable(GL_BLEND);
			glDisable(GL_ALPHA_TEST);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		}
		break;

	}
}