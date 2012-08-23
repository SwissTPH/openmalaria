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


#include "Surface.h"
#include "SurfaceProvider.h"
#include "gl_headers.h"
#include "Line.h"
#include "Font.h"
#include "macros.h"

//		c'tor
Line :: Line(Surface* parent, int h)
:	parent(parent),
	provider(parent->provider),
	textureYOffset(h),
	widthInChars(0),
	texture(parent->getTexture()),
	font(parent->font),
	textureSize(parent->size),
	texCoordWidth(0),
	texCoordHeight((float)parent->font->charSize.y/(float)parent->size),
	texCoordOffsetY((float)h/(float)parent->size),
	cursor(0),
	maxCharWidth(parent->size/parent->font->charSize.x)
{
}

//		d'tor
Line :: ~Line()
{
	setUsed(false);
}

//		clear
void Line :: clear()
{
	widthInChars = 0;
	cursor = 0;
}

//		print	String
void Line :: print(String s)
{
	glBindTexture(GL_TEXTURE_2D, texture);

	int t;
	if (cursor + (int)s.length() >= maxCharWidth) t = maxCharWidth;
	else t = cursor + (int)s.length();

	for (int i = cursor; i < t; i++)
	{  		
			glTexSubImage2D( GL_TEXTURE_2D, 0, font->charSize.x*i, textureYOffset, 
				font->charSize.x, font->charSize.y, GL_RGBA, GL_UNSIGNED_BYTE, font->data[s[i - cursor]]);		
	}
	widthInChars = (float) t;
	cursor = t;
	texCoordWidth = widthInChars*font->charSize.x/(float)textureSize;
}

//		print int
void Line :: print(int i)
{
	if (!i) 
	{
		print("0");
		return;
	}

	int i_pos = i > 0 ? i : -i;
	int digits = (int)(log10((double)i_pos)) + (i > 0 ? 1 : 2);
	
	char* c = new char[digits + 1];
	c[digits] = 0;
	if (i < 0)
	{
		c[0] = '-';	
	}

	int k = digits - 1;
	while (i_pos > 0)
	{
		c[k] = '0' + i_pos%10;
		i_pos /= 10;
		k--;
	}	
	
	print(String(c));  	

	delete c;
}

//		setUsed
void Line :: setUsed(bool whether)
{
	if (whether) 
	{
		provider->usedLines.insert(provider->usedLines.end(),this);
		parent->linesInUse++;
	}
	else 
	{
		provider->usedLines.remove(this);
		provider->freeLines.insert(provider->usedLines.end(),this);
		parent->linesInUse--;
		if (!parent->linesInUse)
		{	
			parent->unload();
			delete parent;
		}
	}
}
//		changeFont
void Line :: changeFont(FontMM* font)
{
	this->font = font;
}

//		render
void Line :: render(float2 charSize, float2 alignment)
{
	float		w = widthInChars*charSize.x,
					h = charSize.y;

	glBindTexture(GL_TEXTURE_2D, texture);
	glBegin(GL_QUADS);

		glTexCoord2f(texCoordWidth, texCoordOffsetY);
		glVertex3f(	(1.0f - alignment.x)*w, alignment.y*h, 0);		

		glTexCoord2f(0, texCoordOffsetY);
		glVertex3f( -alignment.x*w, alignment.y*h, 0);

		glTexCoord2f(0, texCoordOffsetY + texCoordHeight);
		glVertex3f( -alignment.x*w, (alignment.y - 1.0f)*h, 0);

		glTexCoord2f(texCoordWidth, texCoordOffsetY + texCoordHeight);
		glVertex3f(	(1.0f - alignment.x)*w, (alignment.y - 1.0f)*h, 0);

	glEnd();
}

//		render
void Line :: render(float2 charSize, float2 alignment, 
										Color top, Color bottom)
{
	float		w = widthInChars*charSize.x,
					h = charSize.y;

	glBindTexture(GL_TEXTURE_2D, texture);
	glBegin(GL_QUADS);

		top.set();

		glTexCoord2f(texCoordWidth, texCoordOffsetY);
		glVertex3f(	(1.0f - alignment.x)*w, alignment.y*h, 0);		

		glTexCoord2f(0, texCoordOffsetY);
		glVertex3f( -alignment.x*w, alignment.y*h, 0);

		bottom.set();

		glTexCoord2f(0, texCoordOffsetY + texCoordHeight);
		glVertex3f( -alignment.x*w, (alignment.y - 1.0f)*h, 0);

		glTexCoord2f(texCoordWidth, texCoordOffsetY + texCoordHeight);
		glVertex3f(	(1.0f - alignment.x)*w, (alignment.y - 1.0f)*h, 0);

	glEnd();
}
