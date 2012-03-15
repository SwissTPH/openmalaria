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


#include "TextureLoader.h"
#include "boinc_api.h"

//		loadTexture2D
GLuint TextureLoader :: loadTexture2D(std::string filename, TextureType type, 
																			GLint edge)
{
	GLuint index = generateAndBind(GL_TEXTURE_2D);
	loadImage(filename, GL_TEXTURE_2D, type);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);

	return index;
}

//		loadCubeMap
GLuint TextureLoader :: loadCubeMap(std::string prefix, TextureType type, Emptiness emptiness)
{	const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	GLuint index = generateAndBind(GL_TEXTURE_CUBE_MAP);

	if (emptiness == ALL_EQUAL)
	{
		loadImage(prefix, GL_TEXTURE_CUBE_MAP_POSITIVE_X, type);
		loadImage(prefix, GL_TEXTURE_CUBE_MAP_NEGATIVE_X, type);
		loadImage(prefix, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, type);
		loadImage(prefix, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, type);
		loadImage(prefix, GL_TEXTURE_CUBE_MAP_POSITIVE_Z, type);
		loadImage(prefix, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, type);
	}
	else
	{
		retval = boinc_resolve_filename((prefix + "_east.png").c_str(),imagefile,MAX_LENGTH);
		loadImage(imagefile, GL_TEXTURE_CUBE_MAP_POSITIVE_X, type);
		retval = boinc_resolve_filename((prefix + "_west.png").c_str(),imagefile,MAX_LENGTH);
		loadImage(imagefile, GL_TEXTURE_CUBE_MAP_NEGATIVE_X, type);
		switch (emptiness)
		{
			case EMPTY:
					retval = boinc_resolve_filename("blank.png",imagefile,MAX_LENGTH);
					loadImage(imagefile, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, type);
				break;

			case FULL:
					retval = boinc_resolve_filename("full.png",imagefile,MAX_LENGTH);
					loadImage(imagefile, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, type);
				break;

			case DOME:	
					retval = boinc_resolve_filename((prefix + "_top.png").c_str(),imagefile,MAX_LENGTH);
					loadImage(imagefile, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, type);
				break;
		}
		retval = boinc_resolve_filename((prefix + "_floor.png").c_str(),imagefile,MAX_LENGTH);
		loadImage(imagefile, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, type);
		retval = boinc_resolve_filename((prefix + "_north.png").c_str(),imagefile,MAX_LENGTH);
		loadImage(imagefile, GL_TEXTURE_CUBE_MAP_POSITIVE_Z, type);
		retval = boinc_resolve_filename((prefix + "_south.png").c_str(),imagefile,MAX_LENGTH);
		loadImage(imagefile, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, type);
	}

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	return index;
}

//		generateAndBind
GLuint TextureLoader :: generateAndBind(GLenum target)
{
	GLuint index;
	glGenTextures(1, &index);
	glBindTexture(target, index);
	return index;
}

//		loadImage
void TextureLoader :: loadImage(std::string filename, GLenum target, TextureType type)
{
	LOG("loading '" << filename << "'")
	ILuint image;

  ilGenImages(1, &image);
	ilBindImage(image);

 	if (!ilLoadImage(ILstring(filename.c_str())))	
	{
		CRASH( "unable to load " << filename << "!" )	
	}
	
	switch (type)
	{
		case RGBA_TEXTURE: 
		{
			ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);	
			ILubyte* data = ilGetData();

			glTexImage2D(target, 0, GL_RGBA, ilGetInteger(IL_IMAGE_WIDTH), 
					ilGetInteger(IL_IMAGE_HEIGHT), 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
		}
		break;
		case RGB_TEXTURE: 
		{
			ilConvertImage(IL_RGB, IL_UNSIGNED_BYTE);	
			ILubyte* data = ilGetData();

			glTexImage2D(target, 0, GL_RGB, ilGetInteger(IL_IMAGE_WIDTH), 
					ilGetInteger(IL_IMAGE_HEIGHT), 0, GL_RGB, GL_UNSIGNED_BYTE, data);
		}
		break;
		case GRAYSCALE_TEXTURE:
		{ 			
			int H = ilGetInteger(IL_IMAGE_HEIGHT);
			ILubyte	*data = ilGetData();

			glTexImage2D(target, 0, GL_ALPHA8, H, H, 0, GL_ALPHA, GL_UNSIGNED_BYTE, data);
		}
		break;
		case DESATURATED_TEXTURE:
		{ 			
			unsigned int H = ilGetInteger(IL_IMAGE_HEIGHT);
			ilConvertImage(IL_RGBA, IL_UNSIGNED_BYTE);
			ILubyte	*data = ilGetData();
			ILubyte *final = new ILubyte[H*H];

			double r,g,b, lum;
			ITERATE(i,H*H)
			{
				r = (double)data[4*i];
				g = (double)data[4*i + 1];
				b = (double)data[4*i + 2];
				lum = 0.299*r + 0.7512*g + 0.114*b;
				CLAMP(lum,0.0,255.0)
				final[i] = ((int)lum)&0xff;
			}
      			
			glTexImage2D(target, 0, GL_ALPHA8, H, H, 0, GL_ALPHA, GL_UNSIGNED_BYTE, final);
			delete final;
		}
		break;
	}

	ilDeleteImages(1, &image);
}

//		generateCubeMap
GLuint TextureLoader :: generateCubeMap(ProceduralMapType type, unsigned int size)
{
	GLuint index = generateAndBind(GL_TEXTURE_CUBE_MAP);

	generateImage(type, GL_TEXTURE_CUBE_MAP_POSITIVE_X, size);
	generateImage(type, GL_TEXTURE_CUBE_MAP_NEGATIVE_X, size);
	generateImage(type, GL_TEXTURE_CUBE_MAP_POSITIVE_Y, size);
	generateImage(type, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y, size);
	generateImage(type, GL_TEXTURE_CUBE_MAP_POSITIVE_Z, size);
	generateImage(type, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z, size);

	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	return index;
}

//		generateImage
void TextureLoader :: generateImage(ProceduralMapType type, GLenum target, unsigned int size)
{
	switch (type)
	{
		case EMPTY_RGB_MAP :
		{
			GLubyte* data = new GLubyte[size*size*3];
			ITERATE(i,size*size*3)
			{
				data[i] = 0;
			}
			glTexImage2D(target, 0, GL_RGB, size, size, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
			delete[] data;

		}
		break;

		case Z_VALUE_MAP :
		{
			GLubyte* data = new GLubyte[size*size];

			switch(target)
			{
				case GL_TEXTURE_CUBE_MAP_POSITIVE_X:				
				{
					ITERATE(i, size*size)
					{
						int x = (i%size) - size/2, y = (i/size) - size/2;
						double v = 128 - 127.0*(double)x/sqrt((double)(size*size/4 + x*x + y*y));
						CLAMP( v, 0.0, 255.0 )
						data[i] = (GLubyte)(int)v;						
					}
					break;
				}
				case GL_TEXTURE_CUBE_MAP_NEGATIVE_X:				
				{
					ITERATE(i, size*size)
					{
						int x = (i%size) - size/2, y = (i/size) - size/2;
						double v = 128 + 127.0*(double)x/sqrt((double)(size*size/4 + x*x + y*y));
						CLAMP( v, 0.0, 255.0 )
						data[i] = (GLubyte)(int)v;						
					}
					break;
				}
				case GL_TEXTURE_CUBE_MAP_POSITIVE_Y:
				{
					ITERATE(i, size*size)
					{
						int x = (i%size) - size/2, y = (i/size) - size/2;
						double v = 128 + 127.0*(double)y/sqrt((double)(size*size/4 + x*x + y*y));
						CLAMP( v, 0.0, 255.0 )
						data[i] = (GLubyte)(int)v;						
					}
					break;
				}
				case GL_TEXTURE_CUBE_MAP_NEGATIVE_Y:				
				{
					ITERATE(i, size*size)
					{
						int x = (i%size) - size/2, y = (i/size) - size/2;
						double v = 128 - 127.0*(double)y/sqrt((double)(size*size/4 + x*x + y*y));
						CLAMP( v, 0.0, 255.0 )
						data[i] = (GLubyte)(int)v;						
					}
					break;
				}
				case GL_TEXTURE_CUBE_MAP_POSITIVE_Z:				
				{
					ITERATE(i, size*size)
					{
						int x = (i%size) - size/2, y = (i/size) - size/2;
						double v = 128 + 127.0*0.5*(double)size/sqrt((double)(size*size/4 + x*x + y*y));
						CLAMP( v, 0.0, 255.0 )
						data[i] = (GLubyte)(int)v;
					}
					break;
				}
				case GL_TEXTURE_CUBE_MAP_NEGATIVE_Z:				
				{
					ITERATE(i, size*size)
					{
						int x = (i%size) - size/2, y = (i/size) - size/2;
						double v = 128 - 127.0*0.5*(double)size/sqrt((double)(size*size/4 + x*x + y*y));
						CLAMP( v, 0.0, 255.0 )
						data[i] = (GLubyte)(int)v;
					}
					break;
				}
			}

			glTexImage2D(target, 0, GL_INTENSITY8, size, size, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, data);

			break;
		}
	}
}
