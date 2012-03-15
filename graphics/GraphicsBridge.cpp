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


#include "GraphicsBridge.h"
#include "GL_Window.h"
#include "IL/il.h"
#include "IL/ilu.h"
#include "IL/ilut.h"
#include "SurfaceProvider.h"
#include "Font.h"
#include "Anopheles.h"
#include "GL_Window.h"
#include "MouseHandler.h"
#include "KeyHandler.h"
#include "Key.h"
#include <iostream>
#include "boinc_api.h"
#include "Bridge.h"

FieldDisplay*	GraphicsBridge	::	display = 0;
int						GraphicsBridge	::	sampleSize = 0;
int						GraphicsBridge	::	width = -1;
int						GraphicsBridge	::	height = -1;
int						GraphicsBridge	::	preRenderedBoxResolution = -1;
int2					GraphicsBridge	::	mouse = int2(0,0);
String				GraphicsBridge	::	imagePath = "images/";
ProgressBar*	GraphicsBridge	::	progressBar = 0;
GL_Window*		GraphicsBridge	::	window = 0;
MouseHandler* GraphicsBridge	::	mouseHandler = 0;
KeyHandler*		GraphicsBridge	::	keyHandler = 0;
GraphicsBridge	::	KeyMap				GraphicsBridge	::	keyMap;
 
void GraphicsBridge	::	mouseMoved( int x, int y )
{
	int2 newMouse = int2(x,y);
	int2 delta = newMouse - mouse;
	mouse = newMouse;
	mouseHandler->mouseMoved(mouse, delta);
}

Key GraphicsBridge	::	translate( int key )
{ 	
	Key k = Key();
	if ( key > 64 && key < 91 )
	{             		
		k.isCharacter = true;
		k.character = key;
	}
	else
	{
		KeyMap::iterator i = keyMap.find(key);
		k.isCharacter = false;		
		k.specialKey = i->second;		
	}
	return k;
}

void GraphicsBridge	::	init( int bands )
{
	const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	
	display = 0;
	sampleSize = bands;
	imagePath = "images/";

	keyMap.clear();

	keyMap[32] = SPACE;

	keyMap[112] = F1;
	keyMap[113] = F2;
	keyMap[114] = F3;
	keyMap[115] = F4;

	keyMap[116] = F5;
	keyMap[117] = F6;
	keyMap[118] = F7;
	keyMap[119] = F8;

	keyMap[120] = F9;
	keyMap[121] = F10;
	keyMap[122] = F11;
	keyMap[123] = F12;

  
	ilInit();
  iluInit();
	ilutInit();
	retval = boinc_resolve_filename("font_nominal.png",imagefile,MAX_LENGTH);
	FontMM* font = new FontMM(imagefile, int2(22,32), int2(32,32), int2(5,0));	
	SurfaceProvider::init(512, font);

	window = new GL_Window(0);
	//window = new GL_Window(NULL);
	
	mouseHandler = new MouseHandler();
	keyHandler = new KeyHandler();

	mouseHandler->controller = &(window->scene.controller);
	keyHandler->addListener(&(window->scene.controller));
	
	Anopheles::init(); 
}

void add_data(float* data)
{
	if (GraphicsBridge :: display) GraphicsBridge :: display -> addData( data );
}

void add_and_copy_data(float* data)
{
	if (GraphicsBridge :: display) 
	{
		int size = GraphicsBridge :: sampleSize;
		float* copy = new float[size];
		for (int i = 0; i < size; i++) copy[i] = data[i];
		GraphicsBridge :: display -> addData( copy );
	}
}
