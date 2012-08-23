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
#include "Key.h"
#include "KeyHandler.h"
#include "MouseHandler.h"
#include "graphics2.h"
#include "boinc_gl.h"
#include "GL_Window.h"
#include "ShmStruct.h"
#include "config.h"

void boinc_app_key_press(int key, int modifier) 
{
	GraphicsBridge::keyHandler -> keyPressed( GraphicsBridge :: translate(key) );
}

void boinc_app_key_release(int key, int modifier) 
{
	GraphicsBridge::keyHandler -> keyReleased( GraphicsBridge :: translate(key) );
}

void boinc_app_mouse_button(int x, int y, int which, int is_down) 
{ 	
	GraphicsBridge :: mouseHandler -> mouseButtonUsed( which, is_down );
}

void boinc_app_mouse_move(int x, int y, int left, int middle, int right) 
{
	GraphicsBridge :: mouseMoved(x, y);
}

void app_graphics_reread_prefs(){}

void app_graphics_resize(int w, int h) 
{
	GraphicsBridge::window -> resize(w,h);
}

void app_graphics_init() 
{
#if (defined(_GRAPHICS_6))
	GraphicsBridge::init(KappaArraySize);
#else
	GraphicsBridge::init(KappaArraySize);
#endif
}



