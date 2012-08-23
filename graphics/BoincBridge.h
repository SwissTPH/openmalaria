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


#ifndef BOINC_BRIDGE
#define BOINC_BRIDGE

#include "Bridge.h"
#include "int2.h"
#include "Key.h"
#include "KeyHandler.h"
#include "MouseHandler.h"
#include <map>

class BoincTranslator;

class BoincBridge
{
	public:

		BoincBridge();

		static bool		initialized;

		GL_Window*		getWindow(int w, int h, WindowMode mode = FULLSCREEN);
		void					initialize();
		void					render();
		void					pollEvents();


	private:

		void					buildKeyMap();

		GL_Window*		window;
		int2					windowSize, windowLocation;    		
		int2					absolute, relative; /*	mouse coordinates */ 

		MouseHandler	mouseHandler;
		KeyHandler		keyHandler;

		SDL_Translator* translator;
};

#endif