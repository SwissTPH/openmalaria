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


#ifndef MY_GL_WINDOW_H
#define MY_GL_WINDOW_H

#include "Scene.h"
#include "Overlay.h"

class Bridge;

class GL_Window
{
	public:

		GL_Window(Bridge* bridge);

		Scene		scene;

		void		init();
		void		render();
		void		resize(int w, int h);


	private:
    
		int				frameCount;
		Bridge*		bridge;
		Overlay		overlay;
};

#endif
