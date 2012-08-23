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


#include "SystemTimer.h"

int SystemTimer::initialMsecs = 0;
bool SystemTimer::initialized = false;

int SystemTimer 
:: getMsecs()
{
	if (initialized)
	{		
		#ifdef _WIN32
		//the macro CLK_TCK is obsolete (for unix)
		int t = (1000 * clock() / CLK_TCK);
		#else
                int t = (1000 * clock() / CLOCKS_PER_SEC );
		#endif
		int out = t - initialMsecs;
		initialMsecs = t;
		return out;
	}			
	else
	{
		init();
		return 0;
	}
}

void SystemTimer 
:: init()
{	
	initialMsecs = 
	#ifdef _WIN32
	1000 * clock() / CLK_TCK;
	#else
        1000 * clock() / CLOCKS_PER_SEC;
	#endif
	initialized = true;
}

