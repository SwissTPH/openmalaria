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


#include "GL_Window.h"
#include "GraphicsBridge.h"
#include "macros.h"
#include "boinc_api.h" 
//#include "graphics_api.h"
#include "diagnostics.h"
#include "util.h"
#include "filesys.h"

#include "Key.h"
#include "KeyHandler.h"
#include "MouseHandler.h"
#include "config.h"

#include <iostream>

#if (defined(_GRAPHICS_6))
#ifdef _WIN32
#include "boinc_win.h"
#else
#include <math.h>
#include <execinfo.h>
#endif

#include "ShmStruct.h"
#include "graphics2.h"
#include "boinc_gl.h"
#include "app_ipc.h"
UC_SHMEM* shmem = NULL;
//#ifdef _WIN32
struct APP_INIT_DATA dataBOINC;
//#endif
void loadDataFromShm()
{
	if (shmem) {
		add_and_copy_data(shmem->KappaArray);
	}
}

#endif


using namespace std;
class Scene;

#if (!(defined(_GRAPHICS_6)))
extern void worker();
#endif

void app_graphics_render(int xs, int ys, double time_of_day) 
{
#if (defined(_GRAPHICS_6))
	//boinc_graphics_get_shmem() must be called after
	// boinc_parse_init_data_file()
	// Put this in the main loop to allow retries if the
	// worker application has not yet created shared memory
if (shmem == NULL) {
	shmem = (UC_SHMEM*)boinc_graphics_get_shmem("malariacontrol");
	}
	loadDataFromShm();
#endif
	GraphicsBridge::window -> render();
}

#if (!defined(_NO_GRAPHICS))
void AppInvalidParameterHandler(const wchar_t* expression, const wchar_t* function, const wchar_t* file, unsigned int line,	uintptr_t pReserved ) {
	fprintf(
		stderr,
		"Invalid parameter detected in function %s. File: %s Line: %d\n",
		function,
		file,
		line
	);
	fprintf(
		stderr,
		"Expression: %s\n",
		expression
	);
}
int main(int argc, char *argv[])
{   
#if (defined(_GRAPHICS_6))
	int retval;  
	retval = boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED|BOINC_DIAG_REDIRECTSTDERR);
	//Where is this defined?
	#ifdef _WIN32
	_set_invalid_parameter_handler(AppInvalidParameterHandler);
	#endif
	boinc_parse_init_data_file();
	boinc_get_init_data(dataBOINC);
	boinc_graphics_loop(argc, argv);	//GRAPHICS2_WIN.C or GRAPHICS2_unix.C. header: graphics2.h
	boinc_finish_diag();
#else
	int retval;  
	retval = boinc_init_diagnostics(BOINC_DIAG_DUMPCALLSTACKENABLED|BOINC_DIAG_REDIRECTSTDERR);
	// Every once and awhile something looks at a std::vector or some other
	// CRT/STL construct that throws an exception when an invalid parameter
	// is detected.  In this case we should dump whatever information we 
	// can and then bail.  When we bail we should dump as much data as 
	// possible.
	_set_invalid_parameter_handler(AppInvalidParameterHandler);
	//boinc_init_graphics(worker);	//WINDOWS_OPENGL.C

	retval = boinc_finish(0);	
	if (retval){
		exit(retval);
	}

	//Unreachable
#endif
	return 0;
}
#endif
