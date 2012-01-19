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


#include "Overlay.h"
#include "Scene.h"
#include "SkyBox.h"
#include "GraphicsBridge.h"
#include "Font.h"
#include "boinc_api.h"
#include "config.h"


#if (defined(_GRAPHICS_6))
#include "ShmStruct.h"
extern UC_SHMEM* shmem;
extern APP_INIT_DATA dataBOINC;
#else
extern float fdone;
APP_INIT_DATA dataBOINC;
#endif


Overlay :: Overlay()
{
    const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	retval = boinc_resolve_filename("font_outlined.png",imagefile,MAX_LENGTH);
	nameFont = new FontMM(imagefile, int2(22,32), int2(32,32), int2(5,0));
	name = SurfaceProvider :: getInstance() -> getLine();	
	name -> changeFont(nameFont);
	//name -> print("Test Name");
	name -> print(dataBOINC.user_name);
	done = SurfaceProvider :: getInstance() -> getLine();
#if (defined(_GRAPHICS_6))
	//done -> print((float)(shmem->fraction_done));
#else
	done -> print(fdone);
#endif
	sezCredits = SurfaceProvider :: getInstance() -> getLine();	
	sezCredits -> print("credits:");
	credits = SurfaceProvider :: getInstance() -> getLine();	
	credits -> print(dataBOINC.user_total_credit);
	TextureLoader tl;
	
	retval = boinc_resolve_filename("scrollBarInside.png",imagefile,MAX_LENGTH);
	unsigned int in = tl.loadTexture2D(imagefile, GRAYSCALE_TEXTURE);
	
	retval = boinc_resolve_filename("scrollBarOutside.png",imagefile,MAX_LENGTH);
	unsigned int out = tl.loadTexture2D(imagefile, GRAYSCALE_TEXTURE);

	progressBar = new ProgressBar(in, out);
}

Overlay :: ~Overlay()
{
	delete name;
}

void Overlay :: render()
{
	glDisable(GL_DEPTH_TEST);
	float opacity = scene->overlayPresence;
	float2	smallPrint = float2(0.075f, 0.09f),
					fatPrint = float2(0.17f, 0.21f);

	glLoadIdentity();
	float g = 0.1f;
	Color		sun = opacity*(scene->skyBox->sunlightColor + Color(g,g,g)),
					shade = opacity*(scene->skyBox->ambientColor + Color(g,g,g));
	shade.a = 1.0f;
	sun.a = 1.0f;
	glTranslatef(0,  1.0f, -0.5f*opacity - 0.5f);	
	name->render(fatPrint, float2(0.888f, 0.93f), shade, sun);
	glTranslatef(1.9f,0,0);
	sezCredits->render(smallPrint, float2(1.2f, 0.0f));
	credits->render(smallPrint, float2(0.0f, 0.0f));

	glLoadIdentity();
	done->clear();
	done->print((int)(progressBar->value*100.0f));
	done->print("%");
	glTranslatef(0, -1.07f, -0.5f*opacity - 0.5f);
	progressBar->render(shade, sun);
	done->render(smallPrint, float2(0.5f,1.2f));
	glEnable(GL_DEPTH_TEST);
}

