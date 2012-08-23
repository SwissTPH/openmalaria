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


#ifndef FIELD_DISPLAY_H
#define FIELD_DISPLAY_H

#include "ChartRenderer.h"
#include "Display.h"
#include "TextureLoader.h"
#include "Line.h"

class SkyBox;
class Scene;

class FieldDisplay
{
	friend class ChartRenderer;

	public:

		FieldDisplay(DisplayMM* display, int maxSampleCount, int sampleSize, float3 dim);
    
		bool	soft;

		void	render();
		void	update(float deltaT);
		void	addData(Sample newData);	// USED BY OTHER THREAD

	private:

		bool					dataReady, dataRead;
		int						freshness;
		float					dataOffset, dataThroughput,
									timeSinceData, width, depth, height,
									arrowWidth, arrowThickness;
		float3				timeTip, ageTip, infectTip;
		DisplayMM*			display;
		SkyBox*				skyBox;
		Scene*				scene;
		Line					*time, *age, *infectiousness;
		TextureLoader textureLoader;
		ChartRenderer	chartRenderer;
		unsigned int	maxSampleCount, newSamples, sampleSize,
									arrowTex;									
		SampleList		data0, data1, data;
		SampleList		*frontBuffer, *backBuffer;

		void	renderArrows(double phi);
};

#endif

