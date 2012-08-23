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


#ifndef ENVIRONMENT_CONTROLLER_H
#define ENVIRONMENT_CONTROLLER_H

#include "Environment.h"
#include <map>

typedef std::map<float, Environment> EnvironmentMap;
typedef enum {RAINY, FOGGY, SUNNY} WeatherType;
class SkyBox; 

class EnvironmentController
{
	public:

		EnvironmentController(SkyBox* skyBox);

		void setEnvironment(float t);

		void setHumidity(float h)
		{
			humidity = h;
		}

		void setCloudCoverage(float c)
		{
			cloudedness = c;
		}

		float getHumidity()
		{
			return humidity;
		}

		float getCloudCoverage()
		{
			return cloudedness;
		}

		void				addEnvironment( Environment t, WeatherType weather, float time);
		Environment extractEnvironment(float t, const EnvironmentMap& map);

	private:

		SkyBox*					skyBox;
		float					humidity, cloudedness;
		EnvironmentMap	rainy, foggy, sunny;
};

#endif
