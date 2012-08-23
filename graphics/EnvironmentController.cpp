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


#include "EnvironmentController.h"
#include "SkyBox.h"

EnvironmentController :: EnvironmentController(SkyBox* skyBox)
:		skyBox(skyBox), humidity(0.0f), cloudedness(0.0f)
{

		/*	Environment constructor arguments:

							sunlight, sky,
							sun, ambient,
							shadow, haze
		*/

			Environment dawn = Environment(
				Color(0.25f, 0.17f, 0.00f, 1.0f),		Color(0.4f,0.0f,0.3f,1.0f),
				Color(1,0.6f,0.2f,1),								Color(0.1f,0.1f,0.18f,1.0f),
				Color(0.02f,0.02f,0.02f,1.0f),			Color(0.1f,0.14f,0.24f,0.7f)	);

			Environment morning = Environment(
				Color(0.58f,0.53f,0.41f,1.0f),			Color(0.2f,0.3f,0.6f,1.0f),
				Color(1,1,0.4f,1),									Color(0.23f,0.27f,0.32f,0.5f),
				Color(0.05f,0.07f,0.06f,1.0f),			Color(0.9f,0.9f,0.6f,0.25f)	);

			Environment foggyMorning = Environment(
				Color(0.78f,0.73f,0.68f,1.0f),			Color(0.2f,0.4f,0.8f,1.0f),
				Color(1,1,0.4f,1),										Color(0.23f,0.27f,0.29f,1.0f),
				Color(0.05f,0.07f,0.06f,1.0f),			Color(0.4f,0.5f,0.8f,0.35f)	);

			Environment lateMorning = Environment(
				Color(0.93f,0.9f,0.8f,1.0f),				Color(0.2f,0.4f,0.8f,1.0f),
				Color(1,1,0.4f,1),										Color(0.33f,0.33f,0.33f,1.0f),
				Color(0.03f,0.03f,0.03f,1.0f),			Color(0.2f,0.4f,1,0.3f)	);

			Environment noon = Environment(
				Color(1.0f,1.0f,1.0f,0.4f),				Color(0.1f,0.24f,0.6f,1.0f),
				Color(1.0f, 0.9f, 0.7f, 1.0f),		Color(0.3f,0.3f,0.3f,1.0f),
				Color(0.0f,0.0f,0.0f,1.0f),				Color(0.6f,0.8f,1,0.0f)	);

			Environment afternoon = Environment(
				Color(0.7f,0.7f,0.7f,1.0f), Color(0.1f,0.2f,0.8f,1.0f),
				Color(1.0f, 0.9f, 0.7f, 1.0f), Color(0.2f,0.21f,0.23f,1.0f),
				Color(0.08f,0.1f,0.09f,1.0f), Color(0.6f,0.8f,1,0.1f)	);

			Environment dusk = Environment(
				Color(0.8f, 0.65f, 0.2f, 1.0f), Color(0.0f, 0.0f, 0.4f, 1.0f),
				Color(1.0f, 0.7f, 0.4f, 1.0f), Color(0.2f, 0.15f, 0.15f, 1.0f),
				Color(0.0f,0.0f,0.0f,1.0f), Color(0.95f, 0.4f, 0.05f, 0.25f)	);		
			
			Environment foggyDusk = Environment(
				Color(0.6f, 0.65f, 0.4f, 1.0f), Color(0.2f, 0.2f, 0.2f, 1.0f),
				Color(0.9f, 0.9f, 0.8f, 1.0f), Color(0.15f, 0.15f, 0.18f, 0.3f),
				Color(0.2f,0.2f,0.25f,1.0f), Color(0.4f, 0.4f, 0.0f, 0.6f)	);

			Environment sunset = Environment(				
				Color(0.6f, 0.25f, 0.1f, 1.0f), Color(0.0f, 0.0f, 0.0f, 1.0f),
				Color(1.0f, 0.75f, 0.5f, 1.0f), Color(0.3f, 0.25f, 0.25f, 1.0f),
				Color(0.0f,0.0f,0.0f,1.0f), Color(0.85f, 0.4f, 0.05f, 0.65f)	);

 			Environment fog = Environment(
				Color(0.5f, 0.5f, 0.5f, 1.0f), Color(0.4f,0.4f,0.44f,1.0f),
				Color(0.5f, 0.5f, 0.5f, 1.0f), Color(0.3f,0.3f,0.3f,1.0f),
				Color(0.4f,0.4f,0.4f,1.0f), Color(0.6f,0.8f,1.9f,1.0f)	);
			
			Environment rain = Environment(
				Color(0.29f, 0.27f, 0.35f, 1.0f), Color(0.14f,0.1f,0.17f,1.0f),
				Color(0,0,0,0), Color(0.22f,0.2f,0.25f,0.2f),
				Color(0.22f,0.2f,0.25f,1.0f), Color(0.3f,0.3f,0.4f,0.98f)	);

			Environment night = Environment(
				Color(0.25f, 0.25f, 0.32f, 1.0f), Color(0.0f,0.0f,0.0f,1.0f),
				Color(0.7f,0.7f,0.82f,1.0f), Color(0.2f,0.2f,0.28f,1.0f),
				Color(0.02f,0.02f,0.02f,1.0f), Color(0.1f,0.14f,0.24f,0.3f)	);

			Environment foggyNight = Environment(
				Color(0.15f, 0.15f, 0.22f, 1.0f), Color(0.0f,0.0f,0.0f,1.0f),
				Color(0.1f,0.1f,0.116f,1.0f), Color(0.13f,0.13f,0.2f,1.0f),
				Color(0.02f,0.02f,0.02f,1.0f), Color(0.1f,0.12f,0.26f,0.9f)	);

			Environment rainyNight = Environment(
				Color(0.05f, 0.05f, 0.11f, 1.0f), Color(0.0f,0.0f,0.0f,1.0f),
				Color(0.0f,0.0f,0.0f,1.0f), Color(0.14f,0.14f,0.23f,1.0f),
				Color(0.02f,0.02f,0.06f,1.0f), Color(0.18f,0.18f,0.2f,0.9f)	);

	addEnvironment(night, SUNNY, -1000.0f);
	addEnvironment(night, SUNNY, -0.1f);
	addEnvironment(dawn, SUNNY, 0.0f);
	addEnvironment(morning, SUNNY, 0.1f);
	addEnvironment(lateMorning, SUNNY, 0.4f);
	addEnvironment(noon, SUNNY, 0.55f);
	addEnvironment(afternoon, SUNNY, 0.85f);
	addEnvironment(dusk, SUNNY, 1.0f);
	addEnvironment(night, SUNNY, 1.06f);
	addEnvironment(night, SUNNY, 1000.0f);

	addEnvironment(foggyNight, FOGGY, -1000.0f);
	addEnvironment(foggyNight, FOGGY, 0.0f);
	addEnvironment(foggyMorning, FOGGY, 0.05f);
	addEnvironment(fog, FOGGY, 0.15f);
	addEnvironment(fog, FOGGY, 0.75f);
	addEnvironment(foggyDusk, FOGGY, 1.0f);
	addEnvironment(foggyNight, FOGGY, 1.06f);
	addEnvironment(foggyNight, FOGGY, 1000.0f);

	addEnvironment(rainyNight, RAINY, -1000.0f);
	addEnvironment(rainyNight, RAINY, 0.08f);
	addEnvironment(rain, RAINY, 0.12f);
	addEnvironment(rain, RAINY, 0.88f);
	addEnvironment(rainyNight, RAINY, 0.95f);
	addEnvironment(rainyNight, RAINY, 1000.0f);
		
}

void EnvironmentController :: setEnvironment(float t)
{
	Environment a, b, c;
	CLAMP(humidity,0.0f,1.0f)

	if (humidity >= 1.0f)
	{
		c = extractEnvironment(t, rainy);
		c.sun.a = 0;
		skyBox->afterglowColor.set(Color(0,0,0,0));
	}
	else if (humidity > 0.5)
	{
		a = extractEnvironment(t, foggy);
		b = extractEnvironment(t, rainy);
		c = Environment::interpolate(a,b,2.0f*(humidity-0.5f));
		c.sun.a = c.sun.a*(0.8f - humidity)/0.3f;
		if (c.sun.a < 0) c.sun.a = 0;
		skyBox->afterglowColor.set(c.sun.a*4.0f*(1.0f - humidity)*(1.0f - humidity)*Color(1.0f,1.0f,0,1.0f));
	}
	else if (humidity == 0.5)
	{
		c = extractEnvironment(t, foggy);
		skyBox->afterglowColor.set(Color(1.0f,1.0f,0,1.0f));
	}
	else if (humidity > 0.0f)
	{
		a = extractEnvironment(t, sunny);
		b = extractEnvironment(t, foggy);
		c = Environment::interpolate(a,b,2.0f*humidity);
		skyBox->afterglowColor.set(
				(1.0f - 2.0f*(humidity))	* Color(1.0f,0.2f,0,1.0f) 
			+ 2.0f*(humidity)					* Color(1.0f,1.0f,0,1.0f));
	}
	else
	{
		c = extractEnvironment(t, sunny);
		skyBox->afterglowColor.set(Color(1.0f,0.2f,0,1.0f));
	}
	skyBox->assumeEnvironment(c);	
}

Environment EnvironmentController :: extractEnvironment(float t, const EnvironmentMap& map)
{
	float ta, tb, dt;
	Environment a, b, out;
	EnvironmentMap::const_iterator it;
	it = map.upper_bound(t);
	b = it->second;
	tb = it->first;
	it--;
	a = it->second;
	ta = it->first;
	if (tb != ta)
	{
		dt = (t - ta)/(tb - ta);
		out = Environment::interpolate(a,b,dt);
	}
	else return a;
	return out;
}

void EnvironmentController :: addEnvironment( 
	Environment t, WeatherType weather, float time)
{
	switch (weather)
	{
		case RAINY:	
			{
				rainy[time] = t;
			}
			break;
		case FOGGY:	
			{
				foggy[time] = t;
			}
			break;
		case SUNNY:	
			{
				sunny[time] = t;
			}
			break;
	}
}
