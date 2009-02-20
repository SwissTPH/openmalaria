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



#ifndef CONSTANT_H
#define CONSTANT_H

//Error constants
enum error{
	MISSING_VALUE = -99999
};

//Constants defining the different ITN types
enum InterventionType{
	NO_INTERVENTION = 0,
	IRS_INTERVENTION = 1,
	MDA_INTERVENTION = 2,
	VACCINE_INTERVENTION = 3,
	CHANGE_EIR_INTERVENTION = 4,
	CHANGE_HS_INTERVENTION = 5,
	IPTi_INTERVENTION = 6
};

//Constants used for the different ages of interventions
enum ageMDA{	
	MAX_AGE = 999,
	MIN_AGE = 0,
	COMPLIANCE = 1
};

#endif
