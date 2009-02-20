#include <stdlib.h>
/*

 This file is part of OpenMalaria.
 
 Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

*/

#include <cstdlib>
#include "global.h"
#include "inputData.h"
#include <cmath>

/*
Contains global variables and constants and utility functions that are used in different modules.

Constants (parameter)

*/

int modelVersion;
int interval;
int intervalsPerYear;
int maxAgeIntervals;
int simulationMode;
double latentp;

vector<int> infantIntervalsAtRisk;
int infantIntervalsAtRiskSize;
vector<int> infantDeaths;
int infantDeathsSize;

void initGlobal () {

    modelVersion=get_model_version();
    interval=get_interval();
    //exit if days in year is not divisible by interval
    if ( mymodf(daysInYear, interval) !=  0) {
        exit(-12);
    }
    intervalsPerYear=daysInYear/interval;
    infantDeathsSize = intervalsPerYear;
    infantDeaths.resize(infantDeathsSize);
    infantIntervalsAtRiskSize = intervalsPerYear;
    infantIntervalsAtRisk.resize(infantIntervalsAtRiskSize);
    latentp=get_latentp();
    maxAgeIntervals=(int)get_maximum_ageyrs()*intervalsPerYear;
}

void clearGlobalParameters () {
}

int modIntervalsPerYear (int i) {

    int valmodIntervalsPerYear;
    valmodIntervalsPerYear= i % intervalsPerYear;
    if ( valmodIntervalsPerYear ==  0) {
        valmodIntervalsPerYear=intervalsPerYear;
    }
    return valmodIntervalsPerYear;
}

double mymodf(double d1, double d2) {
    //return modf(d1, &d2);
    return (int)d1 % (int)d2;
}

#ifdef _WIN32
int nearbyint(double x){
	return (int) x;
}
int round(double x){
	double intpart;
	if (modf(x,&intpart)>=.5)
		return x>=0.0?(int)ceil(x):(int)floor(x);
	else
		return x<0.0?(int)ceil(x):(int)floor(x);
}
#endif

