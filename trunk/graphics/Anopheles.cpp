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


#include "Anopheles.h"
#include "ObjReader.h"
#include "boinc_api.h"

Mesh*			Anopheles :: mesh = 0;
Segment *	Anopheles :: abdomen = 0;
Segment	*	Anopheles :: head = 0;
Segment	*	Anopheles :: torso = 0; 
Segment * Anopheles :: leftWing = 0; 
Segment	*	Anopheles :: rightWing = 0;

void Anopheles :: init()
{
	const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	retval = boinc_resolve_filename("anopheles_004.obj",imagefile,MAX_LENGTH);
	ObjReader orMM = ObjReader(imagefile);
	mesh = orMM.readMesh(".", 0.023f);
	
	SegmentList& segments = mesh->getSegments();
	SegmentList wings, body;
	ITERATE_LIST(Segment*, segments, i, t,
		switch ((*i)->name[0]){	
			//abdomen
			case 'A':{
				abdomen = *i;
				body.push_back(*i); break;
			}
			//head
			case 'H':{
				head = *i;
				body.push_back(*i); break;
			}
			//leftWing
			case 'L':{
				leftWing = *i;
				wings.push_back(*i); break;
			}
			//rightWing
			case 'R':{
				rightWing = *i;
				wings.push_back(*i); break;
			}
			//torso
			case 'T':{
				torso = *i;
				body.push_back(*i); break;
	}	}	)
	segments.clear();
	ITERATE_LIST(Segment*, body, ib, tb, segments.push_back(*ib);)
	ITERATE_LIST(Segment*, wings, iw, tw, segments.push_back(*iw);)
}	
