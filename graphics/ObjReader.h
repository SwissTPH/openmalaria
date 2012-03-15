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


#ifndef OBJ_READER_H
#define OBJ_READER_H

#include "StringWrapper.h"
#include "Triangle.h"
#include "StringBuffer.h"
#include "Mesh.h"
#include "math_headers.h"
#include <fstream>
#include <vector>
#include <list>
#include <map>

class Mesh;



typedef std::vector<float3> F3list;
typedef std::vector<float2> F2list;
typedef std::map<std::string, Segment*> SegmentMap;

class ObjReader
{
	public:

		ObjReader(String filename);
		~ObjReader();
    		
		String	string;
		String	line;

		Mesh*		readMesh(String textureDirectory, float scale);
		void		test();


	private:

		char*					index;
		unsigned int	lineIndex;
		float					scale;
		F3list				vertices, normals;
		F2list				texCoords;
		Segment*			currentSegment;
		SegmentMap		segments;

		void		readLine();
		String  readString(int index);
		float3	parseVector3();
		float2	parseVector2();
		Triplet	parseTriplet();
		double	parseDouble();
		int			parseInt();
		void		handleLine();
};

#endif
