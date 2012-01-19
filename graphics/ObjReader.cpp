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


#include "ObjReader.h"
#include "macros.h" 
#include <iostream>

ObjReader :: ObjReader(String filename)
{
	StringBuffer sb = StringBuffer();
	std::ifstream file(	filename.c_str() );
	file.seekg (0, std::ios::end);
	unsigned int length = file.tellg();
	file.seekg (0, std::ios::beg);	
	ITERATE(i, length) sb.append((char)file.get());
	file.close();	
	string = sb.getString();
}

ObjReader :: ~ObjReader()
{
}

Mesh*		
ObjReader :: readMesh(String textureDirectory, float scale)
{
	this->scale = scale;
	index = (char*) string.c_str();
	char* end = index + string.length();
	while(index < end)
	{
		readLine();
		handleLine();
	}
	Mesh* m = new Mesh(textureDirectory);
	SegmentMap::iterator i = segments.begin();
	SegmentMap::iterator t = segments.end();
	while (i != t)
	{
		m->addSegment(i->second);
		i++;
	}
	return m;
}

void
ObjReader :: handleLine()
{ 	
	switch (line[0])
	{
		case '#': 
			break;
		case 'v': 
		{
			if (line[1] == ' ') 
			{
				lineIndex = 2;
				vertices.push_back(parseVector3());
			}
			else if (line[1] == 'n') 
			{
				lineIndex = 3;
				normals.push_back(parseVector3());
			}
			else if (line[1] == 't') 
			{
				lineIndex = 3;
				texCoords.push_back(parseVector2());
			}
			break;
		}			
		case 'o': 
		{
			String name = readString(2);
			currentSegment = new Segment();
			segments[name] = currentSegment;
			currentSegment->name = name;
			currentSegment->texture = "*";
			break;
		}
		case 'f': 
		{
			int triplets = 0;
			ITERATE(i, line.length())
			{
				if (line[i] == '/') triplets++;
			}
			triplets /= 2;

			Triplet t0, t1, t2, t3;
			lineIndex = 2;

			if (triplets == 3)
			{
				t0 = parseTriplet();
				t1 = parseTriplet();
				t2 = parseTriplet();
				currentSegment->triangles.push_back(new Triangle(t0,t1,t2));
			}
			else if (triplets == 4)
			{
				t0 = parseTriplet();
				t1 = parseTriplet();
				t2 = parseTriplet();
				t3 = parseTriplet();
				currentSegment->triangles.push_back(new Triangle(t0,t1,t2));
				currentSegment->triangles.push_back(new Triangle(t0,t2,t3));
			}			
			break;
		}	
		case 'u':
		{
			if (	 line[1] == 's' && line[2] == 'e' && line[3] == 'm' 
					&& line[4] == 'a' && line[5] == 'p' && line[6] == ' ')
			{
				currentSegment->texture = readString(7);
			}
			break;
		}
	}
}   

void  
ObjReader :: readLine()
{
	StringBuffer sb = StringBuffer();
	char c = *index; index++;
	sb.append(c);
	while (c != '\n')
	{
		c = *index; index++;
		sb.append(c);
	}
	sb.undo();
	line = sb.getString();
	lineIndex = 0;
}

String  
ObjReader :: readString(int index)
{
	StringBuffer sb = StringBuffer();
	for (unsigned int i = index; i < line.length(); i++)
		sb.append(line[i]);
	return sb.getString();
}

float3	
ObjReader :: parseVector3()
{
	double x, y, z;
	x = parseDouble();
	z = -parseDouble();
	y = parseDouble(); 	
	return float3(x,y,z);
}

float2	
ObjReader :: parseVector2()
{
	double x, y;
	x = parseDouble();
	y = parseDouble();	
	return float2(x,y);
}

Triplet	
ObjReader :: parseTriplet()
{
	int v,t,n;
	v = parseInt();
	t = parseInt();
	n = parseInt();	
	Triplet tpt;
	tpt.vertex = scale*vertices[v - 1];
	tpt.texture = texCoords[t - 1];
	tpt.texture.y = 1.0f - tpt.texture.y;
	tpt.normal = normals[n - 1];
	return tpt;
}

double	
ObjReader :: parseDouble()
{
	int afterComma = 0;
	double out = 0.0, sign = 1.0, exp = 10.0; 
	char c = line[lineIndex];

	if (c == '-') 
	{
		sign = -1.0;
		lineIndex++;
		c = line[lineIndex];
	}	

	while(c != ' ' && c != 0)
	{ 		
		if (c == '.')
		{
			afterComma = true;
		}
		else if (afterComma)
		{
			out += (double)(line[lineIndex] - '0')/exp;
			exp *= 10.0;
		}
		else
		{
			out *= 10.0;
			out += (double)(line[lineIndex] - '0');
		}
		lineIndex++;
		c = line[lineIndex];
	}
	lineIndex++;
	return out*sign;
}

int			
ObjReader :: parseInt()
{
	int out = 0;	
	char c = line[lineIndex];

	while(c >= '0' && c <= '9')
	{ 		
		out *= 10;
		out += (int)(line[lineIndex] - '0');
		lineIndex++;
		c = line[lineIndex];
	}
	lineIndex++;
	return out;
}

void		
ObjReader :: test()
{
}
