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


#include "ChartRenderer.h"
#include "FieldDisplay.h"
#include "SkyBox.h"
#include "Scene.h"
#include "GL/glu.h"

//	c'tor
ChartRenderer :: ChartRenderer(FieldDisplay* chart)
:	chart(chart)
{}

//	prepVertices
void ChartRenderer :: prepVertices
	(SampleList::iterator start, SampleList::iterator end, bool soft)
{
	vertices.clear();
	normals.clear();

	Sample s0, s1;
	s0 = *start; start++; s1 = *start;
	float a, b, d; 	

	dMin = 1000000.0f;
	dMax = 0.0f;

	ITERATE(i, sampleSize)
	{
		a = s0[i]; b = s1[i];	
		a3 = float3(w0,			 h0 + h*a, d0 + i*dd);
		b3 = float3(w0 + dw, h0 + h*b, d0 + i*dd);	
		a3 = o*b3 + (1.0f - o)*a3; 	 
		vertices.push_back(a3);
		d = (a3 - origin).length();
		if (d < dMin) dMin = d;
		else if (d > dMax) dMax = d;
	}
	w0 += dw; 	

	while (start != end)
	{
		s0 = *start;
		ITERATE(i, sampleSize)
		{
			a = s0[i];
			a3 = float3(w0,	h0 + h*a, d0 + i*dd);
			vertices.push_back(a3);		
			d = (a3 - origin).length();
			if (d < dMin) dMin = d;
			else if (d > dMax) dMax = d;
		}
		start++;
		w0 += dw;
	}	 
	
	s1 = *end;
  ITERATE(i, sampleSize)
	{
		a = s0[i]; b = s1[i];  
		a3 = float3(w0 - dw, h0 + h*a, d0 + i*dd);
		b3 = float3(w0, h0 + h*b, d0 + i*dd); 
		a3 = o*b3 + (1.0f - o)*a3;
		vertices.push_back(a3);		
		d = (a3 - origin).length();
		if (d < dMin) dMin = d;
		else if (d > dMax) dMax = d;
	}

	dSpan = dMax - dMin; 
}

//	prepNormals
void ChartRenderer :: prepNormals(unsigned int wdth)
{
	PUSH_NORMAL(sampleSize,0,1,0);	
	for (unsigned int i = 1; i < sampleSize - 1; i++)
		{PUSH_NORMAL(sampleSize+i, i, i+1, i-1); }
	PUSH_NORMAL(2*sampleSize - 1, sampleSize - 1, sampleSize-1, sampleSize-2);

	for (unsigned int x = 1; x < wdth - 1; x++)
	{
		PUSH_NORMAL((x+1)*sampleSize, (x-1)*sampleSize, x*sampleSize+1, x*sampleSize);
		for (unsigned int i = 1; i < sampleSize - 1; i++)
			{PUSH_NORMAL((x+1)*sampleSize+i, (x-1)*sampleSize+i, x*sampleSize+i+1, x*sampleSize+i-1);}
		PUSH_NORMAL((x+2)*sampleSize-1, x*sampleSize-1, (x+1)*sampleSize-1, (x+1)*sampleSize-2);
	}
	
	PUSH_NORMAL(sampleSize*(wdth - 1), sampleSize*(wdth - 2), sampleSize*(wdth - 1)+1, sampleSize*(wdth - 1));
	for (unsigned int i = 1; i < sampleSize - 1; i++)
		{PUSH_NORMAL(sampleSize*(wdth-1)+i, sampleSize*(wdth-2)+i, sampleSize*(wdth-1)+i+1, sampleSize*(wdth-1)+i-1);}
	PUSH_NORMAL(sampleSize*wdth-1, sampleSize*(wdth-1)-1, sampleSize*wdth-1, sampleSize*wdth-2);
}

//	renderHard
void ChartRenderer :: renderHard(unsigned int wdth)
{
	float3 m0, m1, m2;
	float t, a, b;

	ITERATE(x, wdth - 1)
	{
		a3 = vertices[x*sampleSize];
		b3 = vertices[(x+1)*sampleSize];

		m0 = float3(a3.x, h0, a3.z);
		m1 = float3(b3.x, h0, b3.z);				
		
		renderTriangle(origin, m0, b3, a3);
		renderTriangle(origin, m0, m1, b3);

		a3 = vertices[x*sampleSize + sampleSize - 1];
		b3 = vertices[(x+1)*sampleSize + sampleSize - 1];
		c3 = float3(b3.x, b3.y, b3.z + dd);
		d3 = float3(a3.x, a3.y, a3.z + dd);

		m0 = float3(d3.x, h0, d3.z);
		m1 = float3(c3.x, h0, c3.z);

		renderTriangle(origin, a3, b3, c3);
		renderTriangle(origin, a3, c3, d3);

		an = float3(0,0,1);

		renderTriangle(origin, m0, d3, c3, an);
		renderTriangle(origin, m1, m0, c3, an);

		ITERATE(y, sampleSize - 1)
		{ 			
			a3 = vertices[x*sampleSize + y];
			b3 = vertices[(x+1)*sampleSize + y];
			c3 = vertices[(x+1)*sampleSize + y + 1];
			d3 = vertices[x*sampleSize + y + 1];

			m0 = float3(a3.x, a3.y, d3.z);
			m1 = float3(b3.x, b3.y, c3.z);

			renderTriangle(origin, a3, b3, m0);
			renderTriangle(origin, b3, m1, m0);

			a = d3.y - a3.y;
			b = c3.y - b3.y;

			if (a*b >= 0.0f)
			{				
				renderTriangle(origin, m0, m1, c3);
				renderTriangle(origin, m0, c3, d3);
			}
			else
			{				
				t = a / (b - a);
				m2 = m0 - t*(m1 - m0);			
				renderTriangle(origin, m2, m1, c3);
				renderTriangle(origin, m0, m2, d3);
			}
		}
	}

	ITERATE(y, sampleSize - 1)
	{
		a3 = vertices[y];
		b3 = float3(a3.x, h0, a3.z);		
		d3 = vertices[y + 1];
		c3 = float3(d3.x, h0, d3.z);
		d3.y = a3.y;
		an = float3(-1.0f, 0, 0);

		renderTriangle(origin, b3, a3, c3, an);
		renderTriangle(origin, d3, c3, a3, an);

		a3 = vertices[y + sampleSize*(wdth - 1)];
		b3 = float3(a3.x, h0, a3.z);		
		d3 = vertices[y + sampleSize*(wdth - 1) + 1];
		c3 = float3(d3.x, h0, d3.z);
		d3.y = a3.y;
		an = float3(1.0f, 0, 0);

		renderTriangle(origin, a3, b3, c3, an);
		renderTriangle(origin, c3, d3, a3, an);
	}

	a3 = vertices[sampleSize - 1];
	b3 = a3; d3 = a3;
	b3.y = h0; d3.z += dd;
	c3 = float3(d3.x, b3.y, d3.z);

	an = float3(-1.0f, 0, 0);

	renderTriangle(origin, a3, c3, b3, an);
	renderTriangle(origin, d3, c3, a3, an);		

	a3 = vertices[wdth*sampleSize - 1];
	b3 = a3; d3 = a3;
	b3.y = h0; d3.z += dd;
	c3 = float3(d3.x, b3.y, d3.z);

	an.x = 1.0f;

	renderTriangle(origin, c3, a3, b3, an);
	renderTriangle(origin, c3, d3, a3, an);

	an.x = 0.0f; an.y = -1.0f;
	a3 = float3(-width/2.0f, h0, -depth/2.0f);
	c3 = float3( width/2.0f, h0,  depth/2.0f);
	b3 = float3(a3.x, h0, c3.z);
	d3 = float3(c3.x, h0, a3.z);

	renderTriangle(origin, a3, b3, c3, an);
	renderTriangle(origin, a3, c3, d3, an);
}
	
//	interpolate
void ChartRenderer :: interpolate(int lod, unsigned int wdth)
{
	triangles.clear();

	float3	xa, xb, xc, xd,
					nxa, nxb, nxc, nxd,
					xab, xbc, xcd, xda,
					nxab, nxbc, nxcd, nxda;

	ITERATE(y, sampleSize - 1)
	{
		a3 = vertices[y];
		b3 = float3(a3.x, h0, a3.z);		
		d3 = vertices[y + 1];
		c3 = float3(d3.x, h0, d3.z);
		s3 = (a3 + d3)/2.0f;
		an = float3(-1.0f, 0, 0);

		triangles.push_back(Triangle(b3, a3, s3, an, an, an));
		triangles.push_back(Triangle(c3, b3, s3, an, an, an));
		triangles.push_back(Triangle(d3, c3, s3, an, an, an));

		a3 = vertices[y + sampleSize*(wdth - 1)];
		b3 = float3(a3.x, h0, a3.z);		
		d3 = vertices[y + sampleSize*(wdth - 1) + 1];
		c3 = float3(d3.x, h0, d3.z);
		s3 = (a3 + d3)/2.0f;
		an = float3(1.0f, 0, 0);

		triangles.push_back(Triangle(a3, b3, s3, an, an, an));
		triangles.push_back(Triangle(b3, c3, s3, an, an, an));
		triangles.push_back(Triangle(c3, d3, s3, an, an, an));
	}

	ITERATE(x, wdth - 1)
	{
		a3 = vertices[x*sampleSize];
		b3 = float3(a3.x, h0, a3.z);		
		d3 = vertices[(x+1)*sampleSize];
		c3 = float3(d3.x, h0, d3.z);
		s3 = (a3 + d3)/2.0f;
		an = float3(0, 0, -1.0f);

		triangles.push_back(Triangle(a3, b3, s3, an, an, an));
		triangles.push_back(Triangle(b3, c3, s3, an, an, an));
		triangles.push_back(Triangle(c3, d3, s3, an, an, an));  		

		a3 = vertices[(x+1)*sampleSize - 1];
		b3 = float3(a3.x, h0, a3.z);		
		d3 = vertices[(x+2)*sampleSize - 1];
		c3 = float3(d3.x, h0, d3.z);
		s3 = (a3 + d3)/2.0f;
		an = float3(0, 0, 1.0f);

		triangles.push_back(Triangle(b3, a3, s3, an, an, an));
		triangles.push_back(Triangle(c3, b3, s3, an, an, an));
		triangles.push_back(Triangle(d3, c3, s3, an, an, an));

		ITERATE(y, sampleSize - 1)
		{ 			
			a3 = vertices[x*sampleSize + y];
			b3 = vertices[(x+1)*sampleSize + y];
			c3 = vertices[(x+1)*sampleSize + y + 1];
			d3 = vertices[x*sampleSize + y + 1];

			an = normals[x*sampleSize + y];
			bn = normals[(x+1)*sampleSize + y];
			cn = normals[(x+1)*sampleSize + y + 1];
			dn = normals[x*sampleSize + y + 1];

			s3 = (a3 + b3 + c3 + d3)/4.0f;
			sn = ((an + bn + cn + dn)/4.0f).direction();

			xa = (d3 + a3 + b3) / 3.0f;
			xb = (a3 + b3 + c3) / 3.0f;
			xc = (b3 + c3 + d3) / 3.0f;
			xd = (c3 + d3 + a3) / 3.0f;

			nxa = ((dn + an + bn) / 3.0f).direction();
			nxb = ((an + bn + cn) / 3.0f).direction();
			nxc = ((bn + cn + dn) / 3.0f).direction();
			nxd = ((cn + dn + an) / 3.0f).direction();

			xab = (a3 + b3) / 2.0f;
			xbc = (b3 + c3) / 2.0f;
			xcd = (c3 + d3) / 2.0f;
			xda = (d3 + a3) / 2.0f;

			nxab = ((an + bn) / 2.0f).direction();
			nxbc = ((bn + cn) / 2.0f).direction();
			nxcd = ((cn + dn) / 2.0f).direction();
			nxda = ((dn + an) / 2.0f).direction();

			triangles.push_back(Triangle(xa, xb, s3, nxa, nxb, sn));
			triangles.push_back(Triangle(xb, xc, s3, nxb, nxc, sn));
			triangles.push_back(Triangle(xc, xd, s3, nxc, nxd, sn));
			triangles.push_back(Triangle(xd, xa, s3, nxd, nxa, sn));

			triangles.push_back(Triangle(xb, xa, xab, nxb, nxa, nxab));
			triangles.push_back(Triangle(xc, xb, xbc, nxc, nxb, nxbc));
			triangles.push_back(Triangle(xd, xc, xcd, nxd, nxc, nxcd));
			triangles.push_back(Triangle(xa, xd, xda, nxa, nxd, nxda));

			triangles.push_back(Triangle(xa, a3, xab, nxa, an, nxab));
			triangles.push_back(Triangle(xb, b3, xbc, nxb, bn, nxbc));
			triangles.push_back(Triangle(xc, c3, xcd, nxc, cn, nxcd));
			triangles.push_back(Triangle(xd, d3, xda, nxd, dn, nxda));	

			triangles.push_back(Triangle(b3, xb, xab, bn, nxb, nxab));
			triangles.push_back(Triangle(c3, xc, xbc, cn, nxc, nxbc));
			triangles.push_back(Triangle(d3, xd, xcd, dn, nxd, nxcd));
			triangles.push_back(Triangle(a3, xa, xda, an, nxa, nxda));
		}
	}

	an = float3(0.0f, -1.0f, 0.0f);
	a3 = float3(-width/2.0f, h0, -depth/2.0f);
	c3 = float3( width/2.0f, h0,  depth/2.0f);
	b3 = float3(a3.x, h0, c3.z);
	d3 = float3(c3.x, h0, a3.z);

	triangles.push_back(Triangle(a3, b3, c3, an, an, an));
	triangles.push_back(Triangle(a3, c3, d3, an, an, an));
}

//	renderSoft
void ChartRenderer :: renderSoft(unsigned int wdth)
{
	unsigned int s = (unsigned int) triangles.size();
	ITERATE(i, s)
	{
		Triangle& t = triangles[i];
		renderTriangle(origin, t.a, t.b, t.c, t.na, t.nb, t.nc);
	}
}
