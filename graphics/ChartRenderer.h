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


#ifndef CHART_RENDERER_H
#define CHART_RENDERER_H

#include "gl_headers.h"
#include "math_headers.h"
#include "ChartRendererMacros.h"
#include "Color.h"
#include "macros.h"
#include <vector>
#include <list>

typedef float* Sample; 
typedef std::list<Sample> SampleList;
typedef std::vector<float3> F3Vector;
typedef std::vector<float> FVector;
						 
class FieldDisplay;

class ChartRenderer
{   
	public:

		typedef enum {	FRONTSIDE_PASS, BACKSIDE_PASS, 
										COLOR_PASS, DEPTH_PASS, ALL_IN_ONE	} Pass;

		class Triangle 
		{
			public:

				float3	a, b, c, na, nb, nc;

				Triangle(float3	a, float3	b, float3	c, 
					float3	na, float3	nb, float3	nc)
					:	a(a), b(b), c(c), 
						na(na), nb(nb), nc(nc)
				{}
		};  

		typedef std::vector<Triangle> TriVector;

		ChartRenderer(FieldDisplay* chart);

		bool					specularPass;

		unsigned int		maxSampleCount, sampleSize,
										specularTex, graphTex, environmentTex;

		Pass						pass;

		float 					depth, width,
										dd, dw, 
										w0, d0, h0, 
										h, o,
										z_scale,
										dMin, dMax, dSpan;

		float3					origin, 
										a3, b3, c3, d3, s3, 
										an, bn, cn, dn, sn,
										sun, lightX, lightY;

		Color						specularColor, diffuseColor, ambientColor;

		F3Vector				vertices, normals;
		TriVector				triangles;
		FieldDisplay*		chart;

		void				prepVertices(SampleList::iterator s, SampleList::iterator end, bool soft);
		void				prepNormals(unsigned int wdth);
		void				interpolate(int lod, unsigned int wdth);

		void				renderHard(unsigned int wdth);
		void				renderSoft(unsigned int wdth);

		inline void	renderTriangle(	const float3 & origin, const float3 & a, 
													const float3 & b, const float3 & c)
							{
									float3	u = c - a, v = b - a;
									float3	cross = u.cross(v);
									float len = cross.length();

									if (len == 0.0f) return;

									float3	n = cross/len;
									renderTriangle( origin, a, b, c, n, n, n);
							}

		inline void	renderTriangle(	const float3 & origin, const float3 & a, 
													const float3 & b, const float3 & c, const float3 & n)
							{ 
								renderTriangle( origin, a, b, c, n, n, n);
							}
													

		inline void	renderTriangle(	const float3 & origin, 
													const float3 & a, const float3 & b, const float3 & c,
													const float3 & na, const float3 & nb, const float3 & nc)
							{
									switch (pass)
								{	
									float calib, dist, clamp;
									case FRONTSIDE_PASS:
									{	
										calib = 0.06f;	clamp = 0.15f;
										glBegin(GL_TRIANGLES);
											RENDER_FRONT(b);
											RENDER_FRONT(a);
											RENDER_FRONT(c);
										glEnd();

									}	break;

									case BACKSIDE_PASS:
									{
										calib = 0.06f;	clamp = 0.2f;
										glBegin(GL_TRIANGLES);
											RENDER_BACK(b);
											RENDER_BACK(a);
											RENDER_BACK(c);
										glEnd();

									} break;

									case DEPTH_PASS:
									{
										glBegin(GL_TRIANGLES);
											glVertex3f(F3_ARG(b));
											glVertex3f(F3_ARG(a));
											glVertex3f(F3_ARG(c));
										glEnd();

									} break;

									case COLOR_PASS:
									{
										float3	in, out;
										float2	texCoord;
										float		facing, fresnel;
										glBindTexture(GL_TEXTURE_2D, graphTex);
										glBlendFunc(GL_DST_ALPHA, GL_ONE_MINUS_DST_ALPHA);
										glColor4f(diffuseColor.r, diffuseColor.g, diffuseColor.b, 1.0f);
										glDisable(GL_TEXTURE_2D);
										
											glBegin(GL_TRIANGLES);
												RENDER_DIFFUSE(b, nb);
												RENDER_DIFFUSE(a, na);
												RENDER_DIFFUSE(c, nc);
												RENDER_DIFFUSE(b, nb);
												RENDER_DIFFUSE(a, na);	
												RENDER_DIFFUSE(c, nc);
											glEnd();

									if (specularPass)
									{
										glEnable(GL_TEXTURE_2D);
										glBindTexture(GL_TEXTURE_2D, specularTex);										
										glBlendFunc(GL_ONE, GL_ONE);//_MINUS_SRC_COLOR);

											glBegin(GL_TRIANGLES);
												RENDER_SPECULAR(b, nb);
												RENDER_SPECULAR(a, na);	
												RENDER_SPECULAR(c, nc);
											glEnd();	
									}

										glDisable(GL_TEXTURE_2D);
										glEnable(GL_TEXTURE_CUBE_MAP);
										glBindTexture(GL_TEXTURE_CUBE_MAP, environmentTex);										
										glBlendFunc(GL_SRC_ALPHA, GL_ONE);

											glBegin(GL_TRIANGLES);
												RENDER_ENVIRONMENT(b, nb);
												RENDER_ENVIRONMENT(a, na);	
												RENDER_ENVIRONMENT(c, nc);
											glEnd();
										
										glEnable(GL_TEXTURE_2D);
										glDisable(GL_TEXTURE_CUBE_MAP);

									} break;

									case ALL_IN_ONE:
									{
										float3	in, out;
										float2	texCoord;
										float		facing, fresnel;

										glEnable(GL_TEXTURE_2D);
										glBindTexture(GL_TEXTURE_2D, specularTex);
										glDisable(GL_BLEND);

											glBegin(GL_TRIANGLES);
												RENDER_SPECULAR(b, nb);
												RENDER_SPECULAR(a, na);	
												RENDER_SPECULAR(c, nc);
											glEnd();  
										
										glBindTexture(GL_TEXTURE_2D, graphTex);
										glEnable(GL_BLEND);
										glBlendFunc(GL_ONE, GL_ONE);										

											glBegin(GL_TRIANGLES);													
												RENDER_DIFFUSE(b, nb);
												RENDER_DIFFUSE(a, na);
												RENDER_DIFFUSE(c, nc);
											glEnd();

									} break;
								}
							}
};

#endif
