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


#include "FieldDisplay.h"
#include "SkyBox.h"
#include "Scene.h"
#include "Debug.h"
#include "GraphicsBridge.h"
#include "GL/glu.h"
#include "boinc_api.h"

#define RENDER_ARROW(MID, TIP, WID)	\
	glTexCoord2f(0,0.5f); \
	chartRenderer.ambientColor.setTransparent(); \
	glVertex3f( mid.x - WID.x, mid.y, mid.z - WID.z );\
	glTexCoord2f(0,0.5f);\
	chartRenderer.ambientColor.setOpaque(); \
	glVertex3f( MID.x - WID.x, MID.y, MID.z - WID.z );\
	glTexCoord2f(1,0.5f);\
	glVertex3f( MID.x + WID.x, MID.y, MID.z + WID.z );\
	glTexCoord2f(1,0.5f);\
	chartRenderer.ambientColor.setTransparent(); \
	glVertex3f( mid.x + WID.x, mid.y, mid.z + WID.z );\
	glTexCoord2f(0,0.5f);\
	chartRenderer.ambientColor.setOpaque(); \
	glVertex3f( MID.x - WID.x, MID.y, MID.z - WID.z );\
	glTexCoord2f(0,1);\
	glVertex3f( TIP.x - WID.x, TIP.y, TIP.z - WID.z );\
	glTexCoord2f(1,1);\
	glVertex3f( TIP.x + WID.x, TIP.y, TIP.z + WID.z );\
	glTexCoord2f(1,0.5f);\
	glVertex3f( MID.x + WID.x, MID.y, MID.z + WID.z );


//	c'tor
FieldDisplay :: FieldDisplay(DisplayMM* display, int maxSampleCount, int sampleSize, float3 dim)
:	display(display),
	maxSampleCount(maxSampleCount),
	skyBox(display->skyBox),
	scene(display->scene),
	width(dim.x),	depth(dim.y), height(dim.z),
	newSamples(0), sampleSize(sampleSize),
	dataReady(true), dataRead(true), soft(true),
	timeSinceData(0.0f), dataOffset(0.0f),
	dataThroughput(1.45f), freshness(5),
	arrowWidth(0.02f), arrowThickness(0.004f),
	chartRenderer(this)
{	const int MAX_LENGTH = 1000;
	char imagefile [MAX_LENGTH];
    int retval;

	retval = boinc_resolve_filename("specular.png",imagefile,MAX_LENGTH);
	chartRenderer.specularTex = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE, GL_CLAMP_TO_EDGE);
	retval = boinc_resolve_filename("diagram.png",imagefile,MAX_LENGTH);
	chartRenderer.graphTex = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE, GL_CLAMP_TO_EDGE);
	retval = boinc_resolve_filename("arrowhead.png",imagefile,MAX_LENGTH);
	arrowTex = textureLoader.loadTexture2D(imagefile, RGBA_TEXTURE, GL_CLAMP_TO_EDGE);
	chartRenderer.environmentTex = skyBox->gray;
	
	time = SurfaceProvider :: getInstance() -> getLine();
	time->print("time");
	age = SurfaceProvider :: getInstance() -> getLine();	
	age->print("age");
	infectiousness = SurfaceProvider :: getInstance() -> getLine();	
	infectiousness->print("infectiousness");

	backBuffer = &data0;
	frontBuffer = &data1;
	float* f;
	LOG("sampleSize: " << sampleSize)
	for (int i = 0; i < maxSampleCount + 1; i++)
	{
		f = new float[sampleSize];
		for (int j = 0; j < sampleSize; j++) f[j] = 0.1f;
		data.insert(data.end(), f);
	}

	chartRenderer.sampleSize = sampleSize;
	chartRenderer.maxSampleCount = maxSampleCount;
}
    
//	render
void	FieldDisplay :: render()
{ 
	double mvm[16], prm[16];
	int vPort[4] = { 0, 0, 100, 100};

	glGetDoublev(GL_MODELVIEW_MATRIX, mvm);
	//glGetDoublev(GL_PROJECTION_MATRIX, prm);	
	for (int i = 0; i < 16; i++)
		prm[i] = i%4 == i/4 ? 1.0f : 0.0f;

	double ox, oy, oz;
	gluUnProject(50.0, 50.0, 1.0, mvm, prm, vPort, &ox, &oy, &oz);
	//gluProject(0,0,0, mvm, prm, vPort, &ox, &oy, &oz);

	float dataBegin = dataOffset;
	CLAMP(dataBegin, -1.0f, 1.0f);
	SampleList::iterator s = data.begin(), end = data.end();
	end--; end--;

	double phiDeg = Debug :: doubles[0];
	double phi = PI*phiDeg/180.0;
	double sinPhi = sin(phi), cosPhi = cos(phi);

	chartRenderer.origin = float3(ox*cosPhi - oz*sinPhi, oy, oz*cosPhi + ox*sinPhi);	
	chartRenderer.dd = depth / (float)(soft? sampleSize - 1 : sampleSize),
	chartRenderer.dw = width / (float)(maxSampleCount - 2),
	chartRenderer.w0 = -width/2.0f - dataBegin*chartRenderer.dw, 
	chartRenderer.d0 = -depth/2.0f, 
	chartRenderer.h0 = -height/2.0f,
	chartRenderer.h = height,
	chartRenderer.o = dataBegin;
	chartRenderer.depth = depth;
	chartRenderer.width = width;	
	chartRenderer.sun = float3(
		skyBox->sunPosition.x*cosPhi - skyBox->sunPosition.z*sinPhi, 
		skyBox->sunPosition.y, 
		skyBox->sunPosition.z*cosPhi + skyBox->sunPosition.x*sinPhi);
	chartRenderer.lightX = float3(
		skyBox->sunX.x*cosPhi - skyBox->sunX.z*sinPhi, 
		skyBox->sunX.y, 
		skyBox->sunX.z*cosPhi + skyBox->sunX.x*sinPhi);
	chartRenderer.lightY = float3(
		skyBox->sunY.x*cosPhi - skyBox->sunY.z*sinPhi, 
		skyBox->sunY.y, 
		skyBox->sunY.z*cosPhi + skyBox->sunY.x*sinPhi);
	chartRenderer.ambientColor = 1.9f*skyBox->ambientColor;
	chartRenderer.specularColor = skyBox->sunColor*scene->occlusion,
	chartRenderer.diffuseColor = Color(0.01f, 0.04f, 0.08f);//(skyBox->ambientColor + skyBox->sunlightColor)/2.0f;
	if (skyBox->box) 
	{
		chartRenderer.environmentTex = skyBox->box;
		chartRenderer.z_scale = -1.0f;
		chartRenderer.specularPass = false;
	}
	else 
	{
		chartRenderer.environmentTex = skyBox->gray;
		chartRenderer.z_scale = +1.0f;
		chartRenderer.specularPass = true;
	}

	glPushMatrix();
	glRotatef(phiDeg,0,1,0);
	chartRenderer.prepVertices(s, end, soft);	
	unsigned int wdth = (unsigned int)chartRenderer.vertices.size()/sampleSize;
  	
		glEnable(GL_CULL_FACE);
		glEnable(GL_DEPTH_TEST);	
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE);
		glDisable(GL_TEXTURE_2D);
		glDisable(GL_BLEND);
		glDepthFunc(GL_GREATER);
		glDepthRange(0.0, 1.0);
		chartRenderer.pass = ChartRenderer::BACKSIDE_PASS;
		glCullFace(GL_FRONT);
    
	if (soft) 
	{
		chartRenderer.prepNormals(wdth);
		chartRenderer.interpolate(0, wdth);
		chartRenderer.renderSoft(wdth);		
			
			renderArrows(phiDeg);			
    		
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthFunc(GL_LEQUAL);
		glCullFace(GL_BACK);
		chartRenderer.pass = ChartRenderer::DEPTH_PASS;
		chartRenderer.renderSoft(wdth);
		
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE);
		glDepthFunc(GL_EQUAL);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		
		chartRenderer.pass = ChartRenderer::FRONTSIDE_PASS;
		chartRenderer.renderSoft(wdth);
		
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE); 
		glEnable(GL_TEXTURE_2D);		
		chartRenderer.pass = ChartRenderer::COLOR_PASS;				
		chartRenderer.renderSoft(wdth);	
	}
	else
	{			
		chartRenderer.renderHard(wdth);

		renderArrows(phiDeg);
    		
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
		glDepthFunc(GL_LEQUAL);
		glCullFace(GL_BACK);
		glDepthRange(0.0, 1.0);
		chartRenderer.pass = ChartRenderer::DEPTH_PASS;
		chartRenderer.renderHard(wdth);
		
		glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_TRUE);
		glDepthFunc(GL_EQUAL);
		glEnable(GL_BLEND);
		glBlendFunc(GL_ONE_MINUS_DST_ALPHA, GL_ONE_MINUS_SRC_ALPHA);		
		chartRenderer.pass = ChartRenderer::FRONTSIDE_PASS;
		chartRenderer.renderHard(wdth);
		
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE); 
		glEnable(GL_TEXTURE_2D);		
		chartRenderer.pass = ChartRenderer::COLOR_PASS;				
		chartRenderer.renderHard(wdth);
	}

		glDepthFunc(GL_LEQUAL);
		glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glPopMatrix();
	glColor3f(1,0.9f,0.0f);
}

//	update
void FieldDisplay :: update(float deltaT)
{
	timeSinceData += deltaT;	
	if (dataOffset < 1.0f) dataOffset += deltaT*dataThroughput;

	SampleList* buffer;
	
	if (frontBuffer->empty()) return;
	
//	GraphicsBridge :: lock();
	buffer = frontBuffer;
	frontBuffer = backBuffer;
	backBuffer = buffer;
//	GraphicsBridge :: unlock();

	SampleList::iterator	b = backBuffer->begin(), e = backBuffer->end();
	newSamples = (int) backBuffer->size();

	if (timeSinceData > 0.0f)
		dataThroughput = (float)newSamples/timeSinceData;
	
	timeSinceData = 0.0f;
	Sample terminated;

	while (b != e)
	{
		data.insert(data.end(), *b);
		terminated = *data.begin();
		data.pop_front();
		delete[] terminated;
		dataOffset -= 1.0f;
		if (dataOffset < 0.0f) dataOffset = 0.0f;
		b++;
	}
	backBuffer->clear();	
}

//	addData
void FieldDisplay :: addData(Sample s)	// USED BY OTHER THREAD
{ 	
//	GraphicsBridge :: lock();
	frontBuffer->insert(frontBuffer->end(), s);
//	GraphicsBridge :: unlock();
}

//	renderArrows
void	FieldDisplay :: renderArrows(double phi)
{
	skyBox->ambientColor.setOpaque();	
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_ALPHA_TEST);
	glAlphaFunc(GL_GREATER, 0.7f);
	glDisable(GL_CULL_FACE);
	glClearDepth(1.0f);
	glClear(GL_DEPTH_BUFFER_BIT);
	glDepthRange(0.0f, 1.0f);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_FALSE);
	glEnable(GL_BLEND);
	glDepthFunc(GL_LESS);

	float dh = 0.03f;
	float h0 = chartRenderer.h0 - dh;

	float3 mid = float3(-width/2.0f, h0, -depth/2.0f );

	float3 tTip, aTip, iTip;
	tTip = float3( width, h0, -depth/2.0f);
	aTip = float3(-width/2.0f, h0, depth);
	iTip = float3(-width/2.0f, h0 + 1.3f*height, -depth/2.0f);

	float k0 = 1.0f/20.0f, k1 = 1.0f/16.0f, k2 = 1.0f / 0.3f;
	float k1a = 1.0f - k1;

	float3	tMid = k1a*tTip + k1*mid, 
					aMid = k1a*aTip + k1*mid, 
					iMid = k1a*iTip + k1*mid;

	float3	tWid = k0*(aTip - mid),
					aWid = k0*(tTip - mid),
					iWid = k0*(iTip - mid).cross(chartRenderer.origin - mid).direction()*k2;	

	glBindTexture(GL_TEXTURE_2D, arrowTex); 
	glBegin(GL_QUADS);
	 
		//	time
		RENDER_ARROW(tMid, tTip, tWid)

		//	age
		RENDER_ARROW(aMid, aTip, aWid)

    //	infectiousness
		glTexCoord2f(0,0.5f);
		chartRenderer.ambientColor.setTransparent();
		glVertex3f( mid.x - iWid.x, h0, mid.z - iWid.z );
		glTexCoord2f(0,0.5f);
		chartRenderer.ambientColor.setOpaque();
		glVertex3f( iMid.x - iWid.x, iMid.y, iMid.z - iWid.z );

		glTexCoord2f(1,0.5f);
		glVertex3f( iMid.x + iWid.x, iMid.y, iMid.z + iWid.z );
		glTexCoord2f(1,0.5f);
		chartRenderer.ambientColor.setTransparent();
		glVertex3f( mid.x + iWid.x, h0, mid.z + iWid.z );  		
		
		glTexCoord2f(0,0.5f);
		chartRenderer.ambientColor.setOpaque();
		glVertex3f( iMid.x - iWid.x, iMid.y, iMid.z - iWid.z );
		glTexCoord2f(0,1);
		glVertex3f( iTip.x - iWid.x, iTip.y, iTip.z - iWid.z );

		glTexCoord2f(1,1);
		glVertex3f( iTip.x + iWid.x, iTip.y, iTip.z + iWid.z );
		glTexCoord2f(1,0.5f);
		glVertex3f( iMid.x + iWid.x, iMid.y, iMid.z + iWid.z );		

	glEnd();

	float2 letterSize = float2(0.34f, 0.46f);

	glPushMatrix();		
		glTranslatef(F3_ARG(tTip));
		scene->viewController->unrotate(phi);
		time->render(letterSize, float2(0.5f,1.5f));
	glPopMatrix();

	glPushMatrix();		
		glTranslatef(F3_ARG(aTip));
		scene->viewController->unrotate(phi);
		age->render(letterSize, float2(0.5f,1.5f));
	glPopMatrix();

	glPushMatrix();		
		glTranslatef(F3_ARG(iTip));
		scene->viewController->unrotate(phi);
		infectiousness->render(letterSize, float2(0.5f,1.5f));
	glPopMatrix();

	glDisable(GL_BLEND);
	glDisable(GL_TEXTURE_2D);
	glDisable(GL_ALPHA_TEST);
	glEnable(GL_CULL_FACE);
}

