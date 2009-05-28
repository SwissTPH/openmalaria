# Copyright © 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
# Licence: GNU General Public Licence version 2 or later (see COPYING)

noinst_PROGRAMS = malariaGraphics

malariaGraphics_SOURCES = Anopheles.cpp AlphaConfiguration.cpp ChartRenderer.cpp DataConfiguration.cpp DemoConfiguration.cpp \
				Display.cpp FieldDisplay.cpp LineChart.cpp Debug.cpp Color.cpp TextureLoader.cpp \
				MathFunctions.cpp double3.cpp float3.cpp int2.cpp int3.cpp plane.cpp Mesh.cpp \
				NormalMapBuilder.cpp ObjReader.cpp Segment.cpp Overlay.cpp ProgressBar.cpp AnophelesRenderer.cpp \
				Environment.cpp SceneController.cpp Scene.cpp ViewController.cpp CubeRenderer.cpp \
				EnvironmentController.cpp MultiTextureHelper.cpp PreRenderedBox.cpp SkyBox.cpp BoincBridge.cpp \
				CallBackBridge.cpp GL_Window.cpp GraphicsBridge.cpp KeyHandler.cpp main.cpp MouseHandler.cpp \
				SystemTimer.cpp Font.cpp Label.cpp Line.cpp Surface.cpp SurfaceProvider.cpp

malariaGraphics_LDADD =  -lboinc_api -lboinc -lboinc_graphics2 -lglut -lGL -lGLU -lm -lxerces-c -lIL -lILU -lILUT  

malariaGraphics_LDFLAGS = @BOINC_LDFLAGS@

AM_CPPFLAGS = @CPPFLAGS@ @BOINC_CFLAGS@ -I$(top_srcdir)/include -I$(top_srcdir)/include/glew -I$(top_srcdir)/cpp/include
