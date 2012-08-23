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


#ifndef ALL_MACROS_H
#define ALL_MACROS_H

#define PI 3.14159265359f
#define PI_F 3.1415926536f
#define CRASH(m) std::cout << m << '\n'; exit(-1);
#define LOG(m) std::cout << m << '\n';

#define LOG_V(X) std::cout << "X: " << X << '\n';

#define CLAMP(x,a,b) \
	if ((x) < (a)) (x) = (a); \
	else if ((x) > (b)) (x) = (b);

#define ITERATE_SET(A,S,B,C,D) \
	std::set<A>::iterator B = S.begin(); \
	std::set<A>::iterator C = S.end(); \
	while (B != C){D;B++;}

#define ITERATE_LIST(A,S,B,C,D) \
	std::list<A>::iterator B = S.begin(); \
	std::list<A>::iterator C = S.end(); \
	while (B != C){D B++;}

#define ITERATE(A,L) for (unsigned int A = 0; A < L; A++)
#define RANDOM() (float)rand()/(float)RAND_MAX
#define PLUS_RAND(X) X*(float)rand()/(float)RAND_MAX
#define SYMM_RAND(X) 2.0f*X*(float)rand()/(float)RAND_MAX - X

#define SYMM_RAND_3(X) \
	2.0f*X*(float)rand()/(float)RAND_MAX - X, \
	2.0f*X*(float)rand()/(float)RAND_MAX - X, \
	2.0f*X*(float)rand()/(float)RAND_MAX - X

#define ABS(X) (X) > 0 ? (X) : -(X)

#define COLOR_AS_ARGUMENT(C) C.r, C.g, C.b, C.a
#define COLOR_AS_ARGUMENT_NO_ALPHA(C) C.r, C.g, C.b
#define COLOR_AS_ARRAY(C) {C.r, C.g, C.b, C.a}
#define COLOR_AS_ARRAY_NO_ALPHA(C) {C.r, C.g, C.b, 1.0f}

#define F3_ARG(F) F.x, F.y, F.z
#define F3_VTX(F) glVertex3f(F.x, F.y, F.z)

#define TRANSFORM(V, M, TGT) \
	TGT.x = V.x*M[0] + V.y*M[4] + V.z*M[8] + M[12]; \
	TGT.y = V.x*M[1] + V.y*M[5] + V.z*M[9] + M[13]; \
	TGT.z = V.x*M[2] + V.y*M[6] + V.z*M[10] + M[14];

#define RENDER_IMAGE_QUAD(A,B,C,D) \
	glBegin(GL_QUADS); \
	glTexCoord2f(0,0); \
	glVertex3f(A.x, A.y, A.z); \
	glTexCoord2f(1,0); \
	glVertex3f(B.x, B.y, B.z); \
	glTexCoord2f(1,1); \
	glVertex3f(C.x, C.y, C.z); \
	glTexCoord2f(0,1); \
	glVertex3f(D.x, D.y, D.z); \
	glEnd()

#endif
