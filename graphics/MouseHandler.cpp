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


#include "MouseHandler.h"

MouseHandler :: MouseHandler()
: leftDown(false), 
	midDown(false), 
	rightDown(false),
	controller(0)
{}

void MouseHandler :: mouseMoved(int2 abs, int2 rel)
{
	if (!controller) return;

	if (leftDown) controller->rotate(rel);
	else if (rightDown) controller->zoom(rel);
}

void MouseHandler :: mousePressed(int2 where, MouseButton button)
{
	switch (button)
	{
		case LEFT:		leftDown = true;		break;
		case MID:			midDown = true;			break;
		case RIGHT:		rightDown = true;		break;
	}
}

void MouseHandler :: mouseReleased(int2 where, MouseButton button)
{
	switch (button)
	{
		case LEFT:		leftDown = false;		break;
		case MID:			midDown = false;		break;
		case RIGHT:		rightDown = false;	break;
	}
}

void MouseHandler :: mouseButtonUsed(int which, int state)
{
	switch (which)
	{
		case 0:		leftDown = state;		break;
		case 1:		midDown = state;		break;
		case 2:		rightDown = state;	break;
	}
}
