/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#ifndef Hmod_Output_Continuous
#define Hmod_Output_Continuous

#include "Global.h"
#include <FastDelegate.h>
#include <string>

using fastdelegate::MakeDelegate;

namespace OM { namespace Monitoring {
    
    /** Class to deal with continuous output data.
     *
     * Requirements:
     *  (1) frequency of and which data is output should be controllable
     *  (2) format should be compatible with LiveGraph and (German) Excel.
     */
    class Continuous {
    public:
	/** Load XML description of options. If resuming from a checkpoint,
	 * append to output; if not, make sure it's not there (on boinc we
	 * assume we shouldn't overwrite existing files for security reasons).
	 * 
	 * Callbacks should be registered before init() is called. */
	static void init (bool isCheckpoint);
        
        /** When compiled in BOINC mode, this copies data to the final
         * compressing output file. Otherwise it does nothing. */
        static void finalise ();
	
	/** @brief Checkpointing functions */
	//@{
	static void staticCheckpoint (ostream& stream);
	static void staticCheckpoint (istream& stream);
	//@}
	
	/** Register a callback function which produces output.
	 *
	 * This function will be called to generate output, if enabled in XML.
	 * It may output more than one statistic, if for example vector output
	 * is wanted instead of a single value. It should then title these in
	 * the form "name(index)".
	 *
	 * @param optName	Name of this output, (used for XML on/off options)
	 * @param titles	Titles for the output table; each should be preceeded by a \t
	 * @param outputCb A callback function, which when called, outputs its
	 * data to the passed stream, with each entry preceeded by '\t'.
	 */
	static void registerCallback (string optName, string titles, fastdelegate::FastDelegate1<ostream&>);
	
	/// Generate time-step's output. Called at beginning of timestep.
	static void update ();
    };
} }
#endif