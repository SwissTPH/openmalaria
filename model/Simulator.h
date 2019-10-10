/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_Simulator
#define Hmod_Simulator

#include "Global.h"
#include "Population.h"
#include "Transmission/TransmissionModel.h"
using namespace std;

namespace scnXml{
    class Monitoring;
    class Scenario;
}
namespace OM {
    
//! Main simulation class
class Simulator{
public: 
    //!  Inititalise all step specific constants and variables.
    Simulator( const scnXml::Scenario& scenario );
    
    //! Entry point to simulation.
    void start(const scnXml::Monitoring& monitoring);
    
    /// Return true when this simulation started by loading a checkpoint
    inline static bool isCheckpoint(){ return startedFromCheckpoint; }
    
private:
    /** @brief checkpointing functions
    *
    * readCheckpoint/writeCheckpoint prepare to read/write the file,
    * and read/write read and write the actual data. */
    //@{
    void writeCheckpoint();
    void readCheckpoint();
    
    void checkpoint (istream& stream);
    void checkpoint (ostream& stream);
    //@}
    
    // Data
    SimTime m_phaseEnd;
    SimTime m_estimatedEnd;
    int phase;  // only need be a class member because value is checkpointed
    
    static bool startedFromCheckpoint;
    
    friend class AnophelesModelSuite;
};

}
#endif
