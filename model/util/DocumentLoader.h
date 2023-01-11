/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2022 University of Basel
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

//parse the xml scenario file
//

#ifndef Hmod_util_DocumentLoader
#define Hmod_util_DocumentLoader

#include "Global.h"
#include <schema/scenario.h>
#include <string>
#include <memory>

namespace OM { namespace util {

class DocumentLoader {
public:
    /// Current schema version.
    static const int SCHEMA_VERSION = 45;
    
    DocumentLoader () : documentChanged(false) {}
    
    /** @brief Reads the document in the xmlFile
    * 
    * Throws on failure. */
    void loadDocument(std::string);
    
    /** Save any changes which occurred to the document, if
        * documentChanged is true. */
    void saveDocument();
    
    /** Get the base scenario element.
        *
        * Is an operator for brevity: InputData().getModel()...
        */
    inline const scnXml::Scenario& document() {
        return *scenario;
    }

    /** Get a mutable version of scenario element.
    *
    * This is the only entry point for changing the scenario document.
    * 
    * You should set "documentChanged = true;" if you want your changes saved. */
    inline scnXml::Scenario& getMutableScenario() {
        return *scenario;
    }

    /** Set true if the xml document has been changed and should be saved.
        *
        * Note that the document will be saved between initialisation and
        * running the main simulation, so only changes added during init will
        * be saved. (This avoids worrying about checkpointing.) */
    bool documentChanged;

private:
    /// Sometimes used to save changes to the xml.
    std::string xmlFileName;
    
    /** @brief The xml data structure. */
    unique_ptr<scnXml::Scenario> scenario;
};

} }
#endif
