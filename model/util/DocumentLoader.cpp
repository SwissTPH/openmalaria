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

#include "util/DocumentLoader.h"
#include "util/BoincWrapper.h"
#include "util/errors.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <boost/format.hpp>

namespace OM { namespace util {

Checksum DocumentLoader::loadDocument (std::string lXmlFile){
    xmlFileName = lXmlFile;
    //Parses the document
    
    // Opening by filename causes a schema lookup in the scenario file's dir,
    // which does always work. Opening with a stream uses the working directory.
    
    // Note that the schema location can be set manually by passing properties,
    // but we won't necessarily have the right schema version associated with
    // the XML file in that case.
    ifstream fileStream (lXmlFile.c_str(), ios::binary);
    if (!fileStream.good()){
	string msg = "Error: unable to open "+lXmlFile;
	throw util::xml_scenario_error (msg);
    }
    scenario = scnXml::parseScenario (fileStream);
    util::Checksum cksum = util::Checksum::generate (fileStream);
    fileStream.close ();
    int scenarioVersion = scenario->getSchemaVersion();
    if (scenarioVersion < SCHEMA_VERSION) {
        // Don't bother aborting. Mostly if something really is incompatible
        // loading will not succeed anyway.
        cerr<<"Warning: "<<lXmlFile<<" uses an old schema version (latest is "
            <<SCHEMA_VERSION<<")."<<endl;
    }
    if (scenarioVersion > SCHEMA_VERSION)
        throw util::xml_scenario_error ("Error: new schema version unsupported");
    return cksum;
}

void DocumentLoader::saveDocument()
{
    if (documentChanged) {
        // We don't use this any more, so reduce executable size a bit
        cerr << "Error: document saving support is disabled" << endl;
    }
}


} }
