/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "util/XmlUtils.h"
#include "util/errors.h"
#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace OM { namespace util { namespace XmlUtils {

vector<double> DoubleList2std (const scnXml::DoubleList& list, size_t length) {
  const scnXml::DoubleList::ItemSequence seq = list.getItem();
  if (seq.size() != length)
    throw xml_scenario_error ("readDoubleList: XML list has wrong length");
  vector<double> ret (length);
  for (size_t i = 0; i < length; ++i)
    ret[i] = seq[i];
  return ret;
}

void lboundGroups2map (
    map<double,double>& outp,
    const scnXml::AgeGroupValues::GroupSequence& ageGroups,
    const char* eltName,
    bool addUbound
){
    map<double,double>::iterator pos = outp.begin();
    assert(pos == outp.end());
    
    double greatestLbound = -1.0;
    BOOST_FOREACH( const scnXml::Group& group, ageGroups ){
	double lbound = group.getLowerbound();
	if( lbound > greatestLbound ){
	    greatestLbound = lbound;
	    pos = outp.insert( pos, make_pair( lbound, group.getValue() ) );
	} else {
	    throw util::xml_scenario_error( (
		boost::format("%3%: lower bound %1% not greater than previous %2%")
		%lbound %greatestLbound %eltName
	    ).str() );
	}
    }
    
    // first lower-bound must be 0
    if( outp.begin()->first != 0.0 ){
	throw util::xml_scenario_error( (
	    boost::format("%1%: first lower-bound must be 0")
	    %eltName
	).str() );
    }
    
    if( addUbound ){
	// Useful where interpolation is used.
	outp[ numeric_limits<double>::infinity() ] =
	    outp.rbegin()->second;
    }
}

} } }
