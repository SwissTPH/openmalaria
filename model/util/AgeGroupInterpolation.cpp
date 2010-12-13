/*
 This file is part of OpenMalaria.

 Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "util/AgeGroupInterpolation.h"
#include "util/errors.h"
#include "boost/format.hpp"

namespace OM { namespace util {
    

    /** Dummy object, for initialization. Potentially makes code safer. */
    class AgeGroupDummy : public AgeGroupInterpolation
    {
    public:
        static AgeGroupDummy singleton;
        
        virtual double operator() (double ageYears) const {
            throw logic_error( "AgeGroupDummy" );
        }
        virtual void scale( double factor ) {
            throw logic_error( "AgeGroupDummy" );
        }
        
    protected:
        virtual void checkpoint (ostream& stream) {}
        virtual void checkpoint (istream& stream) {}
    };
    AgeGroupDummy AgeGroupDummy::singleton;

    
    AgeGroupInterpolation* AgeGroupInterpolation::dummyObject(){
        return &AgeGroupDummy::singleton;
    }
    AgeGroupInterpolation* AgeGroupInterpolation::makeObject(
        const scnXml::AgeGroupValues& ageGroups, const char* eltName
    ){
        // Type mostly equivalent to a std::string:
        const scnXml::Interpolation& interp = ageGroups.getInterpolation();
        if( interp == "linear" ){
            return new AgeGroupPiecewiseLinear( ageGroups, eltName );
        }else if( interp == "none" ){
            return new AgeGroupPiecewiseConstant( ageGroups, eltName );
        }else{
            throw util::xml_scenario_error( (boost::format( "age group interpolation %1% not implemented" ) %interp).str() );
        }
    }
    void AgeGroupInterpolation::freeObject( AgeGroupInterpolation* obj ){
        assert( obj != NULL );  // should not do that
        if( obj != &AgeGroupDummy::singleton ){
            delete obj;
            obj = &AgeGroupDummy::singleton;
        }
    }


    AgeGroupPiecewiseConstant::AgeGroupPiecewiseConstant(
        const scnXml::AgeGroupValues& ageGroups,
        const char* eltName
    ){
        if( ageGroups.getGroup().size() == 0 ){
            throw util::xml_scenario_error( (
                boost::format( "%1%: at least one age group required" ) %eltName
            ).str() );
        }
        
        map<double,double>::iterator pos = dataPoints.begin();
        assert(pos == dataPoints.end());
        
        double greatestLbound = -1.0;
        BOOST_FOREACH( const scnXml::Group& group, ageGroups.getGroup() ){
            double lbound = group.getLowerbound();
            if( lbound >= greatestLbound ){
                greatestLbound = lbound;
                pos = dataPoints.insert( pos, make_pair( lbound, group.getValue() ) );
            } else {
                throw util::xml_scenario_error( (
                    boost::format("%3%: lower bound %1% not greater than previous %2%")
                    %lbound %greatestLbound %eltName
                ).str() );
            }
        }
        
        // first lower-bound must be 0
        if( dataPoints.begin()->first != 0.0 ){
            throw util::xml_scenario_error( (
                boost::format("%1%: first lower-bound must be 0")
                %eltName
            ).str() );
        }
    }
    
    double AgeGroupPiecewiseConstant::operator() (double ageYears) const {
        assert ( ageYears >= 0.0 );
        map<double,double>::const_iterator it = dataPoints.upper_bound( ageYears );
        --it;   // previous point when ordered by age
        return it->second;
    }
    
    void AgeGroupPiecewiseConstant::scale( double factor ){
        for( map<double,double>::iterator it = dataPoints.begin(); it != dataPoints.end(); ++it ){
            it->second *= factor;
        }
    }
    
    void AgeGroupPiecewiseConstant::checkpoint (ostream& stream){
        dataPoints & stream;
    }
    void AgeGroupPiecewiseConstant::checkpoint (istream& stream){
        dataPoints & stream;
    }

    
    AgeGroupPiecewiseLinear::AgeGroupPiecewiseLinear(
        const scnXml::AgeGroupValues& ageGroups,
        const char* eltName
    ){
        // Our read iterator
        typedef scnXml::AgeGroupValues::GroupConstIterator GroupIt;
        GroupIt it = ageGroups.getGroup().begin();
        
        if( it == ageGroups.getGroup().end() ){
            throw util::xml_scenario_error( (
                boost::format( "%1%: at least one age group required" ) %eltName
            ).str() );
        }
        // first lower-bound must be 0
        if( it->getLowerbound() != 0.0 ){
            throw util::xml_scenario_error( (
                boost::format("%1%: first lower-bound must be 0")
                %eltName
            ).str() );
        }
        
        // Our insert iterator (used for performance)
        map<double,double>::iterator pos = dataPoints.begin();
        assert(pos == dataPoints.end());
        
        // Add first point at zero
        pos = dataPoints.insert( pos, make_pair( 0.0, it->getValue() ) );
        
        double greatestLbound = 0.0;    // also last lbound: distribution must be weakly monotonic
        for( ; it!=ageGroups.getGroup().end(); ++it ){
            double lbound = it->getLowerbound();
            if( lbound >= greatestLbound ){
                double insBound = 0.5*(greatestLbound + lbound);
                pos = dataPoints.insert( pos, make_pair( insBound, it->getValue() ) );
                greatestLbound = lbound;
            } else {
                throw util::xml_scenario_error( (
                    boost::format("%3%: lower bound %1% less than previous %2%")
                    %lbound %greatestLbound %eltName
                ).str() );
            }
        }
        
        dataPoints[ numeric_limits<double>::infinity() ] = dataPoints.rbegin()->second;
    }
    
    double AgeGroupPiecewiseLinear::operator() (double ageYears) const {
        assert ( ageYears >= 0.0 );
        // Get the first point with age greater than ageYears:
        map<double,double>::const_iterator it = dataPoints.upper_bound( ageYears );
        assert( it != dataPoints.end() );
        double a1 = it->first;  // a1 > ageYears
        double f1 = it->second;
        --it;   // previous point when ordered by age
        double a0 = it->first;  // a0 ≤ ageYears
        double f0 = it->second;
        return (ageYears - a0) / (a1 - a0) * (f1 - f0) + f0;
    }
    
    void AgeGroupPiecewiseLinear::scale( double factor ){
        for( map<double,double>::iterator it = dataPoints.begin(); it != dataPoints.end(); ++it ){
            it->second *= factor;
        }
    }
    
    void AgeGroupPiecewiseLinear::checkpoint (ostream& stream){
        dataPoints & stream;
    }
    void AgeGroupPiecewiseLinear::checkpoint (istream& stream){
        dataPoints & stream;
    }
} }
