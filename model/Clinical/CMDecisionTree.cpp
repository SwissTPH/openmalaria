/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#include "Clinical/CMDecisionTree.h"
#include "WithinHost/WHInterface.h"
#include "util/random.h"
#include "Monitoring/Survey.h"
#include "util/errors.h"

#include <limits>
#include <list>
#include <boost/format.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/ptr_container/ptr_map.hpp>

namespace OM { namespace Clinical {

using Monitoring::Survey;
using namespace OM::util;
using namespace boost::assign;
using boost::ptr_map;

/**
 * Should the patient recieve first or second line treatment?
 */
class CMDTCaseType : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create( const ::scnXml::DTCaseType& node );
    
protected:
    virtual void exec( CMHostData hostData ) const{
        //TODO: clarify whether or not this works for complicated cases. Currently probably not.
        assert( (hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED) );
        
        if( hostData.pgState & Episode::SECOND_CASE ) secondLine->exec( hostData );
        else firstLine->exec( hostData );
    }
    
private:
    CMDTCaseType( auto_ptr<CMDecisionTree> firstLine,
                  auto_ptr<CMDecisionTree> secondLine ) :
                  firstLine(firstLine), secondLine(secondLine)
                  {}
    
    const auto_ptr<CMDecisionTree> firstLine;
    const auto_ptr<CMDecisionTree> secondLine;
};

/**
 * Use a diagnostic, and call another decision tree based on the outcome.
 */
class CMDTDiagnostic : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create( const ::scnXml::DTDiagnostic& node );
    
protected:
    virtual void exec( CMHostData hostData ) const{
        if( isRDT )     Survey::current().report_Clinical_RDTs (1);
        else            Survey::current().report_Clinical_Microscopy (1);
        
        double dens = hostData.withinHost.getTotalDensity ();
        double pPositive = 1.0 + specificity * (dens / (dens + dens_50) - 1.0);
        if( random::bernoulli(pPositive) ) positive->exec( hostData );
        else negative->exec( hostData );
    }
    
private:
    CMDTDiagnostic( const string& type,
        auto_ptr<CMDecisionTree> positive,
        auto_ptr<CMDecisionTree> negative ) :
        positive(positive), negative(negative)
    {
        //TODO: these values should not be hard-coded and maybe we should allow
        // more than two types of diagnostic.
        if( type == "microscopy" ){
            // Microscopy sensitivity/specificity data in Africa;
            // Source: expert opinion — Allan Schapira
            dens_50 = 20.0;
            specificity = .75;
            isRDT = false;
        }else{
            assert( type == "RDT" );
            // RDT sensitivity/specificity for Plasmodium falciparum in Africa
            // Source: Murray et al (Clinical Microbiological Reviews, Jan. 2008)
            dens_50 = 50.0;
            specificity = .942;
            isRDT = true;
        }
    }
    
    const double dens_50;     // parasite density giving 50% chance of positive outcome
    const double specificity; // chance of a negative outcome given no parasites
    const bool isRDT;   // for reporting
    
    const auto_ptr<CMDecisionTree> positive;
    const auto_ptr<CMDecisionTree> negative;
};

/**
 * Choose a branch randomly.
 * 
 * TODO: can we re-order the trees such that this can be abstracted out or
 * used only once?
 */
class CMDTRandom : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create( const ::scnXml::DTRandom& node );
    
protected:
    virtual void exec( CMHostData hostData ) const{
        branches.upper_bound( random::uniform_01() )->second->exec( hostData );
    }
    
private:
    CMDTRandom( ptr_map<double,CMDecisionTree>& branches ){
        this.branches.swap( branches );
    }
    
    // keys are cumulative probabilities; last entry should equal 1
    const ptr_map<double,CMDecisionTree> branches;
};


// ———  static functions  ———

auto_ptr<CMDecisionTree> CMDecisionTree::create( const scnXml::DecisionTree& node ){
    if( node.getCaseType().present() ) return CMDTCaseType::create( node.getCaseType().get() );
    if( node.getDiagnostic().present() ) return CMDTDiagnostic::create( node.getDiagnostic().get() );
    if( node.getRandom().present() ) return CMDTRandom::create( node.getRandom().get() );
    assert(false);
}

auto_ptr<CMDecisionTree> CMDTCaseType::create( const scnXml::DTCaseType& node ){
    return new CMDTCaseType(
        CMDecisionTree::create( node.getFirstLine() ),
        CMDecisionTree::create( node.getSecondLine() )
    );
}

auto_ptr<CMDecisionTree> CMDTDiagnostic::create( const scnXml::DTDiagnostic& node ){
    return new CMDTDiagnostic(
        CMDecisionTree::create( node.getPositive() ),
        CMDecisionTree::create( node.getNegative() )
    );
}

auto_ptr< CMDecisionTree > CMDTRandom::create( const scnXml::DTRandom& node ){
    vector<double> cumProbs;
    cumProbs.reserve( node.getOutcome().size() );
    ptr_vector<CMDecisionTree> branches;
    branches.reserve( node.getOutcome().size() );
    
    double cum_p = 0.0;
    BOOST_FOREACH( const scnXml::Outcome& outcome, node.getOutcome() ){
        cum_p += outcome.getP();
        branches.push_back( CMDecisionTree::create( outcome ) );
    }
    
    // Test cum_p is approx. 1.0 in case the input tree is wrong. In any case,
    // we force probabilities to add to 1.0. At least one branch is required!
    if (cum_p < 0.999 || cum_p > 1.001){
        throw util::xml_scenario_error ( (
            boost::format("decision tree (random node): expected probability sum to be 1.0 but found %2%")
            %cum_p
        ).str() );
    }
    // In theory if the last value is exactly 1 it should always be greater
    // than a uniform_01() sample; this just gives us a stronger guarantee of
    // the same thing.
    cumProbs.back() = numeric_limits<double>::infinity();

    return new CMDTRandom( cumProbs, branches );
}

} }
