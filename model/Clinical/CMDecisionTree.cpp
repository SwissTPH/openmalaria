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
#include "PkPd/LSTMTreatments.h"
#include "util/random.h"
#include "Monitoring/Survey.h"
#include "util/errors.h"

#include <limits>
#include <list>
#include <boost/format.hpp>
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/ptr_map.hpp>

namespace OM { namespace Clinical {

using Monitoring::Survey;
using namespace OM::util;
using namespace boost::assign;
using boost::ptr_vector;
using boost::ptr_map;


// ———  special 'multiple' node  ———

/**
 * Branch out to multiple descendants.
 */
class CMDTMultiple : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create( const ::scnXml::DTMultiple& node );
    
protected:
    virtual CMDTOut exec( CMHostData hostData ) const{
        CMDTOut result(false);
        for( ptr_vector<CMDecisionTree>::const_iterator it = children.begin(),
            end = children.end(); it != end; ++it )
        {
            CMDTOut r2 = it->exec( hostData );
            result.treated = result.treated || r2.treated;
        }
        return result;
    }
    
private:
    CMDTMultiple( /*size_t capacity*/ ){
//         children.reserve( capacity );
    }
    
    ptr_vector<CMDecisionTree> children;
};

// ———  branching nodes  ———

/**
 * Should the patient recieve first or second line treatment?
 */
class CMDTCaseType : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create( const ::scnXml::DTCaseType& node );
    
protected:
    virtual CMDTOut exec( CMHostData hostData ) const{
        //TODO: clarify whether or not this works for complicated cases. Currently probably not.
        assert( (hostData.pgState & Episode::SICK) && !(hostData.pgState & Episode::COMPLICATED) );
        
        if( hostData.pgState & Episode::SECOND_CASE ) return secondLine->exec( hostData );
        else return firstLine->exec( hostData );
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
    virtual CMDTOut exec( CMHostData hostData ) const{
        if( isRDT )     Survey::current().report_Clinical_RDTs (1);
        else            Survey::current().report_Clinical_Microscopy (1);
        
        double dens = hostData.withinHost.getTotalDensity ();
        double pPositive = 1.0 + specificity * (dens / (dens + dens_50) - 1.0);
        if( random::bernoulli(pPositive) ) return positive->exec( hostData );
        else return negative->exec( hostData );
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
    
    //NOTE: could be const if constructed differently:
    double dens_50;     // parasite density giving 50% chance of positive outcome
    double specificity; // chance of a negative outcome given no parasites
    bool isRDT;   // for reporting
    
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
    virtual CMDTOut exec( CMHostData hostData ) const{
        ptr_map<double,CMDecisionTree>::const_iterator it =
            branches.upper_bound( random::uniform_01() );
        assert( it != branches.end() );
        return it->second->exec( hostData );
    }
    
private:
    CMDTRandom( ptr_map<double,CMDecisionTree>& branches ){
        this->branches.swap( branches );
    }
    
    // keys are cumulative probabilities; last entry should equal 1
    //NOTE: could be const if initialised differently
    ptr_map<double,CMDecisionTree> branches;
};

// ———  action nodes  ———

/** Do nothing. **/
class CMDTNoAction : public CMDecisionTree {
protected:
    virtual CMDTOut exec( CMHostData hostData ) const{
        return CMDTOut(false);
    }
};

/**
 * Deliver one or more treatments via the PK/PD model.
 */
class CMDTTreatPKPD : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create(
        const ::scnXml::DecisionTree::TreatPKPDSequence& seq );
    
protected:
    virtual CMDTOut exec( CMHostData hostData ) const{
        foreach( const TreatInfo& treatment, treatments ){
            hostData.withinHost.treatPkPd( treatment.schedule, treatment.dosage, hostData.ageYears );
        }
        return CMDTOut(true);
    }
    
private:
    CMDTTreatPKPD( const scnXml::DecisionTree::TreatPKPDSequence& seq ){
        treatments.reserve( seq.size() );
        foreach( const scnXml::DTTreatPKPD& treatElt, seq ){
            treatments.push_back( TreatInfo(
                treatElt.getSchedule(),
                treatElt.getDosage(),
                treatElt.getDelay_h()
            ) );
        }
        assert( treatments.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
    
    struct TreatInfo{
        TreatInfo( const string& s, const string& d, double h ) :
            schedule(PkPd::LSTMTreatments::findSchedule(s)),
            dosage(PkPd::LSTMTreatments::findDosages(d)),
            delay_h(h) {}
        size_t schedule;        // index of the schedule
        size_t dosage;          // name of the dosage table
        double delay_h;         // delay in hours
    };
    vector<TreatInfo> treatments;
};

/**
 * Deliver one or more simple treatments
 */
class CMDTTreatSimple : public CMDecisionTree {
public:
    static auto_ptr<CMDecisionTree> create(
        const ::scnXml::DecisionTree::TreatSimpleSequence& seq );
    
protected:
    virtual CMDTOut exec( CMHostData hostData ) const{
        foreach( const TreatInfo& treatment, treatments ){
            hostData.withinHost.treatSimple( treatment.tsLiver, treatment.tsBlood );
        }
        return CMDTOut(true);
    }
    
private:
    CMDTTreatSimple( const scnXml::DecisionTree::TreatSimpleSequence& seq ){
        treatments.reserve( seq.size() );
        foreach( const scnXml::DTTreatSimple& treatElt, seq ){
            treatments.push_back( TreatInfo(
                treatElt.getTimestepsLiver(), treatElt.getTimestepsBlood()
            ) );
        }
        assert( treatments.size() > 0 );        // CMDTTreatPKPD should not be used in this case
    }
    
    struct TreatInfo{
        TreatInfo( int timestepsLiver, int timestepsBlood ) :
            tsLiver(timestepsLiver), tsBlood(timestepsBlood) {}
        TimeStep tsLiver, tsBlood;
    };
    vector<TreatInfo> treatments;
};


// ———  static functions  ———

auto_ptr<CMDecisionTree> CMDecisionTree::create( const scnXml::DecisionTree& node ){
    // branching nodes
    if( node.getCaseType().present() ) return CMDTCaseType::create( node.getCaseType().get() );
    if( node.getDiagnostic().present() ) return CMDTDiagnostic::create( node.getDiagnostic().get() );
    if( node.getRandom().present() ) return CMDTRandom::create( node.getRandom().get() );
    // action nodes
    if( node.getNoAction().present() ) return auto_ptr<CMDecisionTree>( new CMDTNoAction() );
    if( node.getTreatPKPD().size() ) return CMDTTreatPKPD::create( node.getTreatPKPD() );
    throw xml_scenario_error( "unterminated decision tree" );
}

auto_ptr<CMDecisionTree> CMDTMultiple::create( const scnXml::DTMultiple& node ){
    auto_ptr<CMDTMultiple> self( new CMDTMultiple() );
    foreach( const scnXml::DTCaseType& sn, node.getCaseType() ){
        self->children.push_back( CMDTCaseType::create(sn).release() );
    }
    foreach( const scnXml::DTDiagnostic& sn, node.getDiagnostic() ){
        self->children.push_back( CMDTDiagnostic::create(sn).release() );
    }
    foreach( const scnXml::DTRandom& sn, node.getRandom() ){
        self->children.push_back( CMDTRandom::create(sn).release() );
    }
    if( node.getTreatPKPD().size() ){
        self->children.push_back( CMDTTreatPKPD::create(node.getTreatPKPD()).release() );
    }
    return auto_ptr<CMDecisionTree>( self.release() );
}

auto_ptr<CMDecisionTree> CMDTCaseType::create( const scnXml::DTCaseType& node ){
    return auto_ptr<CMDecisionTree>( new CMDTCaseType(
        CMDecisionTree::create( node.getFirstLine() ),
        CMDecisionTree::create( node.getSecondLine() )
    ) );
}

auto_ptr<CMDecisionTree> CMDTDiagnostic::create( const scnXml::DTDiagnostic& node ){
    return auto_ptr<CMDecisionTree>( new CMDTDiagnostic(
        node.getType(),
        CMDecisionTree::create( node.getPositive() ),
        CMDecisionTree::create( node.getNegative() )
    ) );
}

auto_ptr<CMDecisionTree> CMDTTreatPKPD::create(
    const scnXml::DecisionTree::TreatPKPDSequence& seq )
{
    return auto_ptr<CMDecisionTree>( new CMDTTreatPKPD( seq ) );
}

auto_ptr< CMDecisionTree > CMDTRandom::create(
    const scnXml::DTRandom& node )
{
    ptr_map<double,CMDecisionTree> branches;
    
    double cum_p = 0.0;
    BOOST_FOREACH( const scnXml::Outcome& outcome, node.getOutcome() ){
        cum_p += outcome.getP();
        branches.insert( cum_p, CMDecisionTree::create( outcome ) );
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
    //TODO: set index for last entry to something > 1
//     branches.back()->first = numeric_limits<double>::infinity();

    return auto_ptr<CMDecisionTree>(
        new CMDTRandom( branches ) );
}

} }
