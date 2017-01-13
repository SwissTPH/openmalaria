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

#ifndef Hmod_WithinHost_Vivax
#define Hmod_WithinHost_Vivax

#include "Global.h"
#include "WithinHost/WHInterface.h"

#include <list>
#include <memory>

using namespace std;

class UnittestUtil;

namespace scnXml{
    class LiverStageDrug;
}
namespace OM {
namespace WithinHost {
namespace Pathogenesis {
    class PathogenesisModel;
}

class WHVivax;
/**
 * A brood is the set of hypnozoites resulting from an innoculation, plus an
 * associated blood stage.
 * 
 * In this model, if a hypnozoite releases while a blood stage infection
 * initiated by another hypnozoite from the same brood is active, the newly
 * released hypnozoite does nothing, however, blood stage infections from other
 * broods have no effect.
 */
class VivaxBrood{
public:
    /** Create.
     * 
     * @param host      The host creating this. Not really needed, except to
     * prevent VivaxBrood being default-constructed by its container.
     */
    VivaxBrood( WHVivax *host );
    ~VivaxBrood();
    /** Save a checkpoint. */
    void checkpoint( ostream& stream );
    /** Create from checkpoint. */
    VivaxBrood( istream& stream );
    
    struct UpdResult{
        UpdResult() : newPrimaryBS(false), newRelapseBS(false), newBS(false) {}
        bool newPrimaryBS, newRelapseBS, newBS, isFinished;
    };
    /**
     * Do per time step update: remove finished blood stage infections and act
     * on newly releasing hypnozoites.
     * 
     * @return pair describing whether the infection is finished (no more blood
     *  or liver stages) and whether a new primary blood-stage has started (to
     *  update cumPrimInf)
     */
    UpdResult update();
    
    inline void setHadEvent( bool hadEvent ){ this->hadEvent = hadEvent; }
    inline bool hasHadEvent()const{ return hadEvent; }
    inline void setHadRelapse( bool hadRelapse ){ this->hadRelapse = hadRelapse; }
    inline bool hasHadRelapse()const{ return hadRelapse; }
    
    /** Equivalent to a blood stage existing. We do not model incidence of
     * gametocytes independently, thus this also tests existance of
     * gametocytes. */
    inline bool isPatent() const{
        return bloodStageClearDate > sim::latestTs0();
    }
    
    /** Fully clear blood stage parasites. */
    void treatmentBS();
    
    /** Fully clear liver stage parasites. */
    void treatmentLS();
    
private:
    VivaxBrood() {}     // not default constructible
    
    // list of times at which the merozoite and hypnozoites release, ordered by
    // time of release, soonest last (i.e. last element is next one to release)
    vector<SimTime> releaseDates;
    
    // Either SimTime::never() (no blood stage) or the start of the time step on
    // which the blood stage will clear.
    SimTime bloodStageClearDate;
    
    // Whether the primary blood stage infection has started
    bool primaryHasStarted;

    // Whether the relapse blood stage infection has started
    bool relapseHasStarted;

    // Whether any clinical event has been triggered by this brood
    bool hadEvent;

    // Whether a relapse event has been triggered by this brood
    bool hadRelapse;
};

/**
 * Implementation of a very basic vivax model.
 * 
 * This is for tropical Vivax (low transmission) and where there is little
 * immunity.
 */
class WHVivax : public WHInterface {
public:
    /// @brief Static methods
    //@{
    /// Initialise static parameters
    static void init( const OM::Parameters& parameters, const scnXml::Model& model );
    
    /** Called when health system parameters are loaded.
     * 
     * If no "LiverStageDrug" parameters are present, this is still called but with null pointer.
     */
    static void setHSParameters( const scnXml::LiverStageDrug* );
    //@}

    /// @brief Constructors, destructors and checkpointing functions
    //@{
    WHVivax( double comorbidityFactor );
    virtual ~WHVivax();
    //@}
    
    virtual double probTransmissionToMosquito( double tbvFactor, double *sumX )const;
    virtual double pTransGenotype( double pTrans, double sumX, size_t genotype );
    
    virtual bool summarize(const Host::Human& human) const;
    
    virtual void importInfection();
    
    virtual void update(int nNewInfs, vector<double>& genotype_weights,
            double ageInYears, double bsvFactor);
    
    virtual bool diagnosticResult( const Diagnostic& diagnostic ) const;

    virtual Pathogenesis::StatePair determineMorbidity( Host::Human& human, double ageYears, bool );
    
    virtual void clearImmunity();
    
protected:
    virtual void treatment( Host::Human& human, TreatmentId treatId );
    virtual bool treatSimple( const Host::Human& human, SimTime timeLiver, SimTime timeBlood );
    virtual void optionalPqTreatment( const Host::Human& human );
    
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
    // None of these do anything in this model:
    virtual void treatPkPd(size_t schedule, size_t dosages, double age);
    virtual double getTotalDensity() const;
    virtual double getCumulative_h() const;
    virtual double getCumulative_Y() const;
    
private:
    WHVivax( const WHVivax& ) {}        // not copy constructible
    
    list<VivaxBrood> infections;
    
    /* Is flagged as never getting PQ: this is a heteogeneity factor. Example:
     * Set to zero if everyone can get PQ, 0.5 if females can't get PQ and
     * males aren't tested (i.e. all can get it) or (1+p)/2 where p is the
     * chance a male being tested and found to be G6PD deficient. */
    bool noPQ;
    
    Pathogenesis::State morbidity;
    
    // The number of primary blood stage infections to have started since this
    // human was born.
    uint32_t cumPrimInf;
    
    // Either never or expiry time of treatment
    SimTime treatExpiryLiver, treatExpiryBlood;

    // The probability of a clinical event from a first infection
    double pEvent;
    // The probability of a clinical event from a relapse
    double pFirstRelapseEvent;
    
    // Used for reporting; updated by update()
    double pSevere;

    friend class ::UnittestUtil;
};

}
}
#endif
