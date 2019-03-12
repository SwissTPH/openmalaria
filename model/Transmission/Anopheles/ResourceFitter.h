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

#ifndef Hmod_Anopheles_ResourceFitter
#define Hmod_Anopheles_ResourceFitter

#include "Global.h"
#include "Transmission/Anopheles/MosqTransmission.h"
#include "Transmission/Anopheles/LifeCycle.h"
#include <gsl/gsl_vector.h>

namespace OM {
namespace Transmission {
namespace Anopheles {

// forward declarations
// these are private functions called by GSL (which cannot call member functions)
int ResourceFitter_rootfind_sampler( const gsl_vector *x, void *params, gsl_vector *f );
double ResourceFitter_minimise_sampler( const gsl_vector *x, void *params );

/// Container to run life-cycle model with fixed human inputs and run fitting
/// algorithms.
struct ResourceFitter {
    /** Store fixed parameters.
     * 
     * @param mosqTrans MosqTransmission object with parameters
     * @param lcParams Life cycle parameters. invLarvalResources member is set.
     * @param PA Average P_A value (assumed constant)
     * @param Pdf Average P_df value (assumed constant)
     * @param iNvSv Scale factor for S_v to give an initial estimate of N_v
     * @param iOvSv Scale factor for S_v to give an initial estimate of O_v
     * @param forceDebug If true, turn debugging on regardless of command-line
     */
    ResourceFitter( const MosqTransmission& mosqTrans,
                    LifeCycleParams& lcParams,
                    double PA, double Pdf,
                    double iNvSv, double iOvSv,
                    bool forceDebug = false
                  );
    ~ResourceFitter();
    
    // Prevent copying
    ResourceFitter( const ResourceFitter& ) = delete;
    ResourceFitter& operator= ( const ResourceFitter& ) = delete;
    
    // debugging function
    void printState();
    
    /** Set S_v as target and store P_dif needed to calculate S_v.
     *
     * @param S_v Densities of infectious biting mosquitoes over year; has
     *  length 365.
     * @param sampledP_dif P_dif calculated per time-step over the last n years;
     *  has length 365 * n for some integer n.
     */
    void targetS_vWithP_dif( const vector<double>& S_v,
                             const vector<double>& sampledP_dif );
    
    /** Set emergence-rate as target. */
    void targetEmergenceRate( const vector<double>& emergeRate );
    
    /** Run fitting algorithms (root-finding or minimisation). */
    void fit();
    
private:
    enum FitMethod {
        minimise, find_root
    };
    /** Run root-finding/minimisation algorithm.
     *
     * @param order The number of dimensions to try fitting at once. Probably
     * only 365 is valid for root-finding (since otherwise a root may not exist);
     * in any case the last time this is called must be with order 365.
     * @param method Either minimise or find_root.
     * @param maxIter The maximum number of iterations for the algorithm.
     */
    void fit( size_t order, FitMethod method, size_t maxIter );
    
    /** Run captive model. Uses a warmup length of lcParams.getTotalDuration(),
     * then a sampling duration of 365 intervals. */
    void simulate1Year();
    
public: //TODO(vec lifecycle): only public for debugging
    /** Sample: given descriptor for resource availability x, calculate
     * resultant emergence rate.
     * 
     * Function sets buf to (sampled emergence rate) - (target emergence rate).
     * 
     * @param x Periodic descriptor of annual resource availability. If not of
     * length 365 days, vector will be scaled.
     */
    void sampler( const gsl_vector *x );
private:
    
    /** Copy the input to invLarvalResources, with appropriate transformations. */
    void copyToLarvalResources( const gsl_vector* input );
    
    vector<double>& invLarvalResources;
    MosqTransmission transmission;
    double P_A, P_df;
    double initNvFromSv, initOvFromSv;
    enum FitTarget {
        FT_NONE,        // none set yet
        FT_EMERGENCE,
        FT_S_V
    } fitTarget;
    // A vector which is filled with sampled values. Must be of length 1 year;
    // first value corresponds to emergence sampled for 1st day of year.
    std::vector<double> samples;
    std::vector<double> annualP_dif;
    std::vector<double> target;
    gsl_vector *initial_guess;
    gsl_vector *buf;    // working memory, of length 365
    
    friend int ResourceFitter_rootfind_sampler( const gsl_vector *x, void *params, gsl_vector *f );
    friend double ResourceFitter_minimise_sampler( const gsl_vector *x, void *params );
};

}
}
}
#endif
