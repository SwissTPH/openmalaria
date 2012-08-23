/*
This file is part of OpenMalaria.

Copyright (C) 2005-2010 Swiss Tropical Institute and Liverpool School Of Tropical Medicine

OpenMalaria is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*This model refers to the paper:
 *L. MOLINEAUX, H. H. DIEBNER, M. EICHNER, W. E. COLLINS, G. M. JEFFERY and K. DIETZ (2001). Plasmodium falciparum parasitaemia described by a new mathematical model. Parasitology, 122 , pp 379-391
 *doi:10.1017/S0031182001007533*/

#ifndef Hmod_MOLINEAUXINFECTION_H
#define Hmod_MOLINEAUXINFECTION_H

#include "WithinHost/Infection/CommonInfection.h"

namespace OM { namespace WithinHost {

class MolineauxInfection : public CommonInfection {

public:
	MolineauxInfection (istream& stream);
	    //! Constructor
	MolineauxInfection(uint32_t protID);

	virtual ~MolineauxInfection () {};
	static void init();
	virtual bool updateDensity(double survivalFactor, TimeStep ageOfInfection);

protected:
    virtual void checkpoint (ostream& stream);

private:
    double getVariantSpecificSummation(int i, double P_current);
    double getVariantTranscendingSummation();

    /**This function adapt the growth rate.
     * We can't use the Molineaux as it, since
     * this model is a two day timestep model.
     * the density p(t+1) is then extrapolated.
     *
     */
    void updateGrowthRateMultiplier();

	// v: number of variants per clone (one infection = one new clone)
	static const size_t v = 50;
	// taus: used for the variantTranscending and variantSpecific array, 4 Molineaux timesteps = 8 days
	static const size_t taus = 4;

	///@brief static variables red from parameters
	//@{
	/// depends on lognormal/gamma distribution for first_local_max
	static bool first_local_maximum_gamma;
	/// mean/shape and sd/scale of the first local maximum density for lognormal/gamma distribution
	static double mean_shape_first_local_max, sd_scale_first_local_max;
	
	/// depends on between lognormal/gamma distribution for mean_duration diff_pos_days
	static bool mean_duration_gamma;
	/// mean/shape and sd/scale of the difference between last positive and first positive days for lognormal/gamma distribution
	static double mean_shape_diff_pos_days, sd_scale_diff_pos_days;
	
	/// boolean choosing between gamma and lognormal distribution for equation 11
	static bool multi_factor_gamma;
	//@}
	
	/** @brief q^(i+1) array
	 *
	 * all the values of q^1... q^50 are stored in this array.
	 * this prevent the recalculation of those values on every two timesteps. */
	static double qPow[v];
	
	// m[i]: Multiplication factor, per two-day cycle of variant i
	float m[v];
        // variantTranscendingSummation: See Molineaux paper, equation 7
        float variantTranscendingSummation;
        float laggedPc[taus];
        /* Pstar_c, Pstar_m: two host-specific critical densities...
         * Those two values depend on the first local maximum or the difference
         * between the last positive day and the first positive day. */
        float Pstar_c, Pstar_m;
	struct Variant {
            // growthRate[i]: variant's i growthRate
            float growthRate;
            // P[i]: variant's i density
            float P;
            // variantSpecificSummation: See Molineaux paper, equation 6
            float variantSpecificSummation;
            // initP[i]: Density of in t+2 emerging variant i
            float initP;
	    float laggedP[taus];
	    
	    Variant ();
	    
	    /// Checkpointing
	    void operator& (ostream& stream);
	    void operator& (istream& stream);
	    
	    void updateGrowthRateMultiplier( double pd, double immune_response_escape );
	    double updateDensity (double survivalFactor, TimeStep ageOfInfection);
	    double getVariantSpecificSummation();
	};
	vector<Variant> variants;
};

}}
#endif
