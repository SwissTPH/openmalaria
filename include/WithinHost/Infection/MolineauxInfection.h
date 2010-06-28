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
	static void initParameters();
	virtual bool updateDensity(double survivalFactor, int ageOfInfection);

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
	static const int v = 50;
	// taus: used for the variantTranscending and variantSpecific array, 4 Molineaux timesteps = 8 days
	static const int taus = 4;

	///@brief static variables red from parameters
	//@{
	/// mean and sd of the first local maximum density
	static double mean_first_local_max, sd_first_local_max;
	/// mean and sd of the difference between last positive and first positive days
	static double mean_diff_pos_days, sd_diff_pos_days;
	//@}
	
	/** @brief The static variables (double)
	 *
	 * sProb: fraction of parasites switching among variants per two-day cycle
	 * q: Parameter of the geometric distribution of switching probabilities
	 * k_c,k_m: constants allowing calculation of Pstar_c and Pstar_m from host-specific data
	 * Pstar_v: critical density of a variant, common to all variants
	 * kappa_c, kappa_m, kappa_v: Stiffness parameters for saturation of immune responses
	 * C: Maximum daily antigenic stimulus, per mul, of the acquired variant-transcending immune response
	 * sigma, rho: decay parameters, per day, of the acquired variant-specific and variant-transcending immune responses
	 * beta: Minimum value of the probability that a parasite escape control by the acquired and variant-transcending immune response
	 * mu_m, sigma_m: Mean and standard deviation to use for the normal distribution setting the variant specific multiplication factor.
	 */
	//@{
	static const double sigma=0.02;
	static const double rho=0.0;
	static const double beta=0.01;
	static const double sProb=0.02;
	static const double q=0.3;
	static const double mu_m=16.0;
	static const double sigma_m=10.4;
	static const double k_c=0.2;
	static const double k_m=0.04;
	static const double Pstar_v=30.0;
	static const double kappa_c=3.0;
	static const double kappa_m=1.0;
	static const double kappa_v=3.0;
	static const double C=1.0;
	//@}

	/*
	 * The dynamic variables:
	 * m[i]: Multiplication factor, per two-day cycle of variant i
	 * growthRate[i]: variant's i growthRate
	 * P[i]: variant's i density
	 * initP[i]: Density of in t+2 emerging variant i
	 * Pstar_c, Pstar_m: two host-specific critical densities... Those two values depend on the first local maximum or the difference
	 * between the last positive day and the first positive day.
	 * variantTranscendingSummation: See Molineaux paper, equation 7
	 * variantSpecificSummation: See Molineaux paper, equation 6
	 */
	double m[v],variantTranscendingSummation, laggedPc[taus], Pstar_c, Pstar_m;
	struct Variant {
	    double growthRate, P, variantSpecificSummation, initP;
	    double laggedP[taus];
	    
	    Variant ();
	    
	    /// Checkpointing
	    template<class S>
	    void operator& (S& stream) {
		growthRate & stream;
		P & stream;
		variantSpecificSummation & stream;
		initP & stream;
		for(int i = 0; i < taus; ++i)
		    laggedP[i] & stream;
	    }
	    
	    void updateGrowthRateMultiplier( double pd, double immune_response_escape );
	    double updateDensity (double survivalFactor, int ageOfInfection);
	    double getVariantSpecificSummation();
	};
	vector<Variant> variants;
};

}}
#endif
