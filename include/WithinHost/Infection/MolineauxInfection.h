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

	static const int v = 50;
	static const int taus = 4;

	static int delta;
	static double C, sigma, rho, beta, sProb, q, mu_m, sigma_m, k_c, k_m, Pstar_v, kappa_c, kappa_m, kappa_v;
	double m[v],variantTranscendingSummation, growthRate[v], P[v], variantSpecificSummation[v], laggedP[taus][v], laggedPc[taus], Pstar_c, Pstar_m, initP[v];
};

}}
#endif
