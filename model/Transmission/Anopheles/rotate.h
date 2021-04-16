/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2021 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2021 University of Basel
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

#ifndef Hmod_rotate_H
#define Hmod_rotate_H

#include "Global.h"
#include <vector>
#include <limits>

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;

inline int argmax(const vecDay<double> &vec)
{
    int imax = 0;
    double max = 0.0;
    for(int i=0; i<vec.size().inDays(); i++)
    {
        double v = vec[SimTime::fromDays(i)];
        if(v >= max)
        {
            max = v;
            imax = i;
        }
    }
    return imax;
}

inline double findAngle(const double EIRRotageAngle, const vector<double> & FSCoeffic, const vecDay<double> &sim)
{
    vecDay<double> temp(sim.size(), 0.0);

    double delta = 2.0 * M_PI / 365.0;

    double min = std::numeric_limits<double>::infinity();
    double minAngle = 0.0;
    for(double angle=-M_PI; angle<M_PI; angle+=delta)
    {
        vectors::expIDFT(temp, FSCoeffic, EIRRotageAngle + angle);

        // Minimize l1-norm
        double sum = 0.0;
        for(SimTime i=SimTime::zero(); i<SimTime::oneYear(); i+=SimTime::oneDay())
        {
            double v = fabs(temp[i] - sim[i]);
            sum += v*v;
        }

        sum = sqrtf(sum);
        if(sum < min)
        {
            min = sum;
            minAngle = angle;
            // cout << angle << " " << min << " " << sum << endl;
        }

        // Or minimize peaks offset
        // int m1 = argmax(temp);
        // int m2 = argmax(sim);
        // int offset = abs(m1-m2);

        // if(offset < min)
        // {
        //     min = offset;
        //     minAngle = angle;
        // }

    }
    return minAngle;
}

}
}
}

#endif