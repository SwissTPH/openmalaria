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

inline double findAngle(double EIRRotageAngle, const vector<double> & FSCoeffic, vecDay<double> &sim)
{
    vecDay<double> temp(sim.size(), 0.0);

    double delta = 2 * M_PI / 365.0;

    double min = std::numeric_limits<double>::infinity();
    double minAngle = 0.0;
    for(double angle=-M_PI; angle<M_PI; angle+=delta)
    {
        vectors::expIDFT(temp, FSCoeffic, EIRRotageAngle + angle);

        // Minimize l1-norm
        double sum = 0.0;
        for(SimTime i=SimTime::zero(); i<SimTime::oneYear(); i+=SimTime::oneDay())
            sum += fabs(temp[i] - sim[i]);

        if(sum < min)
        {
            min = sum;
            minAngle = angle;
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