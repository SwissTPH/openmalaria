#ifndef Hmod_rotate_H
#define Hmod_rotate_H

#include "Global.h"
#include <vector>
#include <limits>

namespace OM {
namespace Transmission {
namespace Anopheles {
using namespace OM::util;

inline int argmax(const std::vector<double> &vec)
{
    int imax = 0;
    double max = 0.0;
    for(size_t i=0; i<vec.size(); i++)
    {
        double v = vec[i];
        if(v >= max)
        {
            max = v;
            imax = i;
        }
    }
    return imax;
}

inline double findAngle(const double EIRRotageAngle, const vector<double> & FSCoeffic, const std::vector<double> &sim)
{
    std::vector<double> temp(sim.size(), 0.0);

    double delta = 2.0 * M_PI / 365.0;

    double min = std::numeric_limits<double>::infinity();
    double minAngle = 0.0;
    for(double angle=-M_PI; angle<M_PI; angle+=delta)
    {
        vectors::expIDFT(temp, FSCoeffic, EIRRotageAngle + angle);

        // Minimize l1-norm
        double sum = 0.0;
        for(int i=0; i<SimTime::oneYear().inDays(); i++)
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