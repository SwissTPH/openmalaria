// fpclassify_test.cpp

// Copyright (c) 2006 Johan Rade

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifdef _MSC_VER
#   pragma warning(disable: 4127 4702)
#endif

//-------------------------------------

#include <limits>
#include <boost/test/auto_unit_test.hpp>
#include "../../../boost/math/fpclassify.hpp"
#include "../../../boost/math/signbit.hpp"

namespace {

// the anonymous namespace resolves ambiguities on platforms
// with fpclassify etc functions at global scope

using boost::math::fpclassify;
using boost::math::isfinite;
using boost::math::isnormal;
using boost::math::isinf;
using boost::math::isnan;

using boost::math::changesign;

//------------------------------------------------------------------------------

template<class T> void test();

//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(fpclassify_test)
{
    test<float>();
    test<double>();
    test<long double>();
}

//------------------------------------------------------------------------------

template<class T> void test()
{
    T x;

    x = static_cast<T>(0);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK(!(isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_ZERO);

    x = (changesign)(static_cast<T>(0));
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK(!(isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_ZERO);

    if(std::numeric_limits<T>::has_denorm) {

        x = std::numeric_limits<T>::denorm_min();
        BOOST_CHECK((isfinite)(x));
        BOOST_CHECK(!(isnormal)(x));
        BOOST_CHECK(!(isinf)(x));
        BOOST_CHECK(!(isnan)(x));
        BOOST_CHECK((fpclassify)(x) == FP_SUBNORMAL);

        x = -std::numeric_limits<T>::denorm_min();
        // some platforms may round off subnormals to 0
        if(x != 0) {
            BOOST_CHECK((isfinite)(x));
            BOOST_CHECK(!(isnormal)(x));
            BOOST_CHECK(!(isinf)(x));
            BOOST_CHECK(!(isnan)(x));
            BOOST_CHECK((fpclassify)(x) == FP_SUBNORMAL);
        }
    }

    // some platforms, for instance Apple / GCC 4.0 / PowerPC,
    // have a broken std::numeric_limits<long double>::min
    if(sizeof(T) <= 8 || (std::numeric_limits<long double>::min)()
            < 1e15 * (std::numeric_limits<double>::min())) {

        x = (std::numeric_limits<T>::min)() / 2;
        // some platforms may round off subnormals to 0
        if(x != 0) {
            BOOST_CHECK((isfinite)(x));
            BOOST_CHECK(!(isnormal)(x));
            BOOST_CHECK(!(isinf)(x));
            BOOST_CHECK(!(isnan)(x));
            BOOST_CHECK((fpclassify)(x) == FP_SUBNORMAL);
        }

        x = -(std::numeric_limits<T>::min)() / 2;
        // some platforms may round off subnormals to 0
        if(x != 0) {
            BOOST_CHECK((isfinite)(x));
            BOOST_CHECK(!(isnormal)(x));
            BOOST_CHECK(!(isinf)(x));
            BOOST_CHECK(!(isnan)(x));
            BOOST_CHECK((fpclassify)(x) == FP_SUBNORMAL);
        }
    }

    x = (std::numeric_limits<T>::min)();
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = -(std::numeric_limits<T>::min)();
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    //..........................................................................

    x = static_cast<T>(1);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = static_cast<T>(-1);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = static_cast<T>(123.456);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = static_cast<T>(-123.456);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = static_cast<T>(3407);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = static_cast<T>(-3407);
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = (std::numeric_limits<T>::max)();
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    x = -(std::numeric_limits<T>::max)();
    BOOST_CHECK((isfinite)(x));
    BOOST_CHECK((isnormal)(x));
    BOOST_CHECK(!(isinf)(x));
    BOOST_CHECK(!(isnan)(x));
    BOOST_CHECK((fpclassify)(x) == FP_NORMAL);

    if(std::numeric_limits<T>::has_infinity) {

        x = std::numeric_limits<T>::infinity();
        BOOST_CHECK(!(isfinite)(x));
        BOOST_CHECK(!(isnormal)(x));
        BOOST_CHECK((isinf)(x));
        BOOST_CHECK(!(isnan)(x));
        BOOST_CHECK((fpclassify)(x) == FP_INFINITE);

        x = (changesign)(std::numeric_limits<T>::infinity());
        BOOST_CHECK(!(isfinite)(x));
        BOOST_CHECK(!(isnormal)(x));
        BOOST_CHECK((isinf)(x));
        BOOST_CHECK(!(isnan)(x));
        BOOST_CHECK((fpclassify)(x) == FP_INFINITE);
    }

    if(std::numeric_limits<T>::has_quiet_NaN) {

        x = std::numeric_limits<T>::quiet_NaN();
        BOOST_CHECK(!(isfinite)(x));
        BOOST_CHECK(!(isnormal)(x));
        BOOST_CHECK(!(isinf)(x));
        BOOST_CHECK((isnan)(x));
        BOOST_CHECK((fpclassify)(x) == FP_NAN);

        x = (changesign)(std::numeric_limits<T>::quiet_NaN());
        BOOST_CHECK(!(isfinite)(x));
        BOOST_CHECK(!(isnormal)(x));
        BOOST_CHECK(!(isinf)(x));
        BOOST_CHECK((isnan)(x));
        BOOST_CHECK((fpclassify)(x) == FP_NAN);
    }

    if(std::numeric_limits<T>::has_signaling_NaN) {

        // Intel 7 and VC 6 have broken numeric_limits<T>::signaling_NaN
        if(std::numeric_limits<T>::signaling_NaN()
                != -std::numeric_limits<T>::infinity()) {

            x = std::numeric_limits<T>::signaling_NaN();
            BOOST_CHECK(!(isfinite)(x));
            BOOST_CHECK(!(isnormal)(x));
            BOOST_CHECK(!(isinf)(x));
            BOOST_CHECK((isnan)(x));
            BOOST_CHECK((fpclassify)(x) == FP_NAN);
            
            x = (changesign)(std::numeric_limits<T>::signaling_NaN());
            BOOST_CHECK(!(isfinite)(x));
            BOOST_CHECK(!(isnormal)(x));
            BOOST_CHECK(!(isinf)(x));
            BOOST_CHECK((isnan)(x));
            BOOST_CHECK((fpclassify)(x) == FP_NAN);
        }
    }
}

//------------------------------------------------------------------------------

}   // namespace
