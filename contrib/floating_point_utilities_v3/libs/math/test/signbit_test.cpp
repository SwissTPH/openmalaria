// signbit_test.cpp

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
#include "../../../boost/math/signbit.hpp"

namespace {

// the anonymous namespace resolves ambiguities on platforms
// with fpclassify etc functions at global scope

using boost::math::signbit;
using boost::math::copysign;
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
    T x, y;

    x = static_cast<T>(0);
    BOOST_CHECK(!(signbit)(x));

    x = (changesign)(static_cast<T>(0));
    BOOST_CHECK(x == 0);
    BOOST_CHECK((signbit)(x));

    if(std::numeric_limits<T>::has_denorm) {

        x = std::numeric_limits<T>::denorm_min();
        BOOST_CHECK(!(signbit)(x));

        x = -std::numeric_limits<T>::denorm_min();
        // some platforms may round off subnormals to 0
        if(x != 0)
            BOOST_CHECK((signbit)(x));
    }

    x = (std::numeric_limits<T>::min)();
    BOOST_CHECK(!(signbit)(x));

    x = -(std::numeric_limits<T>::min)();
    BOOST_CHECK((signbit)(x));

    x = static_cast<T>(1);
    BOOST_CHECK(!(signbit)(x));

    x = static_cast<T>(-1);
    BOOST_CHECK((signbit)(x));

    x = static_cast<T>(123.456);
    BOOST_CHECK(!(signbit)(x));

    x = static_cast<T>(-123.456);
    BOOST_CHECK((signbit)(x));
    
    x = static_cast<T>(3407);
    BOOST_CHECK(!(signbit)(x));
    BOOST_CHECK((changesign)(x) == -3407);
    BOOST_CHECK((copysign)(x, 1) == 3407);
    BOOST_CHECK((copysign)(x, -1) == -3407);
    
    x = static_cast<T>(-3407);
    BOOST_CHECK((signbit)(x));
    BOOST_CHECK((changesign)(x) == 3407);
    BOOST_CHECK((copysign)(x, 1) == 3407);
    BOOST_CHECK((copysign)(x, -1) == -3407);

    x = (std::numeric_limits<T>::max)();
    BOOST_CHECK(!(signbit)(x));

    x = -(std::numeric_limits<T>::max)();
    BOOST_CHECK((signbit)(x));
    
    if(std::numeric_limits<T>::has_infinity) {

        x = std::numeric_limits<T>::infinity();
        BOOST_CHECK(!(signbit)(x));

        x = (changesign)(std::numeric_limits<T>::infinity());
        BOOST_CHECK((signbit)(x));
    }

    if(std::numeric_limits<T>::has_quiet_NaN) {

        x = std::numeric_limits<T>::quiet_NaN();
        y = (changesign)(x);
        BOOST_CHECK((signbit)(x) ^ (signbit)(y));
    }

    if(std::numeric_limits<T>::has_signaling_NaN) {

        x = std::numeric_limits<T>::signaling_NaN();
        y = (changesign)(x);
        BOOST_CHECK((signbit)(x) ^ (signbit)(y));
    }
}

//------------------------------------------------------------------------------

}   // namespace
