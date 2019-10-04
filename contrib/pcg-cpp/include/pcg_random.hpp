/*
 * PCG Random Number Generation for C++
 *
 * Copyright 2014-2019 Melissa O'Neill <oneill@pcg-random.org>,
 *                     and the PCG Project contributors.
 *
 * SPDX-License-Identifier: (Apache-2.0 OR MIT)
 *
 * Licensed under the Apache License, Version 2.0 (provided in
 * LICENSE-APACHE.txt and at http://www.apache.org/licenses/LICENSE-2.0)
 * or under the MIT license (provided in LICENSE-MIT.txt and at
 * http://opensource.org/licenses/MIT), at your option. This file may not
 * be copied, modified, or distributed except according to those terms.
 *
 * Distributed on an "AS IS" BASIS, WITHOUT WARRANTY OF ANY KIND, either
 * express or implied.  See your chosen license for details.
 *
 * For additional information about the PCG random number generation scheme,
 * visit http://www.pcg-random.org/.
 */

/*
 * This code provides the reference implementation of the PCG family of
 * random number generators.  The code is complex because it implements
 *
 *      - several members of the PCG family, specifically members corresponding
 *        to the output functions:
 *             - XSH RR         (good for 64-bit state, 32-bit output)
 *             - XSH RS         (good for 64-bit state, 32-bit output)
 *             - XSL RR         (good for 128-bit state, 64-bit output)
 *             - RXS M XS       (statistically most powerful generator)
 *             - XSL RR RR      (good for 128-bit state, 128-bit output)
 *             - and RXS, RXS M, XSH, XSL       (mostly for testing)
 *      - at potentially *arbitrary* bit sizes
 *      - with four different techniques for random streams (MCG, one-stream
 *        LCG, settable-stream LCG, unique-stream LCG)
 *      - and the extended generation schemes allowing arbitrary periods
 *      - with all features of C++11 random number generation (and more),
 *        some of which are somewhat painful, including
 *            - initializing with a SeedSequence which writes 32-bit values
 *              to memory, even though the state of the generator may not
 *              use 32-bit values (it might use smaller or larger integers)
 *            - I/O for RNGs and a prescribed format, which needs to handle
 *              the issue that 8-bit and 128-bit integers don't have working
 *              I/O routines (e.g., normally 8-bit = char, not integer)
 *            - equality and inequality for RNGs
 *      - and a number of convenience typedefs to mask all the complexity
 *
 * The code employes a fairly heavy level of abstraction, and has to deal
 * with various C++ minutia.  If you're looking to learn about how the PCG
 * scheme works, you're probably best of starting with one of the other
 * codebases (see www.pcg-random.org).  But if you're curious about the
 * constants for the various output functions used in those other, simpler,
 * codebases, this code shows how they are calculated.
 *
 * On the positive side, at least there are convenience typedefs so that you
 * can say
 *
 *      pcg32 myRNG;
 *
 * rather than:
 *
 *      pcg_detail::engine<
 *          uint32_t,                                           // Output Type
 *          uint64_t,                                           // State Type
 *          pcg_detail::xsh_rr_mixin<uint32_t, uint64_t>, true, // Output Func
 *          pcg_detail::specific_stream<uint64_t>,              // Stream Kind
 *          pcg_detail::default_multiplier<uint64_t>            // LCG Mult
 *      > myRNG;
 *
 */

#ifndef PCG_RAND_HPP_INCLUDED
#define PCG_RAND_HPP_INCLUDED 1

#include <algorithm>
#include <cinttypes>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <utility>
#include <locale>
#include <new>
#include <stdexcept>

#ifdef _MSC_VER
    #pragma warning(disable:4146)
#endif

#ifdef _MSC_VER
    #define PCG_ALWAYS_INLINE __forceinline
#elif __GNUC__
    #define PCG_ALWAYS_INLINE __attribute__((always_inline))
#else
    #define PCG_ALWAYS_INLINE inline
#endif

/*
 * The pcg_extras namespace contains some support code that is likley to
 * be useful for a variety of RNGs, including:
 *      - 128-bit int support for platforms where it isn't available natively
 *      - bit twiddling operations
 *      - I/O of 128-bit and 8-bit integers
 *      - Handling the evilness of SeedSeq
 *      - Support for efficiently producing random numbers less than a given
 *        bound
 */

#include "pcg_extras.hpp"

namespace pcg_detail {

using namespace pcg_extras;

/*
 * The LCG generators need some constants to function.  This code lets you
 * look up the constant by *type*.  For example
 *
 *      default_multiplier<uint32_t>::multiplier()
 *
 * gives you the default multipler for 32-bit integers.  We use the name
 * of the constant and not a generic word like value to allow these classes
 * to be used as mixins.
 */

template <typename T>
struct default_multiplier {
    // Not defined for an arbitrary type
};

template <typename T>
struct default_increment {
    // Not defined for an arbitrary type
};

#define PCG_DEFINE_CONSTANT(type, what, kind, constant) \
        template <>                                     \
        struct what ## _ ## kind<type> {                \
            static constexpr type kind() {              \
                return constant;                        \
            }                                           \
        };

PCG_DEFINE_CONSTANT(uint8_t,  default, multiplier, 141U)
PCG_DEFINE_CONSTANT(uint8_t,  default, increment,  77U)

PCG_DEFINE_CONSTANT(uint16_t, default, multiplier, 12829U)
PCG_DEFINE_CONSTANT(uint16_t, default, increment,  47989U)

PCG_DEFINE_CONSTANT(uint32_t, default, multiplier, 747796405U)
PCG_DEFINE_CONSTANT(uint32_t, default, increment,  2891336453U)

PCG_DEFINE_CONSTANT(uint64_t, default, multiplier, 6364136223846793005ULL)
PCG_DEFINE_CONSTANT(uint64_t, default, increment,  1442695040888963407ULL)

PCG_DEFINE_CONSTANT(pcg128_t, default, multiplier,
        PCG_128BIT_CONSTANT(2549297995355413924ULL,4865540595714422341ULL))
PCG_DEFINE_CONSTANT(pcg128_t, default, increment,
        PCG_128BIT_CONSTANT(6364136223846793005ULL,1442695040888963407ULL))


/*
 * Each PCG generator is available in four variants, based on how it applies
 * the additive constant for its underlying LCG; the variations are:
 *
 *     single stream   - all instances use the same fixed constant, thus
 *                       the RNG always somewhere in same sequence
 *     mcg             - adds zero, resulting in a single stream and reduced
 *                       period
 *     specific stream - the constant can be changed at any time, selecting
 *                       a different random sequence
 *     unique stream   - the constant is based on the memory address of the
 *                       object, thus every RNG has its own unique sequence
 *
 * This variation is provided though mixin classes which define a function
 * value called increment() that returns the nesessary additive constant.
 */


/*
 * no stream (mcg)
 */

template <typename itype>
class no_stream {
protected:
    static constexpr bool is_mcg = true;

    // Is never called, but is provided for symmetry with specific_stream
    void set_stream(...)
    {
        abort();
    }

public:
    typedef itype state_type;

    static constexpr itype increment() {
        return 0;
    }

    static constexpr bool can_specify_stream = false;

    static constexpr size_t streams_pow2()
    {
        return 0u;
    }

protected:
    constexpr no_stream() = default;
};


/*
 * specific stream
 */

template <typename itype>
class specific_stream {
protected:
    static constexpr bool is_mcg = false;

    itype inc_ = default_increment<itype>::increment();

public:
    typedef itype state_type;
    typedef itype stream_state;

    constexpr itype increment() const {
        return inc_;
    }

    itype stream()
    {
         return inc_ >> 1;
    }

    void set_stream(itype specific_seq)
    {
         inc_ = (specific_seq << 1) | 1;
    }

    static constexpr bool can_specify_stream = true;

    static constexpr size_t streams_pow2()
    {
        return (sizeof(itype)*8) - 1u;
    }

protected:
    specific_stream() = default;

    specific_stream(itype specific_seq)
        : inc_(itype(specific_seq << 1) | itype(1U))
    {
        // Nothing (else) to do.
    }
};


/*
 * This is where it all comes together.  This function joins together three
 * mixin classes which define
 *    - the LCG additive constant (the stream)
 *    - the LCG multiplier
 *    - the output function
 * in addition, we specify the type of the LCG state, and the result type,
 * and whether to use the pre-advance version of the state for the output
 * (increasing instruction-level parallelism) or the post-advance version
 * (reducing register pressure).
 *
 * Given the high level of parameterization, the code has to use some
 * template-metaprogramming tricks to handle some of the suble variations
 * involved.
 */

template <typename xtype, typename itype,
          typename output_mixin,
          bool output_previous = true,
          typename stream_mixin = specific_stream<itype>,
          typename multiplier_mixin = default_multiplier<itype> >
class engine : protected output_mixin,
               public stream_mixin,
               protected multiplier_mixin {
protected:
    itype state_;

    struct can_specify_stream_tag {};
    struct no_specifiable_stream_tag {};

    using stream_mixin::increment;
    using multiplier_mixin::multiplier;

public:
    typedef xtype result_type;
    typedef itype state_type;

    static constexpr size_t period_pow2()
    {
        return sizeof(state_type)*8 - 2*stream_mixin::is_mcg;
    }

    // It would be nice to use std::numeric_limits for these, but
    // we can't be sure that it'd be defined for the 128-bit types.

    static constexpr result_type min()
    {
        return result_type(0UL);
    }

    static constexpr result_type max()
    {
        return result_type(~result_type(0UL));
    }

protected:
    itype bump(itype state)
    {
        return state * multiplier() + increment();
    }

    itype base_generate()
    {
        return state_ = bump(state_);
    }

    itype base_generate0()
    {
        itype old_state = state_;
        state_ = bump(state_);
        return old_state;
    }

public:
    result_type operator()()
    {
        if (output_previous)
            return this->output(base_generate0());
        else
            return this->output(base_generate());
    }

    result_type operator()(result_type upper_bound)
    {
        return bounded_rand(*this, upper_bound);
    }

protected:
    static itype advance(itype state, itype delta,
                         itype cur_mult, itype cur_plus);

    static itype distance(itype cur_state, itype newstate, itype cur_mult,
                          itype cur_plus, itype mask = ~itype(0U));

    itype distance(itype newstate, itype mask = itype(~itype(0U))) const
    {
        return distance(state_, newstate, multiplier(), increment(), mask);
    }

public:
    void advance(itype delta)
    {
        state_ = advance(state_, delta, this->multiplier(), this->increment());
    }

    void backstep(itype delta)
    {
        advance(-delta);
    }

    void discard(itype delta)
    {
        advance(delta);
    }

    bool wrapped()
    {
        if (stream_mixin::is_mcg) {
            // For MCGs, the low order two bits never change. In this
            // implementation, we keep them fixed at 3 to make this test
            // easier.
            return state_ == 3;
        } else {
            return state_ == 0;
        }
    }

    engine(itype state = itype(0xcafef00dd15ea5e5ULL))
        : state_(this->is_mcg ? state|state_type(3U)
                              : bump(state + this->increment()))
    {
        // Nothing else to do.
    }

    // This function may or may not exist.  It thus has to be a template
    // to use SFINAE; users don't have to worry about its template-ness.

    template <typename sm = stream_mixin>
    engine(itype state, typename sm::stream_state stream_seed)
        : stream_mixin(stream_seed),
          state_(this->is_mcg ? state|state_type(3U)
                              : bump(state + this->increment()))
    {
        // Nothing else to do.
    }

    template<typename SeedSeq>
    engine(SeedSeq&& seedSeq, typename std::enable_if<
                  !stream_mixin::can_specify_stream
               && !std::is_convertible<SeedSeq, itype>::value
               && !std::is_convertible<SeedSeq, engine>::value,
               no_specifiable_stream_tag>::type = {})
        : engine(generate_one<itype>(std::forward<SeedSeq>(seedSeq)))
    {
        // Nothing else to do.
    }

    template<typename SeedSeq>
    engine(SeedSeq&& seedSeq, typename std::enable_if<
                   stream_mixin::can_specify_stream
               && !std::is_convertible<SeedSeq, itype>::value
               && !std::is_convertible<SeedSeq, engine>::value,
        can_specify_stream_tag>::type = {})
        : engine(generate_one<itype,1,2>(seedSeq),
                 generate_one<itype,0,2>(seedSeq))
    {
        // Nothing else to do.
    }


    template<typename... Args>
    void seed(Args&&... args)
    {
        new (this) engine(std::forward<Args>(args)...);
    }

    template <typename xtype1, typename itype1,
              typename output_mixin1, bool output_previous1,
              typename stream_mixin_lhs, typename multiplier_mixin_lhs,
              typename stream_mixin_rhs, typename multiplier_mixin_rhs>
    friend bool operator==(const engine<xtype1,itype1,
                                     output_mixin1,output_previous1,
                                     stream_mixin_lhs, multiplier_mixin_lhs>&,
                           const engine<xtype1,itype1,
                                     output_mixin1,output_previous1,
                                     stream_mixin_rhs, multiplier_mixin_rhs>&);

    template <typename xtype1, typename itype1,
              typename output_mixin1, bool output_previous1,
              typename stream_mixin_lhs, typename multiplier_mixin_lhs,
              typename stream_mixin_rhs, typename multiplier_mixin_rhs>
    friend itype1 operator-(const engine<xtype1,itype1,
                                     output_mixin1,output_previous1,
                                     stream_mixin_lhs, multiplier_mixin_lhs>&,
                            const engine<xtype1,itype1,
                                     output_mixin1,output_previous1,
                                     stream_mixin_rhs, multiplier_mixin_rhs>&);

    template <typename CharT, typename Traits,
              typename xtype1, typename itype1,
              typename output_mixin1, bool output_previous1,
              typename stream_mixin1, typename multiplier_mixin1>
    friend std::basic_ostream<CharT,Traits>&
    operator<<(std::basic_ostream<CharT,Traits>& out,
               const engine<xtype1,itype1,
                              output_mixin1,output_previous1,
                              stream_mixin1, multiplier_mixin1>&);

    template <typename CharT, typename Traits,
              typename xtype1, typename itype1,
              typename output_mixin1, bool output_previous1,
              typename stream_mixin1, typename multiplier_mixin1>
    friend std::basic_istream<CharT,Traits>&
    operator>>(std::basic_istream<CharT,Traits>& in,
               engine<xtype1, itype1,
                        output_mixin1, output_previous1,
                        stream_mixin1, multiplier_mixin1>& rng);
};

template <typename CharT, typename Traits,
          typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin, typename multiplier_mixin>
std::basic_ostream<CharT,Traits>&
operator<<(std::basic_ostream<CharT,Traits>& out,
           const engine<xtype,itype,
                          output_mixin,output_previous,
                          stream_mixin, multiplier_mixin>& rng)
{
    using pcg_extras::operator<<;

    auto orig_flags = out.flags(std::ios_base::dec | std::ios_base::left);
    auto space = out.widen(' ');
    auto orig_fill = out.fill();

    pcg_extras::Output(out, rng.multiplier());
    out << space;
    pcg_extras::Output(out, rng.increment());
    out << space;
    pcg_extras::Output(out, rng.state_);
    out << space;

    out.flags(orig_flags);
    out.fill(orig_fill);
    return out;
}

template <typename CharT, typename Traits,
          typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin, typename multiplier_mixin>
std::basic_istream<CharT,Traits>&
operator>>(std::basic_istream<CharT,Traits>& in,
           engine<xtype,itype,
                    output_mixin,output_previous,
                    stream_mixin, multiplier_mixin>& rng)
{
    using pcg_extras::operator>>;

    auto orig_flags = in.flags(std::ios_base::dec | std::ios_base::skipws);

    itype multiplier, increment, state;
    pcg_extras::Input(in, multiplier);
    pcg_extras::Input(in, increment);
    pcg_extras::Input(in, state);
    if (!in.fail()) {
        bool good = true;
        if (multiplier != rng.multiplier()) {
           good = false;
        } else if (rng.can_specify_stream) {
           rng.set_stream(increment >> 1);
        } else if (increment != rng.increment()) {
           good = false;
        }
        if (good) {
            rng.state_ = state;
        } else {
            in.clear(std::ios::failbit);
        }
    }

    in.flags(orig_flags);
    return in;
}


template <typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin, typename multiplier_mixin>
itype engine<xtype,itype,output_mixin,output_previous,stream_mixin,
             multiplier_mixin>::advance(
    itype state, itype delta, itype cur_mult, itype cur_plus)
{
    // The method used here is based on Brown, "Random Number Generation
    // with Arbitrary Stride,", Transactions of the American Nuclear
    // Society (Nov. 1994).  The algorithm is very similar to fast
    // exponentiation.
    //
    // Even though delta is an unsigned integer, we can pass a
    // signed integer to go backwards, it just goes "the long way round".

    constexpr itype ZERO = 0u;  // itype may be a non-trivial types, so
    constexpr itype ONE  = 1u;  // we define some ugly constants.
    itype acc_mult = 1;
    itype acc_plus = 0;
    while (delta > ZERO) {
       if (delta & ONE) {
          acc_mult *= cur_mult;
          acc_plus = acc_plus*cur_mult + cur_plus;
       }
       cur_plus = (cur_mult+ONE)*cur_plus;
       cur_mult *= cur_mult;
       delta >>= 1;
    }
    return acc_mult * state + acc_plus;
}

template <typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin, typename multiplier_mixin>
itype engine<xtype,itype,output_mixin,output_previous,stream_mixin,
               multiplier_mixin>::distance(
    itype cur_state, itype newstate, itype cur_mult, itype cur_plus, itype mask)
{
    constexpr itype ONE  = 1u;  // itype could be weird, so use constant
    bool is_mcg = cur_plus == itype(0);
    itype the_bit = is_mcg ? itype(4u) : itype(1u);
    itype distance = 0u;
    while ((cur_state & mask) != (newstate & mask)) {
       if ((cur_state & the_bit) != (newstate & the_bit)) {
           cur_state = cur_state * cur_mult + cur_plus;
           distance |= the_bit;
       }
       assert((cur_state & the_bit) == (newstate & the_bit));
       the_bit <<= 1;
       cur_plus = (cur_mult+ONE)*cur_plus;
       cur_mult *= cur_mult;
    }
    return is_mcg ? distance >> 2 : distance;
}

template <typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin_lhs, typename multiplier_mixin_lhs,
          typename stream_mixin_rhs, typename multiplier_mixin_rhs>
itype operator-(const engine<xtype,itype,
                               output_mixin,output_previous,
                               stream_mixin_lhs, multiplier_mixin_lhs>& lhs,
               const engine<xtype,itype,
                               output_mixin,output_previous,
                               stream_mixin_rhs, multiplier_mixin_rhs>& rhs)
{
    static_assert(
        std::is_same<stream_mixin_lhs, stream_mixin_rhs>::value &&
            std::is_same<multiplier_mixin_lhs, multiplier_mixin_rhs>::value,
        "Incomparable generators");
    if (lhs.increment() == rhs.increment()) {
       return rhs.distance(lhs.state_);
    } else  {
       constexpr itype ONE = 1u;
       itype lhs_diff = lhs.increment() + (lhs.multiplier()-ONE) * lhs.state_;
       itype rhs_diff = rhs.increment() + (rhs.multiplier()-ONE) * rhs.state_;
       if ((lhs_diff & itype(3u)) != (rhs_diff & itype(3u))) {
           rhs_diff = -rhs_diff;
       }
       return rhs.distance(rhs_diff, lhs_diff, rhs.multiplier(), itype(0u));
    }
}


template <typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin_lhs, typename multiplier_mixin_lhs,
          typename stream_mixin_rhs, typename multiplier_mixin_rhs>
bool operator==(const engine<xtype,itype,
                               output_mixin,output_previous,
                               stream_mixin_lhs, multiplier_mixin_lhs>& lhs,
                const engine<xtype,itype,
                               output_mixin,output_previous,
                               stream_mixin_rhs, multiplier_mixin_rhs>& rhs)
{
    return    (lhs.multiplier() == rhs.multiplier())
           && (lhs.increment()  == rhs.increment())
           && (lhs.state_       == rhs.state_);
}

template <typename xtype, typename itype,
          typename output_mixin, bool output_previous,
          typename stream_mixin_lhs, typename multiplier_mixin_lhs,
          typename stream_mixin_rhs, typename multiplier_mixin_rhs>
inline bool operator!=(const engine<xtype,itype,
                               output_mixin,output_previous,
                               stream_mixin_lhs, multiplier_mixin_lhs>& lhs,
                       const engine<xtype,itype,
                               output_mixin,output_previous,
                               stream_mixin_rhs, multiplier_mixin_rhs>& rhs)
{
    return !operator==(lhs,rhs);
}


template <typename xtype, typename itype,
         template<typename XT,typename IT> class output_mixin,
         bool output_previous = (sizeof(itype) <= 8),
         template<typename IT> class multiplier_mixin = default_multiplier>
using setseq_base = engine<xtype, itype,
                         output_mixin<xtype, itype>, output_previous,
                         specific_stream<itype>,
                         multiplier_mixin<itype> >;

template <typename xtype, typename itype,
         template<typename XT,typename IT> class output_mixin,
         bool output_previous = (sizeof(itype) <= 8),
         template<typename IT> class multiplier_mixin = default_multiplier>
using mcg_base = engine<xtype, itype,
                      output_mixin<xtype, itype>, output_previous,
                      no_stream<itype>,
                      multiplier_mixin<itype> >;

/*
 * OUTPUT FUNCTIONS.
 *
 * These are the core of the PCG generation scheme.  They specify how to
 * turn the base LCG's internal state into the output value of the final
 * generator.
 *
 * They're implemented as mixin classes.
 *
 * All of the classes have code that is written to allow it to be applied
 * at *arbitrary* bit sizes, although in practice they'll only be used at
 * standard sizes supported by C++.
 */

/*
 * XSH RR -- high xorshift, followed by a random rotate
 *
 * Fast.  A good performer.  Slightly better statistically than XSH RS.
 */

template <typename xtype, typename itype>
struct xsh_rr_mixin {
    static xtype output(itype internal)
    {
        constexpr bitcount_t bits        = bitcount_t(sizeof(itype) * 8);
        constexpr bitcount_t xtypebits   = bitcount_t(sizeof(xtype)*8);
        constexpr bitcount_t sparebits   = bits - xtypebits;
        constexpr bitcount_t wantedopbits =
                              xtypebits >= 128 ? 7
                            : xtypebits >=  64 ? 6
                            : xtypebits >=  32 ? 5
                            : xtypebits >=  16 ? 4
                            :                    3;
        constexpr bitcount_t opbits =
                              sparebits >= wantedopbits ? wantedopbits
                                                        : sparebits;
        constexpr bitcount_t amplifier = wantedopbits - opbits;
        constexpr bitcount_t mask = (1 << opbits) - 1;
        constexpr bitcount_t topspare    = opbits;
        constexpr bitcount_t bottomspare = sparebits - topspare;
        constexpr bitcount_t xshift      = (topspare + xtypebits)/2;
        bitcount_t rot = opbits ? bitcount_t(internal >> (bits - opbits)) & mask
                                : 0;
        bitcount_t amprot = (rot << amplifier) & mask;
        internal ^= internal >> xshift;
        xtype result = xtype(internal >> bottomspare);
        result = rotr(result, amprot);
        return result;
    }
};


/*
 * XSL RR -- fixed xorshift (to low bits), random rotate
 *
 * Useful for 128-bit types that are split across two CPU registers.
 */

template <typename xtype, typename itype>
struct xsl_rr_mixin {
    static xtype output(itype internal)
    {
        constexpr bitcount_t xtypebits = bitcount_t(sizeof(xtype) * 8);
        constexpr bitcount_t bits = bitcount_t(sizeof(itype) * 8);
        constexpr bitcount_t sparebits = bits - xtypebits;
        constexpr bitcount_t wantedopbits = xtypebits >= 128 ? 7
                                       : xtypebits >=  64 ? 6
                                       : xtypebits >=  32 ? 5
                                       : xtypebits >=  16 ? 4
                                       :                    3;
        constexpr bitcount_t opbits = sparebits >= wantedopbits ? wantedopbits
                                                             : sparebits;
        constexpr bitcount_t amplifier = wantedopbits - opbits;
        constexpr bitcount_t mask = (1 << opbits) - 1;
        constexpr bitcount_t topspare = sparebits;
        constexpr bitcount_t bottomspare = sparebits - topspare;
        constexpr bitcount_t xshift = (topspare + xtypebits) / 2;

        bitcount_t rot =
            opbits ? bitcount_t(internal >> (bits - opbits)) & mask : 0;
        bitcount_t amprot = (rot << amplifier) & mask;
        internal ^= internal >> xshift;
        xtype result = xtype(internal >> bottomspare);
        result = rotr(result, amprot);
        return result;
    }
};



/* ---- End of Output Functions ---- */

} // namespace pcg_detail

namespace pcg_engines {

using namespace pcg_detail;

/* Predefined types for XSH RR */

typedef setseq_base<uint32_t, uint64_t, xsh_rr_mixin>  setseq_xsh_rr_64_32;


/* Predefined types for XSL RR (only defined for "large" types) */

typedef setseq_base<uint64_t, pcg128_t, xsl_rr_mixin>  setseq_xsl_rr_128_64;
typedef mcg_base<uint64_t, pcg128_t, xsl_rr_mixin>  mcg_xsl_rr_128_64;

} // namespace pcg_engines

typedef pcg_engines::setseq_xsh_rr_64_32        pcg32;

typedef pcg_engines::setseq_xsl_rr_128_64       pcg64;
typedef pcg_engines::mcg_xsl_rr_128_64          pcg64_fast;


#ifdef _MSC_VER
    #pragma warning(default:4146)
#endif

#endif // PCG_RAND_HPP_INCLUDED
