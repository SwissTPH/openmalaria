/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2005-2011 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

// NOTE: I cheat and include definitions as well as declarations in this header.
// For now this is fine since the header is only included once, but if that
// were changed it might be necessary to move definitions to a cpp file.
#ifndef Hmod_GSL_Multidim_solver
#define Hmod_GSL_Multidim_solver

#include "gsl_multiroots.h"
#include "gsl_multimin.h"

namespace OM {
namespace util {

/** Interface for multidimensional minimisation/root-finding algorithms. */
struct MultidimSolver {
    /** Iterate.
     * 
     * @returns GSL_SUCCESS on success, GSL_ENOPROG if no progress is being
     *  made, or another error code.
     */
    virtual int iterate() =0;
    /** Return true if a minimum/root has been found to a sufficient degree of
     * accuracy: e_abs. */
    virtual bool success(double e_abs) =0;
    /** Get the current best estimate of function input. */
    virtual gsl_vector* x() =0;
};

/** Implementation using GSL multidimensional minimisation algorithms. */
struct MultidimMinimiser : public MultidimSolver {
    /** Initialise and set variables.
     * 
     * @param algorithm Algorithm: gsl_multimin_fminimizer_nmsimplex2,
     *  gsl_multimin_fminimizer_nmsimplex or gsl_multimin_fminimizer_nmsimplex2rand
     * @param n Number of dimensions of function input
     * @param function The function to sample from
     * @param params Extra state data passed to the function
     * @param x Initial guess
     * @param step_size Initial step size
     */
    MultidimMinimiser( const gsl_multimin_fminimizer_type *algorithm,
                       size_t n,
                       double (*function) (const gsl_vector *x, void *params),
                       void *params,
                       gsl_vector *x,
                       gsl_vector *step_size
    ){
        gsl_multimin_function func = { function, n, params };
        s = gsl_multimin_fminimizer_alloc( algorithm, n );
        gsl_multimin_fminimizer_set( s, &func, x, step_size );
    }
    ~MultidimMinimiser() {
        gsl_multimin_fminimizer_free (s);
    }
    virtual int iterate() {
        return gsl_multimin_fminimizer_iterate( s );
    }
    virtual bool success(double e_abs) {
        double size = gsl_multimin_fminimizer_size (s);
        return gsl_multimin_test_size (size, e_abs) == GSL_SUCCESS;
    }
    virtual gsl_vector* x() {
        return gsl_multimin_fminimizer_x( s );
    }
private:
    gsl_multimin_fminimizer *s;
};

/** Implementation using GSL multidimensional root-finding algorithms. */
struct MultidimRootFinder : public MultidimSolver {
    /** Initialise and set variables.
     * 
     * @param algorithm Algorithm: gsl_multiroot_fsolver_type, gsl_multiroot_fsolver_hybrid,
     *  gsl_multiroot_fsolver_dnewton or gsl_multiroot_fsolver_broyden
     * @param n Number of dimensions of function input
     * @param function The function to sample from
     * @param params Extra state data passed to the function
     * @param x Initial guess
     */
    MultidimRootFinder( const gsl_multiroot_fsolver_type *algorithm,
                       size_t n,
                       int (*function) (const gsl_vector *x, void *params, gsl_vector *f),
                       void *params,
                       gsl_vector *x
    ){
        gsl_multiroot_function func = { function, n, params };
        s = gsl_multiroot_fsolver_alloc( algorithm, n );
        gsl_multiroot_fsolver_set( s, &func, x );
    }
    ~MultidimRootFinder() {
        gsl_multiroot_fsolver_free (s);
    }
    virtual int iterate() {
        return gsl_multiroot_fsolver_iterate( s );
    }
    virtual bool success( double e_abs) {
        return gsl_multiroot_test_residual (s->f, e_abs) == GSL_SUCCESS;
    }
    virtual gsl_vector* x() {
        return gsl_multiroot_fsolver_root( s );
    }
private:
    gsl_multiroot_fsolver *s;
};

}
}

#endif
