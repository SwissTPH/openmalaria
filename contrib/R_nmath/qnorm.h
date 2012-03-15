/*
 *  Header : A C Library of Special Functions
 *  Copyright (C) 1998       Ross Ihaka
 *  Copyright (C) 2000--2005 The R Development Core Team
 *  based on AS 111 (C) 1977 Royal Statistical Society
 *  and   on AS 241 (C) 1988 Royal Statistical Society
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#ifndef R_NMATH_QNORM
#define R_NMATH_QNORM

namespace R {
/**
 *  SYNOPSIS
 *
 *      double qnorm5(double p, double mu, double sigma,
 *                    int lower_tail, int log_p)
 *            {qnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *      Compute the quantile function for the normal distribution.
 *
 *      For small to moderate probabilities, algorithm referenced
 *      below is used to obtain an initial approximation which is
 *      polished with a final Newton step.
 *
 *      For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *      Beasley, J. D. and S. G. Springer (1977).
 *      Algorithm AS 111: The percentage points of the normal distribution,
 *      Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */
double qnorm5(double p, double mu, double sigma, bool lower_tail, bool log_p);
}

#endif
