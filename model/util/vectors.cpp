/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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

#include "util/vectors.h"
#include <cstring>


void vectors::scale (vector<double>& vec, double a) {
  for (size_t i = 0; i < vec.size(); ++i)
    vec[i] *= a;
}

double vectors::sum (vector<double>& vec) {
  double r = 0.0;
  for (size_t i = 0; i < vec.size(); ++i)
    r += vec[i];
  return r;
}


bool vectors::approxEqual (const double a, const double b) {
  const double LIM = 1e-6;
  return (fabs(a-b) <= max(fabs(a),fabs(b)) * LIM);
}

bool vectors::approxEqual (const vector<double>& vec1, const vector<double>& vec2) {
  if (vec1.size() != vec2.size())
    return false;
  for (size_t i = 0; i < vec1.size(); ++i) {
    if (!approxEqual (vec1[i], vec2[i]))
      return false;
  }
  return true;
}


vector<double> vectors::gsl2std (const gsl_vector* vec) {
  vector<double> ret (vec->size);
  memcpy (&ret[0], vec->data, ret.size()*sizeof(double));
  return ret;
}


gsl_vector* vectors::std2gsl (const vector<double>& vec, size_t length) {
  if (vec.size() != length)
    throw length_error ("vectorStd2Gsl: vec has incorrect length");
  gsl_vector* ret = gsl_vector_alloc (length);
  memcpy (ret->data, &vec[0], length * sizeof(double));
  return ret;
}

gsl_vector* vectors::std2gsl (const double* vec, size_t length) {
  gsl_vector* ret = gsl_vector_alloc (length);
  memcpy (ret->data, vec, length * sizeof(double));
  return ret;
}


vector<double> vectors::DoubleList2std (const scnXml::DoubleList& list, size_t length) {
  const scnXml::DoubleList::ItemSequence seq = list.getItem();
  if (seq.size() != length)
    throw xml_scenario_error ("readDoubleList: XML list has wrong length");
  vector<double> ret (length);
  for (size_t i = 0; i < length; ++i)
    ret[i] = seq[i];
  return ret;
}
