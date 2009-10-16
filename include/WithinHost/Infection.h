/*
 This file is part of OpenMalaria.
 
 Copyright (C) 2005-2009 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
 OpenMalaria is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or (at
 your option) any later version.
 
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef Hmod_Infection
#define Hmod_Infection

#include "Global.h"
#include "proteome.h"


class Infection {
public:
  static float cumulativeYstar; //!< Critical value for immunity trigger (cumulative densities)
  static float cumulativeHstar; //!< Critical value for immunity trigger (cumulative inoculations)
  
  Infection (int now) :
    _startdate(now),
    _density(0.0)
  {}
  Infection (istream& in);
  
  void writeInfection (ostream& out) const;
  
  //! Get proteome
  inline const ProteomeInstance* getProteome() const {
    return _proteome;
  }
  
protected:
  //! Proteome (used in a different situation than genotype) 
  const ProteomeInstance* _proteome; 
  
  //! Start date of the infection
  int _startdate; 
  //! Current density of the infection
  double _density;
};

#endif
