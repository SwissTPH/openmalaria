/*

  This file is part of OpenMalaria.
 
  Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
 
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

#ifndef Hmod_oldwhost
#define Hmod_oldwhost

#include <iostream>
#include "global.h"
#include "withinHostModel.h"

using namespace std;

class Human;
class Infection;

/*! Old Within Host Model class.
 */
class OldWithinHostModel : public WithinHostModel {
public:
  OldWithinHostModel();
  ~OldWithinHostModel();
  
  friend ostream& operator<<(ostream& out, const WithinHostModel &model);
  friend istream& operator>>(istream& in, WithinHostModel &model);


  virtual void update(double age);
  
  virtual void summarize(double age);
  
  //! Create a new infection requires that the human is allocated and current
  virtual void newInfection(int);

  /*!  Clears all infections which have expired (their startdate+duration is less
  than the current time). */
  virtual void clearOldInfections();

  //! Clears all infections in an individual
  virtual void clearAllInfections();
  
  void medicate(string drugName, double qty, int time);

  void calculateDensities(Human&);
  
  void write(ostream& out) const;
  void read(istream& in);

private:
  /*!  SP drug action applies to each infection depending on genotype and when
  the individual had their last dose of SP */
  void SPAction(Human&);

  void calculateDensity(Infection *inf, double);
  
  void treatInfections();
  //! Treats all infections in an individual
  void treatAllInfections();
  
  //! time at which attenuated infection 'would' end if SP present
  int _SPattenuationt;
  double cumulativeY;
  double cumulativeh;
  double timeStepMaxDensity;
  
  //!multiplicity of infection
  int _MOI;
  //!Number of infections with densities above the limit of detection
  int patentInfections;
  
  /// Encapsulates drug code for each human
  DrugProxy _proxy;
  
  /** The list of all infections this human has.
   * 
   * Since infection models and within host models are very much intertwined,
   * the idea is that each WithinHostModel has its own list of infections. */
  std::list<Infection*> infections;
};

#endif
