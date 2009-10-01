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

#ifndef Hmod_dummywithinhost
#define Hmod_dummywithinhost

#include "global.h"
#include "WithinHost/WithinHostModel.h"
#include "WithinHost/DummyInfection.h"
#include "Drug/DrugModel.h"
;

using namespace std;

/*! Dummy Within Host Model class.
 */
class DummyWithinHostModel : public WithinHostModel {
public:
  DummyWithinHostModel();
  DummyWithinHostModel(istream& in);
  ~DummyWithinHostModel();
  

  virtual void update();
  
  virtual void summarize(double age);
  
  //! Create a new infection requires that the human is allocated and current
  virtual void newInfection();

  /*!  Clears all infections which have expired (their startdate+duration is less
  than the current time). */
  virtual void clearOldInfections();

  //! Clears all infections in an individual
  virtual void clearAllInfections();
  
  void medicate(string drugName, double qty, int time, double age);

  void calculateDensities(double ageInYears, double BSVEfficacy);
  
  /*! Until now, this only includes decay of immunity against
  asexual blood stages */
  virtual void updateImmuneStatus();
  
  virtual void immunityPenalisation();
  
  void write(ostream& out) const;
  
  bool parasiteDensityDetectible() const {
    return totalDensity > detectionLimit;
  }
  
private:
  /// Encapsulates drug code for each human
  DrugModel* drugProxy;
  
  //TODO: check why we have 2 cumulativeh and cumulativeY params
  //!Number of infections received since birth
  double _cumulativeh;
  //!Cumulative parasite density since birth
  double _cumulativeY;
  //!cumulativeY from previous timestep
  double _cumulativeYlag;
  
  double cumulativeY;
  double cumulativeh;
  
  //!multiplicity of infection
  int _MOI;
  //!Number of infections with densities above the limit of detection
  int patentInfections;
  
  /** The list of all infections this human has.
   * 
   * Since infection models and within host models are very much intertwined,
   * the idea is that each WithinHostModel has its own list of infections. */
  std::list<DummyInfection> infections;
  
  static const int MAX_INFECTIONS;
};

#endif
