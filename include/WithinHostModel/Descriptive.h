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
#include "WithinHostModel.h"
#include "WithinHostModel/DescriptiveInfection.h"
#include "drug.h"

using namespace std;

class Human;

/*! Old Within Host Model class.
 */
class DescriptiveWithinHostModel : public WithinHostModel {
public:
  DescriptiveWithinHostModel();
  DescriptiveWithinHostModel(istream& in);
  virtual ~DescriptiveWithinHostModel();
  
  virtual void update(double age);
  
  virtual void summarize(double age);
  
  //! Create a new infection requires that the human is allocated and current
  virtual void newInfection();

  /*!  Clears all infections which have expired (their startdate+duration is less
  than the current time). */
  virtual void clearOldInfections();

  //! Clears all infections in an individual
  virtual void clearAllInfections();
  
  void medicate(string drugName, double qty, int time);

  void calculateDensities(Human&);
  
  /*! Until now, this only includes decay of immunity against
  asexual blood stages */
  virtual void updateImmuneStatus();
  
  virtual void immunityPenalisation();
  
  virtual void write(ostream& out) const;
  
protected:
  /*!  SP drug action applies to each infection depending on genotype and when
  the individual had their last dose of SP */
  virtual void SPAction(Human&);
  
  virtual void IPTattenuateAsexualDensity (DescriptiveInfection& infec);
  virtual void IPTattenuateAsexualMinTotalDensity (Human&);
  
  /** @brief Checkpointing of variables in DescriptiveWithinHostModel.
   *
   * Since write() is replaced in sub-classes, having separately
   * named methods allows the same code to be used in both classes. */
  //@{
  void writeDescriptiveWHM(ostream& out) const;
  //@}
  
  //!multiplicity of infection
  int _MOI;
  
  /** The list of all infections this human has.
   * 
   * Since infection models and within host models are very much intertwined,
   * the idea is that each WithinHostModel has its own list of infections. */
  std::list<DescriptiveInfection*> infections;
  
  /** Just push back a new infection of the relevant type, created from in. */
  virtual void loadInfection (istream& in);

  //!Cumulative parasite density since birth
  double _cumulativeY;
  
  /** Used within calculateDensities and other functions, set each call.
   *
   * Doesn't need to be checkpointed. */
  double timeStepMaxDensity;
  
private:
  //TODO: check why we have 2 cumulativeh and cumulativeY params
  //!Number of infections received since birth
  double _cumulativeh;
  //!cumulativeY from previous timestep
  double _cumulativeYlag;
  
  //!innate ability to control parasite densities
  double _innateImmunity;
  
  //!Number of infections with densities above the limit of detection
  int patentInfections;
  
  /// Encapsulates drug code for each human
  DrugProxy _proxy;
};

#endif
