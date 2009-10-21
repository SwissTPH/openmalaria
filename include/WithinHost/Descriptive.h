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

#include "WithinHost/WithinHostModel.h"
#include "WithinHost/DescriptiveInfection.h"

using namespace std;

/*! Old Within Host Model class.
 */
class DescriptiveWithinHostModel : public WithinHostModel {
public:
  /// Create a new WHM
  DescriptiveWithinHostModel();
  /// Load a DescriptiveWithinHostModel, including infections.
  DescriptiveWithinHostModel(istream& in);
  virtual ~DescriptiveWithinHostModel();
  
  /** @brief Checkpointing of variables */
  //@{
  virtual void write(ostream& out) const;
protected:
  /** Special checkpointing constructor for derived use.
   *
   * Same as the other checkpointing constructor except that
   * this one doesn't load infections. */
  DescriptiveWithinHostModel(istream& in, bool derived);
  /// Called by both checkpointing constructors
  void readDescriptiveWHM(istream& in);
  /// Called by write() and derived write() functions.
  void writeDescriptiveWHM(ostream& out) const;
  //@}
  
public:
  virtual void update();
  
  virtual void summarize(double age);
  
  //! Create a new infection requires that the human is allocated and current
  virtual void newInfection();
  
  //! Clears all infections in an individual
  virtual void clearAllInfections();
  
  void calculateDensities(double ageInYears, double BSVEfficacy);
  
  /*! Until now, this only includes decay of immunity against
  asexual blood stages */
  virtual void updateImmuneStatus();
  
  virtual void immunityPenalisation();
  
  bool parasiteDensityDetectible() const {
    return totalDensity > detectionLimit;
  }
  
protected:
  /*!  SP drug action applies to each infection depending on genotype and when
  the individual had their last dose of SP */
  virtual void SPAction();
  
  virtual void IPTattenuateAsexualDensity (DescriptiveInfection& infec);
  virtual void IPTattenuateAsexualMinTotalDensity ();
  
  static const int MAX_INFECTIONS;
  
  //!multiplicity of infection
  int _MOI;
  
  /** The list of all infections this human has.
   * 
   * Since infection models and within host models are very much intertwined,
   * the idea is that each WithinHostModel has its own list of infections. */
  std::list<DescriptiveInfection*> infections;
  
  //!Cumulative parasite density since birth
  double _cumulativeY;
  
private:
  //!Number of infections received since birth
  double _cumulativeh;
  //!cumulativeY from previous timestep
  double _cumulativeYlag;
  
  //!innate ability to control parasite densities
  double _innateImmunity;
  
  //!Number of infections with densities above the limit of detection
  int patentInfections;
};

#endif
