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

#ifndef Hmod_WithinHost_Descriptive
#define Hmod_WithinHost_Descriptive

#include "WithinHost/WithinHostModel.h"
#include "WithinHost/DescriptiveInfection.h"

using namespace std;

namespace OM { namespace WithinHost {

/*! Old Within Host Model class.
 *
 * Note: this implementation has a few bugs with (hopefully) small effect
 * conditionally fixed (see MAX_DENS_CORRECTION and
 * INNATE_MAX_DENS). Thus results can be preserved. */
class DescriptiveWithinHostModel : public WithinHostModel {
public:
  /// Create a new WHM
  DescriptiveWithinHostModel();
  virtual ~DescriptiveWithinHostModel();
  
  virtual void newInfection();
  /// load an infection from a checkpoint
  virtual void loadInfection(istream& stream);
  virtual void clearAllInfections();
  
  virtual void calculateDensities(double ageInYears, double BSVEfficacy);
  
protected:
  virtual int countInfections (int& patentInfections);
  
  ///@brief IPT extensions − empty otherwise
  //@{
  /*!  SP drug action applies to each infection depending on genotype and when
  the individual had their last dose of SP */
  virtual void SPAction() {}
  virtual void IPTattenuateAsexualMinTotalDensity () {}
  virtual void IPTattenuateAsexualDensity (DescriptiveInfection* inf) {}
  //@}
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
  /** The list of all infections this human has.
   * 
   * Since infection models and within host models are very much intertwined,
   * the idea is that each WithinHostModel has its own list of infections. */
  std::list<DescriptiveInfection*> infections;
};

} }
#endif