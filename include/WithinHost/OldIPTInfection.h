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

#ifndef Hmod_OldIPTInfection
#define Hmod_OldIPTInfection
#include "DescriptiveInfection.h"


struct genotype {
  int ID;
  //!In order to save memory, we just define the ID of the genotype. Attributes of the
  //!genotype can be accessed via arrays in mod_intervention.
  //!(e.g. freq = mod_intervention.GenotypeFreq(iTemp%iData%gType%ID)
  //!attributes are:
  //!freq: Probability of being infected by this specific genotype
  //!ACR: Probability of being cured (due to SP)
  //!proph: Prophylactic effect of SP (measured in time steps)
  //!tolperiod: time window of tolerance period
  //!SPattenuation: Factor of how parasites are attenuated  by SP (genotype specific)
};

namespace scnXml {
  class Interventions;	// XML data passed to initParameters
}

/// IPT extension of DescriptiveInfection
class OldIPTInfection : public DescriptiveInfection {
public:
  ///@name Static init/cleanup
  //@{
  static void initParameters(const scnXml::Interventions& xmlInterventions);
  static void clearParameters();
  //@}
  
  ///@name CTOR & DTOR
  //@{
  //! Constructor
  /*! \param lastSPdose Time interval of last SP Dose. */
  OldIPTInfection(int lastSPdose, int simulationTime);
  
  /** Checkpoint-reading constructor */
  OldIPTInfection (istream& in);
  
  /** Destructor */
  virtual ~OldIPTInfection() {}
  //@}
  
  virtual void write (ostream& out) const;
  
  /** @name Getter functions
   *
   * Avoid making these virtual by casting DescriptiveInfection references to
   * OldIPTInfection type. */
  int getGenoTypeID() { return _gType.ID; };
  bool getSPattenuate() { return _SPattenuate; };
  
private:
  //! Genotype responsible for infection
  genotype _gType;
  //! IPTi parameter (indicator for attenuation).
  bool _SPattenuate;
  
  /// @name genotypes
  //@{
  static int numberOfGenoTypes;
  static double *genotypeFreq;
  static int *genotypeTolPeriod;
public:
  static int *genotypeProph;
  static double *genotypeACR;
  static double *genotypeAtten;
  //@}
};

#endif
