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

#ifndef Hmod_DescriptiveIPTInfection
#define Hmod_DescriptiveIPTInfection
#include "DescriptiveInfection.h"

namespace scnXml {
    class Interventions;	// XML data passed to initParameters
}

namespace OM { namespace WithinHost {
    

/** IPT extension of DescriptiveInfection
 *
 * Note: proteome_ID parameter from base Infection is used here to store the genotype.
 *
 * NOTE: This IPT code (this class and DescriptiveIPTWithinHost) are
 * unmaintained in order to keep results comparable with previous experiments
 * run. */
class DescriptiveIPTInfection : public DescriptiveInfection {
public:
  ///@name Static init/cleanup
  //@{
  static void initParameters(const scnXml::Interventions& xmlInterventions);
  //@}
  
  ///@name CTOR & DTOR
  //@{
  //! Constructor
  /*! \param lastSPdose Time interval of last SP Dose. */
  DescriptiveIPTInfection(TimeStep lastSPdose);
  DescriptiveIPTInfection (istream& stream);
  
  /** Destructor */
  virtual ~DescriptiveIPTInfection() {}
  //@}
  
  /** The event that the last SP dose clears parasites. */
  bool eventSPClears (TimeStep _lastSPDose);
  
  /// Return: _SPattenuate == 1. Name by DH.
  bool doSPAttenuation () { return _SPattenuate == 1; }
  double asexualAttenuation ();
  /// Extraction by DH; probably not most accurate name.
  TimeStep getAsexualAttenuationEndDate () {
    return _startdate + TimeStep(static_cast<int>(_duration.asInt() * genotypes[proteome_ID].atten));	//FIXME: should probably add latentp
  }
  
protected:
    virtual void checkpoint (ostream& stream);
    
private:
  //! IPTi parameter (indicator for attenuation).
  bool _SPattenuate;
  
  //!In order to save memory, we just define the ID of the genotype. Attributes of the
  //!genotype can be accessed via arrays in mod_intervention.
  //!(e.g. freq = mod_intervention.GenotypeFreq(iTemp%iData%gType%ID)
  struct GenotypeData {
      GenotypeData(double cF, TimeStep tP, TimeStep p, double acr, double at) :
        cumFreq(cF), tolPeriod(tP), proph(p), ACR(acr), atten(at) {}
    //!freq: Probability of being infected by this specific genotype
    double cumFreq;
    //!tolperiod: time window of tolerance period
    TimeStep tolPeriod;
    //!proph: Prophylactic effect of SP (measured in time steps)
    TimeStep proph;
    //!ACR: Probability of being cured (due to SP)
    double ACR;
    //!SPattenuation: Factor of how parasites are attenuated  by SP (genotype specific)
    double atten;
  };
  /// Per genotype data, set by initParameters
  static vector<GenotypeData> genotypes;
};

} }
#endif