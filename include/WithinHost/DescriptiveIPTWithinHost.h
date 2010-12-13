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

#ifndef Hmod_DescriptiveIPT
#define Hmod_DescriptiveIPT

#include "WithinHost/DescriptiveWithinHost.h"

namespace OM { namespace WithinHost {
    
/** Extension to the DescriptiveWithinHostModel, including IPT (intermittent
 * preventative treatment) using a simple drug-action model (SPAction).
 *
 * NOTE: This IPT code (this class and DescriptiveIPTInfection) are
 * unmaintained in order to keep results comparable with previous experiments
 * run. */
class DescriptiveIPTWithinHost : public DescriptiveWithinHostModel {
public:
  ///@name Static init/cleanup
  //@{
  /** Determines whether IPT is present or not (iptActive), and if so
   * initialises parameters here and in DescriptiveIPTInfection. */
  static void init();
  static void cleanup();
  //@}
  
  DescriptiveIPTWithinHost ();
  
  //! Create a new infection requires that the human is allocated and current
  virtual void newInfection();
  virtual void loadInfection(istream& stream);
  
  /// Conditionally clear all infections
  virtual void clearInfections (bool isSevere);
  /// Continuous intervention: give an IPTi dose
  virtual void deployIptDose (Monitoring::AgeGroup ageGroup);
  /// Prescribe IPTi with probability compliance. Only called if IPT present.
  virtual void IPTiTreatment (Monitoring::AgeGroup ageGroup);
  
  /// Is IPT present?
  /// set by initParameters
  static bool iptActive;
  
protected:
  virtual bool eventSPClears (DescriptiveInfection* inf);
  virtual void IPTattenuateAsexualMinTotalDensity ();
  virtual void IPTattenuateAsexualDensity (DescriptiveInfection* inf);
  
  virtual void checkpoint (istream& stream);
  virtual void checkpoint (ostream& stream);
  
private:
  //! time at which attenuated infection 'would' end if SP present
  int _SPattenuationt;
  /** Timestep of last SP Dose given (TIMESTEP_NEVER if no SP dose given). */
  int _lastSPDose;
  /// Timestep of last IPTi or placebo dose given (TIMESTEP_NEVER if never given).
  int _lastIptiOrPlacebo;
  
  //!Cumulative number of infections since birth
  //NOTE: we also have _cumulativeh; why both?
  int _cumulativeInfections;
  
  /// @brief Static data set by initParameters
  //@{
  /// Values (codes)
  static int iptiEffect;
  //@}
};

} }
#endif