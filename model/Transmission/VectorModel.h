/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_VectorModel
#define Hmod_VectorModel

#include "Global.h"
#include "Transmission/TransmissionModel.h"
#include "Transmission/Anopheles/AnophelesModel.h"

namespace scnXml {
  class Vector;
}

namespace OM {
    class Population;
namespace Transmission {
    using Anopheles::AnophelesModel;
    
/** Transmission models, Chitnis et al.
 * 
 * This class contains code for species-independent components. Per-species
 * code is in the Anopheles directory and namespace. */
class VectorModel : public TransmissionModel {
public:
  VectorModel(const scnXml::Entomology& entoData, const scnXml::Vector vectorData, int populationSize);
  virtual ~VectorModel();
  
  /** Extra initialisation when not loading from a checkpoint, requiring
   * information from the human population structure. */
  virtual void init2 ();
  
  virtual void initVectorInterv( const scnXml::Description::AnophelesSequence& list,
        size_t instance, const string& name );
  virtual void initVectorTrap( const scnXml::VectorTrap::DescriptionSequence list,
        size_t instance, const scnXml::VectorTrap::NameOptional name );
  
  virtual void scaleEIR (double factor);
//   virtual void scaleXML_EIR (scnXml::Entomology&, double factor) const;
  
  virtual SimTime minPreinitDuration ();
  virtual SimTime expectedInitDuration ();
  virtual SimTime initIterate ();
  
  virtual void vectorUpdate ();
  virtual void update ();

  virtual double calculateEIR( Host::Human& human, double ageYears,
        vector<double>& EIR );
  
  virtual const map<string,size_t>& getSpeciesIndexMap();
  virtual void deployVectorPopInterv (size_t instance);
  virtual void deployVectorTrap( size_t instance, double number, SimTime lifespan );
  virtual void uninfectVectors();
  
  virtual void summarize ();
  
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
  void ctsCbN_v0 (ostream& stream);
  void ctsCbP_A (ostream& stream);
  void ctsCbP_df (ostream& stream);
  void ctsCbP_dif (ostream& stream);
  void ctsCbN_v (ostream& stream);
  void ctsCbO_v (ostream& stream);
  void ctsCbS_v (ostream& stream);
  void ctsCbAlpha (const Population& population, ostream& stream);
  void ctsCbP_B (const Population& population, ostream& stream);
  void ctsCbP_CD (const Population& population, ostream& stream);
  void ctsNetInsecticideContent (const Population& population, ostream& stream);
  void ctsIRSInsecticideContent (const Population& population, ostream& stream);
  void ctsIRSEffects (const Population& population, ostream& stream);
  void ctsCbResAvailability (ostream& stream);
  void ctsCbResRequirements (ostream& stream);
  
  
    /// Number of iterations performed during initialization.
    /// 
    /// Special cases: 0 during initial one-human-lifespan warmup, <0 after
    /// initialisation.
    int initIterations;
    
  /** @brief Access to per (anopheles) species data.
   *
   * Set by constructor so don't checkpoint. */
  //@{
  /** The number of discrete species of anopheles mosquitos to be modelled.
   *
   * Must be the same as species.size() and at least 1. */
  size_t numSpecies;
  
  /** Per anopheles species data.
   *
   * Array will be recreated by constructor, but some members of AnophelesModel
   * need to be checkpointed. */
  vector<AnophelesModel> species;
  
  /** A map of anopheles species/variant name to an index in species.
   *
   * When the main ento data is read from XML, each anopheles section is read
   * into one index of the species array. The name and this index are added
   * here.
   * 
   * Other data read from XML should look up the name here and use the index
   * found. Doesn't need checkpointing. */
  map<string,size_t> speciesIndex;
  //@}
  
  /// @brief Saved data for use in initialisation / fitting cycle
  //@{
    util::vecDay2D<double> saved_sum_avail;
    util::vecDay2D<double> saved_sigma_df;
    util::vecDay3D<double> saved_sigma_dif;
  //@}
  
  friend class PerHost;
  friend class AnophelesModelSuite;
};

} }
#endif
