/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005,2006,2007,2008 Swiss Tropical Institute and Liverpool School Of Tropical Medicine
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
#include "Transmission/ITN.h"
#include "Transmission/IRS.h"

namespace scnXml {
  class Vector;
}

namespace OM { namespace Transmission {
    
/** Transmission models, Chitnis et al.
 * 
 * This class contains code for species-independent components. Per-species
 * code is in the Vector directory and namespace. */
class VectorModel : public TransmissionModel {
public:
  VectorModel(const scnXml::Vector vectorData, int populationSize);
  virtual ~VectorModel();
  
  /** Extra initialisation when not loading from a checkpoint, requiring
   * information from the human population structure. */
  virtual void init2 (const std::list<Host::Human>& population, int populationSize);
  
  virtual void scaleEIR (double factor);
//   virtual void scaleXML_EIR (scnXml::EntoData&, double factor) const;
  
  virtual TimeStep minPreinitDuration ();
  virtual TimeStep expectedInitDuration ();
  virtual TimeStep initIterate ();
  
  virtual void vectorUpdate (const std::list<Host::Human>& population, int populationSize);
  virtual void update (const std::list<Host::Human>& population, int populationSize);

  virtual double calculateEIR(PerHost& host, double ageYears); 
  
  virtual void setITNDescription ( const scnXml::ITNDescription& elt);
  virtual void setIRSDescription (const scnXml::IRS&);
  virtual void setVADescription (const scnXml::VectorDeterrent&);
  virtual void intervLarviciding (const scnXml::Larviciding::DescriptionType&);
  virtual void uninfectVectors();
  
  inline const ITNParams& getITNParams () const{
      return _ITNParams;
  }
  inline const IRSParams& getIRSParams() const{
      return _IRSParams;
  }
  virtual void summarize (Monitoring::Survey& survey);
  
protected:
    virtual void checkpoint (istream& stream);
    virtual void checkpoint (ostream& stream);
    
private:
  /** Return the index in speciesIndex of mosquito, throwing if not found. */
  size_t getSpeciesIndex (string mosquito)const {
    map<string,size_t>::const_iterator sIndex = speciesIndex.find (mosquito);
    if (sIndex == speciesIndex.end()) {
      ostringstream oss;
      oss << "Intervention description for anopheles species \""
	  << mosquito << "\": species not found in entomology description";
      throw util::xml_scenario_error(oss.str());
    }
    return sIndex->second;
  }
    
    /** Return the mean availability of human population to mosquitoes. */
    static double meanPopAvail (const std::list<Host::Human>& population, int populationSize);
    
    /** Confirm simulationMode allows use of interventions; throw if not. */
    void checkSimMode() const;
    
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
  
  
    /// Number of iterations performed during initialization, or negative when done.
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
   * Array will be recreated by constructor, but some members of VectorAnopheles
   * need to be checkpointed. */
  vector<Anopheles::AnophelesModel> species;
  
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
  
  /** Parameters used by ITN and IRS models. */
  ITNParams _ITNParams;
  /** ditto */
  IRSParams _IRSParams;
  
  friend class PerHost;
  friend class VectorAnophelesSuite;
};

} }
#endif
