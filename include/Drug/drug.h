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

#ifndef Hmod_drug
#define Hmod_drug

#include <string>
#include <list>
#include <map>
#include <vector>
#include "global.h"
#include "proteome.h"

using namespace std;


class Human;

//! Initialises the drug module.
/*!
  \param withinHostTimeStep is the within host timestep in minutes.
  \param simulatorTimeStep is the simulator timestep in minutes.
 */
void initDrugModule(int _withinHostTimeStep, int _simulatorTimeStep);

//! Drug Dose.
class Dose {
  //!In minutes from start of simulator step
  int time;

  //! units?
  double quantity;

  public:
  ~Dose();
  Dose(const Dose &_original);
  Dose(int _time, double _quantity);
  Dose &operator=(const Dose &rhs);
  int operator==(const Dose &rhs) const;
  int operator<(const Dose &rhs) const;
  friend ostream& operator<<(ostream& out, const Dose &dose);
  friend istream& operator>>(istream& out, Dose &dose);
};

//! A class holding drug info.
/*!
  Holds all information pertaining drugs.

  For now there is a single class implementing a simple PK/PD model.

  In the future this will probably become abstract and subclasses be implemented.
 */
class Drug {
  //BEGIN Drug-type fields (same for all drugs of same type)
  //! The drug abbreviated name, used for registry lookups.
  string abbreviation;
  //! The drug name.
  string name; 
  //! Absorption factor.
  /*! Absorption = dose * factor / weight
   */
  double absorptionFactor;
  //! Half-life (in minutes)
  double halfLife;
  //! Pharma dynamic list of parameters.
  /*! A ordered list of required mutations.
   * The parameter value can be found on pdParameters.
   * The order is important, the first one takes precedence
   *   (a map cannot impement this).
   */
  vector< vector<Mutation*> >* requiredMutations;
  //! PD parameters (check requiredMutations)
  vector<double>* pdParameters;
  //! Fast data structure to know the PD param per proteome
  map<int,double>* proteomePDParameters;
  //END


  //Below here, fields should only be instantiated for humans.

  //! The human where it is used
  Human* human;
  //! A list (really a list) of doses.
  vector<Dose*> doses;
  //! Drug concentration (ng/mL ?).
  double _concentration;
  //! Drug concentration on the next cycle.
  double _nextConcentration;
  //!Used on human (reference to original drug structure)
  bool onHuman;


  public:
  Drug(const Drug &_original);
  Drug &operator=(const Drug &rhs);
  int operator==(const Drug &rhs) const;
  int operator<(const Drug &rhs) const;

  Drug(string _name, string _abbreviation,
      double _absorptionFactor, double _halfLife);
  ~Drug();
  friend ostream& operator<<(ostream& out, const Drug& drug);
  friend istream& operator>>(istream& in, Drug& drug);

  string getAbbreviation() const { return abbreviation;}
  double getAbsorptionFactor() const { return absorptionFactor;}
  double getHalfLife() const { return halfLife;}
  void setConcentration(double concentration);
  double getConcentration() const { return _concentration;}
  double getNextConcentration() const { return _nextConcentration;}
  void addConcentration(double concentration);
  double calculateDrugFactor(ProteomeInstance* infProteome) const;
  double calculateDecay(int _time) const;
  void decay();

  //! A new instance is returned for usage
  Drug use();

  //! Adds a PD Rule.
  /*! The order of rule adding is important! The first add should be the
   *  one with most mutations (typically the most resistant), the last
   *  one should be the sensitive (ie vector<mutation> = 0).
   */
  void addPDRule(vector<Mutation*> requiredMutations, double pdFactor);

  //! Parses the proteme instances.
  /*! Creates an association between ProteomeInstance and PD factor.
   *  This is solely for performance purposes.
   */
  void parseProteomeInstances();

};


//! The list of available drugs. Singleton, use getRegistry.
/*! This should really be a pointer to a class,
    and an instancer called on request from getDrug.
 */
class DrugRegistry {
  static vector<Drug*> drugs;
  static DrugRegistry* instance;

  DrugRegistry();
  friend ostream& operator<<(ostream& out, const DrugRegistry& registry);
  friend istream& operator>>(istream& in, DrugRegistry& registry);

  public:
  static DrugRegistry* getRegistry();

  //! Adds a new drug to the list
  static void addDrug(Drug* drug) throw (int);

  //!Returns a drug
  static Drug* getDrug(string abbreviation) throw (int);

};

//! Responsible interactions with the within host and clinical modules.
/*! Acts as a proxy, this has the following benefits:
 *    1. WH module only needs to call PD once (and not once for each drug)
 *    2. Ditto for general human maintenance of PK levels
 *    3. Can decide as to sinergy among drugs
 */
class DrugProxy {
  list<Drug*> _drugs;
  DrugRegistry* registry;
  double weight;	// human's weight

  public:
  DrugProxy();
  /// Destructor. NOTE: could it be a normal dtor?
  void destroy();
  
  //! Medicates an individual.
  /*! \param drugName - The drug abbreviation.
   *  \param qty      - the quantity (which units?).
   *  \param time     - Time in minutes since start of the simulation tStep.
   *
   *  Medicate has to be called in correct time order (ie first lower times).
   */
  void medicate(string _drugAbbrev, double _qty, int _time);
  double calculateDrugsFactor(ProteomeInstance* infProteome);
  void decayDrugs();

  void write (ostream& out) const;
  void read (istream& in);
  
  void setWeight (double w) {
    weight = w;
  }
};

#endif
