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

#ifndef Hmod_proteome
#define Hmod_proteome

#include <string>
#include <vector>
#include <list>
#include "global.h"

using namespace std;

class Mutation;

//! Initialises the proteome module.
void initProteomeModule();

class Mutation;
class ProteinPosition;

/** Only used within proteome code. */
class Protein {
  string name;
  vector<ProteinPosition*> positions;
  
public:
  Protein(string _name);
  /** @brief Create, reading from checkpoint */
  Protein(istream& in);
  ~Protein();
  
  inline string getName() const {
    return name;
  }
  void addPosition(ProteinPosition* _position);
  Mutation* getMutation(int _position, char _allele) throw(int);
  void write (ostream& out);
};

/** Only used within proteome code. */
class ProteinPosition {
  Protein* protein;
  int position;
  char wildType;
  vector<Mutation*> mutations;
  
public:
  ProteinPosition(Protein* _protein, int _position, char _wildType);
  /** @brief Create, reading from checkpoint */
  ProteinPosition(Protein* _protein, istream& in);
  ~ProteinPosition();
  
  void addMutation(Mutation* _mutation);
  inline Protein* getProtein() {
    return protein;
  }
  inline string getProteinName() {
    return protein->getName();
  }
  Mutation* getMutation(char _allele) throw(int);
  inline int getPosition() const {
    return position;
  }
  void write (ostream& out);
};

/** Used by proteome and drug code. */
class Mutation {
  ProteinPosition* position;
  char allele;
  
public:
  Mutation(ProteinPosition* _position, char _allele);
  ~Mutation();
  inline string getProteinName() const {
    return position->getProteinName();
  }
  inline int getPosition() const {
    return position->getPosition();
  }
  inline char getAllele() const {
    return allele;
  }
  int operator==(const Mutation& rhs) const;
};

/** Each infection has an instance of this class. */
class ProteomeInstance {
  static int currentID;
  static vector<ProteomeInstance> instances;
  
  int proteomeID;
  //TODO: Fitness
  // List of mutations. DON'T try using a custom comparator with C++ set/map,
  // when we store pointers, because the comparator can't specify when two
  // elements are equal.
  list<Mutation*> mutations;
  
  // /** @brief Create, reading from checkpoint */
  //ProteomeInstance (istream& in);
  //void write (ostream& out);
  
public:
  /** Creates all unique instances of proteome. */
  static void init (Mutation*);
  
  /** Drug code needs a list of all instances. */
  static inline vector<ProteomeInstance> getInstances() {
    return instances;
  }
  
  //! For a new infection, randomly chooses and returns a proteome.
  static ProteomeInstance* newInfection();
  /// When loading a checkpoint, use the proteome ID to find the original proteome.
  static inline ProteomeInstance* getProteome(int proteome) {
    return &instances[proteome];
  }
  
  
  ProteomeInstance();
  ~ProteomeInstance();
  
  inline int getProteomeID() const {
    return proteomeID;
  }
  /** True if this ProteomeInstance has all mutations in _mutations. */
  bool hasMutations(vector<Mutation*> _mutations) const;
  
  friend class PkPdDrugSuite;
};

/** Methods getInfection and getProteome are used in Infection code. */
class ProteomeManager {
  static vector<Protein*> proteins;
  
public:
  static void addProtein(Protein* _protein);
  static Mutation* getMutation(string _proteinName, int _position, char _allele) throw(int);
  static void write (ostream& out);
  static void read (istream& in);
};
#endif
