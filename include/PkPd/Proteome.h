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

#ifndef Hmod_Proteome
#define Hmod_Proteome

#include "Global.h"

#include <string>
#include <vector>
#include <list>

using namespace std;

namespace OM { namespace PkPd {
    class Mutation;
    class ProteinPosition;

/* NOTE:
Could replace the vectors/lists with maps and rewrite as something like:
    class ProteinSet {
	map<name,ProteinPositionSet> proteins;
	...
    }
Advantage: makes it clear when elements are keys and that there are no duplicate keys
(without duplicating keys in data).

Or could just replace the whole thing with a list of mutations, with position and protein and
mutation names as enums.
*/

// Note: All data here (apart from pointer values) is reproducibly reloadable from the XML data, so
// checkpointing is not needed.

/** Only used within proteome code. */
class Protein {
  string name;
  vector<ProteinPosition*> positions;
  
public:
  Protein(string _name) : name(_name) {}
  ~Protein ();
  
  void addPosition(ProteinPosition* _position);
  Mutation* getMutation(int _position, char _allele);
  inline bool isNamed (const string _name) const {
      return name == _name;
  }
  inline bool operator==(const Protein& rhs) const {
      return name == rhs.name;
  }
};

/** Only used within proteome code. */
class ProteinPosition {
  Protein* protein;
  int position;
  char wildType;
  vector<Mutation*> mutations;
  
public:
  ProteinPosition(Protein* _protein, int _position, char _wildType);
  ~ProteinPosition();
  
  void addMutation(Mutation* _mutation);
  inline Protein* getProtein() {
    return protein;
  }
  Mutation* getMutation(char _allele);
  inline int getPosition() const {
    return position;
  }
  inline bool operator==(const ProteinPosition& rhs) const {
      return position == rhs.position && (*protein) == (*rhs.protein);
  }
};

/** Used by proteome and drug code. */
class Mutation {
  ProteinPosition* position;
  char allele;
  
public:
  Mutation(ProteinPosition* _position, char _allele);
  ~Mutation();
  
  inline int getPosition() const {
    return position->getPosition();
  }
  inline char getAllele() const {
    return allele;
  }
  inline bool operator==(const Mutation& rhs) const {
      return allele == rhs.allele && (*position) == (*rhs.position);
  }
};

/** Each infection has an instance of this class.
 *
 * Static data here currently doesn't need checkpointing. */
class ProteomeInstance {
  static int currentID;
  static vector<ProteomeInstance> instances;
  
  
  uint32_t proteomeID;
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
  static void cleanup ();
  
  /** Drug code needs a list of all instances. */
  static inline vector<ProteomeInstance>& getInstances() {
    return instances;
  }
  
  /** For a new infection, randomly chooses a proteome and returns a reference
   * to it.
   * 
   * Ownership is not passed; i.e. do not try to delete this pointer! */
  static const ProteomeInstance* newInfection();
  /// When loading a checkpoint, use the proteome ID to find the original proteome.
  static inline ProteomeInstance* getProteome(int proteome) {
    return &instances[proteome];
  }
  
  
  ProteomeInstance();
  ~ProteomeInstance();
  
  inline uint32_t getProteomeID() const {
    return proteomeID;
  }
  /** True if this ProteomeInstance has all mutations in _mutations. */
  bool hasMutations(vector<Mutation*> _mutations) const;
};

/** Manages the list of proteins and (through these) the mutations.
 *
 * Static data here is set-up directly from XML and doesn't need checkpointing. */
class ProteomeManager {
  static vector<Protein*> proteins;
  
public:
    //! Initialises the proteome module.
    static void init ();
    static void cleanup ();
    
    static void addProtein(Protein* _protein);
    static Mutation* getMutation(string _proteinName, int _position, char _allele);
};

} }
#endif