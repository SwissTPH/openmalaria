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
#include "global.h"

using namespace std;


//! Initialises the proteome module.
void initProteomeModule();

class ProteinPosition;
class Mutation;

class Protein {
  string name;
  vector<ProteinPosition*> positions;

  public:
  Protein(string _name);
  /** @brief Create, reading from checkpoint */
  Protein(istream& in);
  ~Protein();
  
  string getName() const;
  void addPosition(ProteinPosition* _position);
  Mutation* getMutation(int _position, char _allele) throw(int);
  void write (ostream& out);
};

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
  Protein* getProtein();
  string getProteinName();
  Mutation* getMutation(char _allele) throw(int);
  int getPosition() const;
  void write (ostream& out);
};

class Mutation {
  ProteinPosition* position;
  char allele;

  public:
  Mutation(ProteinPosition* _position, char _allele);
  ~Mutation();
  string getProteinName() const;
  int getPosition() const;
  char getAllele() const;
  int operator==(const Mutation& rhs) const;
};

class ProteomeInstance {
  static int currentID;
  int proteomeID;
  //Fitness to be done
  vector<Mutation*> mutations;

  public:
  ProteomeInstance();
  /** @brief Create, reading from checkpoint */
  ProteomeInstance (istream& in);
  ~ProteomeInstance();
  void addMutation(Mutation* _mutation);
  int getProteomeID();
  bool hasMutations(vector<Mutation*> _mutations);
  void write (ostream& out);
};

class ProteomeManager {
  static vector<ProteomeInstance*> instances;
  static vector<Protein*> proteins;

  public:
  static void addInstance(ProteomeInstance* _instance);
  static void addProtein(Protein* _protein);
  //!returns the proteome of a new infection
  static ProteomeInstance* getInfection();
  static ProteomeInstance* getProteome(int proteome);
  static Mutation* getMutation(string _proteinName, int _position, char _allele) throw(int);
  static vector<ProteomeInstance*> getInstances();
  static void write (ostream& out);
  static void read (istream& in);
};
#endif
