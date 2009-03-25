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

#ifndef Hmod_oldITNWHM
#define Hmod_oldITNWHM

#include "oldWithinHostModel.h"

class OldITNWithinHostModel : public OldWithinHostModel {
public:
  OldITNWithinHostModel () : _SPattenuationt(0) {}
  
protected:
  /*!  SP drug action applies to each infection depending on genotype and when
  the individual had their last dose of SP */
  void SPAction(Human&);
  
  void IPTattenuateAsexualDensity (std::list<DescriptiveInfection>::iterator i);
  void IPTattenuateAsexualMinTotalDensity (Human&);
  
  virtual void write(ostream& out) const;
  virtual void read(istream& in);
  
private:
  //! time at which attenuated infection 'would' end if SP present
  int _SPattenuationt;
};

#endif
