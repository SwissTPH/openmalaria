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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301,
 * USA.
 */

#ifndef Hmod_Transmission
#define Hmod_Transmission

#include "Transmission/NonVectorModel.h"
#include "Transmission/VectorModel.h"

namespace OM
{
namespace Transmission
{
///@brief Creation, destruction and checkpointing
//@{
/// Creates a derived class
static TransmissionModel *createTransmissionModel(const scnXml::Entomology &entoData, int populationSize)
{
  // Entomology contains either a list of at least one anopheles or a list of at
  // least one EIRDaily.
  const scnXml::Entomology::VectorOptional& vectorData = entoData.getVector();

  TransmissionModel *model;
  if (vectorData.present())
    model = new VectorModel(entoData, vectorData.get(), populationSize);
  else {
      const scnXml::Entomology::NonVectorOptional& nonVectorData = entoData.getNonVector();
    if (!nonVectorData.present())       // should be a validation error, but anyway...
      throw util::xml_scenario_error ("Neither vector nor non-vector data present in the XML!");
    if (util::ModelOptions::option( util::VECTOR_LIFE_CYCLE_MODEL ) ||
        util::ModelOptions::option( util::VECTOR_SIMPLE_MPD_MODEL ))
        throw util::xml_scenario_error("VECTOR_*_MODEL is only compatible with the vector model (and non-vector data is present).");
    model = new NonVectorModel(entoData, nonVectorData.get());
  }

  if( entoData.getScaledAnnualEIR().present() ){
      model->scaleEIR( entoData.getScaledAnnualEIR().get() / model->annualEIR );
      assert( util::vectors::approxEqual( model->annualEIR, entoData.getScaledAnnualEIR().get() ) );
  }

  if( util::CommandLine::option( util::CommandLine::PRINT_ANNUAL_EIR ) ){
      //Note: after internal scaling (which doesn't imply exit)
      //but before external scaling.
      cout << "Total annual EIR: "<<model->annualEIR<<endl;
  }

  return model;
}
}
}

#endif