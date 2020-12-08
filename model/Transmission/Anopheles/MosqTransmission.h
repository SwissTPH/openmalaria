// /* This file is part of OpenMalaria.
//  * 
//  * Copyright (C) 2005-2015 Swiss Tropical and Public Health Institute
//  * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
//  * 
//  * OpenMalaria is free software; you can redistribute it and/or modify
//  * it under the terms of the GNU General Public License as published by
//  * the Free Software Foundation; either version 2 of the License, or (at
//  * your option) any later version.
//  * 
//  * This program is distributed in the hope that it will be useful, but
//  * WITHOUT ANY WARRANTY; without even the implied warranty of
//  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  * General Public License for more details.
//  * 
//  * You should have received a copy of the GNU General Public License
//  * along with this program; if not, write to the Free Software
//  * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//  */

// #ifndef Hmod_Anopheles_Transmission
// #define Hmod_Anopheles_Transmission

// #include "Global.h"
// #include "Transmission/Anopheles/EmergenceModel.h"
// #include "schema/entomology.h"

// #include <limits>

// namespace OM {
// namespace Transmission {
// namespace Anopheles {
// using util::vecDay2D;

// // enumeration of gettable stats for cts out
// enum VecStat { PA, PAmu, PA1, PAh, PDF, PDIF, NV, OV, SV };

// /** Encapsulates the central part of the Chitnis et al transmission model:
//  * vector transmission of malaria.
//  * 
//  * This is only part of the model; the VectorModel class is largely just a
//  * wrapper to support multiple mosquito species, and the AnophelesModel
//  * class adds parameter initialisation and intervention support to this class
//  * as well as translating between the (1- or 5-day) time-steps used by the
//  * simulator and the 1-day time-step used by this model.
//  */
// class MosqTransmission {
// public:
//     ///@brief initialisation functions
//     //@{
//     MosqTransmission();
    
//     MosqTransmission (MosqTransmission&& o):
//             emergence(move(o.emergence)),
//             mosqRestDuration(move(o.mosqRestDuration)),
//             EIPDuration(move(o.EIPDuration)),
//             N_v_length(move(o.N_v_length)),
//             minInfectedThreshold(move(o.minInfectedThreshold)),
//             P_A(move(o.P_A)),
//             P_df(move(o.P_df)),
//             P_dif(move(o.P_dif)),
//             P_dff(move(o.P_dff)),
//             N_v(move(o.N_v)),
//             S_v(move(o.S_v)),
//             fArray(move(o.fArray)),
//             ftauArray(move(o.ftauArray)),
//             uninfected_v(move(o.uninfected_v)),
//             timeStep_N_v0(move(o.timeStep_N_v0))
//     {}
    
//     void operator= (MosqTransmission&& o) {
//             emergence = move(o.emergence);
//             mosqRestDuration = move(o.mosqRestDuration);
//             EIPDuration = move(o.EIPDuration);
//             N_v_length = move(o.N_v_length);
//             minInfectedThreshold = move(o.minInfectedThreshold);
//             P_A = move(o.P_A);
//             P_df = move(o.P_df);
//             P_dif = move(o.P_dif);
//             P_dff = move(o.P_dff);
//             N_v = move(o.N_v);
//             S_v = move(o.S_v);
//             fArray = move(o.fArray);
//             ftauArray = move(o.ftauArray);
//             uninfected_v = move(o.uninfected_v);
//             timeStep_N_v0 = move(o.timeStep_N_v0);
//     }
    
//     /** Initialise parameters and variables.
//      * 
//      * This is only a fraction of parameter initialisation; see also
//      * AnophelesModel::initialise. */
//     void initialise ( const scnXml::AnophelesParams::LifeCycleOptional& lcOpt, const scnXml::AnophelesParams::SimpleMPDOptional& simpleMPDOpt, const scnXml::Mosq& mosq );
    
//     /** (Re) allocate and initialise some state variables. Must be called
//      * before model is run. */
//     void initState ( double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah,
//                      double tsP_df, double tsP_dff,
//                      double initNvFromSv, double initOvFromSv,
//                      const vecDay<double>& forcedS_v);
    
//     /// Helper function for initialisation.
//     void initIterateScale ( double factor );
    
//     /** Set up the non-host-specific interventions. */
//     inline void initVectorInterv( const scnXml::VectorSpeciesIntervention& elt, size_t instance ){
//         emergence->initVectorInterv( elt, instance ); }
//     //@}
    
//     /** Update by one day (may be called multiple times for 1 time-step update).
//      * 
//      * @param d0 Time of the start of the day-long update period
//      * @param tsP_A P_A for this time-step
//      * @param tsP_df P_df for this time-step
//      * @param tsP_dif P_dif for this time-step, per parasite genotype
//      * @param tsP_dff P_dff for this time step
//      * @param partialEIR Vector, per genotype; after calculation, the latest
//      *  S_v values are multiplied by EIR_factor and added to this.
//      * @param EIR_factor see parameter partialEIR
//      */
//     void update( SimTime d0, double tsP_A, double tsP_Amu, double tsP_A1, double tsP_Ah, double tsP_df,
//                    const vector<double> tsP_dif, double tsP_dff,
//                    bool isDynamic,
//                    vector<double>& partialEIR, double EIR_factor);
    
//     ///@brief Interventions and reporting
//     //@{
//     void uninfectVectors();
//     //@}
    
//     inline SimTime getEIPDuration() const {
//         return EIPDuration;
//     }
    
//     ///@brief Functions used in reporting
//     //@{
//     /// Reset per-time-step statistics before running time-step updates
//     inline void resetTSStats() {
//         timeStep_N_v0 = 0.0;
//     }
//     /// Get mean emergence per day during last time-step
//     inline double getLastN_v0 () const{
//         return timeStep_N_v0 / SimTime::oneTS().inDays();
//     }
//     /// Get mean P_A/P_df/P_dif/N_v/O_v/S_v during last time-step
//     /// @param vs PA, PDF, PDIF, NV, OV or SV
//     double getLastVecStat( VecStat vs )const;
    
//     inline SimTime getMosqRestDuration() const {
//         return mosqRestDuration;
//     }
    
//     inline double getResAvailability() const{
//         return emergence->getResAvailability();
//     }
//     inline double getResRequirements() const{
//         return emergence->getResRequirements();
//     }
    
//     /// Write some per-species summary information.
//     void summarize( size_t species )const;
//     //@}

    
// private:

// };

// }
// }
// }
// #endif
