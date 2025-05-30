/* This file is part of OpenMalaria.
 *
 * Copyright (C) 2025 The Kids Research Institute Australia
 * Copyright (C) 2025 Swiss Tropical and Public Health Institute
 * Copyright (C) 2025 University of Basel
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

#ifndef Hmod_util_XMLChecker
#define Hmod_util_XMLChecker

#include "Global.h"
#include <schema/scenario.h>

namespace OM
{
    namespace util
    {
        class XMLChecker
        {
        public:
            XMLChecker(){}
            ~XMLChecker(){}

            /*
            * Should be called after the input XML has been loaded and before all subsequent
            * initialisation. I.e. well before the simulation starts running.
            *
            * Performs custom checks on the input XML which are not enforced in the schema itself.
            *
            * Throws iff a check fails.
            *
            * The purpose of these this XML checking is to indentify certain problems in the input XML
            * as early as possible, so that the program can be exited (this is a job for the caller)
            * instead of execution time being wasted on a simulation that will ultimately crash
            * (or worse, produce meaningless results).
            *
            * Of course, it would be better to handle as many such issues as possible in the schema
            * itself, rather than here.  However, some required checks are not possible to do using
            * the version of the XSD spec supported by the XML validation library used by OpenMalaria.
            */
            void PerformPostValidationChecks(const scnXml::Scenario& scenario)
            {
                CheckModelOptionsAndParams();
            }
        private:

            /*
            * Verifies that either:
            *
            *  1) Both model options and parameters are explicitly stated in the input XML but not a
            *     model name.
            *
            *  or
            *
            *  2) At least a model name (e.g. a name referring to the base model) *is* explicitly stated.
            *     (Optionally, model options and/or parameters may also be stated.)
            */
            void CheckModelOptionsAndParams()
            {
                // TODO
            }
        };
    }
}

#endif
