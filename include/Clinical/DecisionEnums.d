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

/* This is BOTH a C(++) header file and a D import file!
 * As such, there are quite a few limitations:
 * Don't use macros
 *   hence no #ifndef HEADER_IDENT ... #endif
 *   which is OK (although potentially inefficient), since nothing is imported
 * No namespaces
 * Don't set the D module name
 * Filename must have .d extension */

/** Enumeration flags specifying decisions within the case-management tree.
 *
 * Not all may be present, if no decision was taken in the relevent part of the
 * tree, or if the information was filtered out. It is not directly possible to
 * tell the difference between a decision whose "flag" has value 0 and that
 * decision not having been taken, but for instance quality and adherence
 * decisions should have taken place except when drug is QN or NO_AM.
 * 
 * Decisions are:
 *  Testing method (or no test),
 *  Result of testing for parasites,
 *  Anti-malarial drug prescribed,
 *  Quality of prescribed AM (anti-malarial),
 *  Adherence of patient taking prescribed AM,
 *  Quality of case management (applicable to severe cases),
 *  Treatment type for severe cases (inc. in/outside a hospital),
 *  Treatment seeking delay (usually 0-2 days).
 * 
 * The MASK values are used to filter bits associated with a decision; e.g.
 * (x & QUALITY_MASK) can only have values QUALITY_BAD or QUALITY_GOOD.
 * 
 * TSDELAY is slightly different; (x & TSDELAY_MASK) << TSDELAY_LSHIFT yields
 * an integer in the range [0,15] which is the delay in days. */
enum DecisionEnums {
    NONE		= 0x0,
    
    TEST_NONE		= 0x0,
    TEST_MICROSCOPY	= 0x1,
    TEST_RDT		= 0x2,
    TEST_MASK		= 0xF,
    
    RESULT_TRUE_POS	= 0x10,
    RESULT_FALSE_NEG	= 0x20,
    RESULT_TRUE_NEG	= 0x30,
    RESULT_FALSE_POS	= 0x40,
    RESULT_TRUE		= 0x50,		///< true positive but malaria not cause of sickness
    RESULT_FALSE	= 0x60,		///< false negative where malaria is not cause
    RESULT_MASK		= 0xF0,
    
    DRUG_NO_AM		= 0x000,	///< no anti-malarial drug
    DRUG_SP		= 0x100,
    DRUG_AL		= 0x200,
    DRUG_QN		= 0x300,	///< value can be changed, but (value-1) won't match desicionID in VC's case management tree
    DRUG_MASK		= 0xF00,
    
    QUALITY_BAD		= 0x0000,
    QUALITY_GOOD	= 0x1000,
    QUALITY_MASK	= 0x1000,
    
    ADHERENCE_BAD	= 0x0000,
    ADHERENCE_GOOD	= 0x2000,
    ADHERENCE_MASK	= 0x2000,
    
    MANAGEMENT_BAD	= 0x0000,
    MANAGEMENT_GOOD	= 0x4000,
    MANAGEMENT_MASK	= 0x4000,
    
    TREATMENT_NO_AM	= 0x0,
    TREATMENT_PREREF	= 0x10000,
    TREATMENT_DELREF	= 0x20000,
    TREATMENT_PREREF_IMMREF	= 0x50000,
    TREATMENT_PREREF_DELREF	= 0x30000,
    TREATMENT_PARENTAL	= 0x80000,
    TREATMENT_MASK	= 0xF0000,
    
    TSDELAY_LSHIFT	= 20,
    TSDELAY_MASK	= 0xF00000,
};
