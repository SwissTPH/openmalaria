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
 *
 * This allows the enum values to be shared so there's no possibility of them
 * getting out-of-sync between the converter and the simulator.
 * 
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
 * Treatment type names aren't yet final.
 * 
 * TSDELAY is slightly different; (x & TSDELAY_MASK) >> TSDELAY_SHIFT yields
 * an integer in the range [0,15] which is the delay in days. */
enum DecisionEnums {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte
     * (note & | ^ are C++'s binary AND, OR and XOR operators). */
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
    
    //TODO: completely rename (w/wo pre-ref, w/wo hospital, hospital delay (0,1 or 2?)
    TREATMENT_NO_AM	= 0x0,
    TREATMENT_PREREF	= 0x10000,
    TREATMENT_DELREF	= 0x20000,
    TREATMENT_PREREF_IMMREF	= 0x30000,
    TREATMENT_PREREF_DELREF	= 0x40000,
    TREATMENT_PARENTAL	= 0x50000,
    TREATMENT_DELREF_GOOD	= 0x60000,
    TREATMENT_PREREF_IMMREF_GOOD	= 0x70000,
    TREATMENT_PREREF_DELREF_GOOD	= 0x80000,
    TREATMENT_PARENTAL_GOOD	= 0x90000,
    TREATMENT_NUM_TYPES	= 10,
    TREATMENT_SHIFT	= 16,
    TREATMENT_MASK	= 0xF0000,
    
    TSDELAY_NUM_MAX	= 3,
    TSDELAY_SHIFT	= 20,
    TSDELAY_MASK	= 0xF00000,
};
