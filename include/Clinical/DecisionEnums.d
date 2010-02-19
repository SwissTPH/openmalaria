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
 *  Result of testing for parasites (outcome of test; may not be accurate),
 *  Anti-malarial drug prescribed,
 *  Adherence of patient taking prescribed AM (UC),
 *  Quality of prescribed AM (anti-malarial) (UC),
 *  Treatment type for severe cases (inc. in/outside a hospital, pre-referals, quality of case management) (severe),
 *  Treatment seeking delay (usually 0-2 days) (UC).
 * 
 * The MASK values are used to filter bits associated with a decision; e.g.
 * (x & QUALITY_MASK) can only have values QUALITY_BAD or QUALITY_GOOD.
 * 
 * TSDELAY is slightly different; (x & TSDELAY_MASK) >> TSDELAY_SHIFT yields
 * an integer in the range [0,15] which is the delay in days. */
/* TODO:
One or two new branch-points in trees.
What of this information do we need to fix though?

Tree path information
-------------------

Needs to be known by code:
UC1/UC2/severe input
age-over-5 input
test type output
result input/output
treatment seeking delay output
in/out of hospital output

Used for medicate-lookup:
Drug type
adherance
quality
treatment delay

Other info, not needed outside tree:
Some of treatment seeking stuff?
*/
enum DecisionEnums {
    /* Values here are written in hexadecimal: http://en.wikipedia.org/wiki/Hexadecimal
     * Many are designed to be "flags", so the value corresponds to a single bit:
     * http://en.wikipedia.org/wiki/Flag_byte
     * (note & | ^ are C++'s binary AND, OR and XOR operators).
     * Maximum value I'd recommend using as a flag: 0x4000_0000 */
    NONE				= 0x0,
    
    /** MUST correspond to Pathogenesis::MORBIDITY_MASK.
     *
     * Pathogenesis constants from Constant.h within this mask may be used. */
    MORBIDITY_MASK		= 0x3F,
    
    AGE_OVER5			= 0x80,
    
    TEST_NONE			= 0x0,
    TEST_MICROSCOPY	= 0x100,
    TEST_RDT			= 0x200,
    TEST_MASK			= 0x300,
    
    RESULT_POSITIVE	= 0x400,
    RESULT_NEGATIVE	= 0x800,
    RESULT_DETERMINE	= 0xC00,
    RESULT_MASK		= 0xC00,
    
    DRUG_NO_AM		= 0x0000,	///< no anti-malarial drug
    DRUG_SP			= 0x1000,
    DRUG_AL			= 0x2000,
    DRUG_IV_QN		= 0x8000,	///< Quinine
    DRUG_IV_AS		= 0x9000,	///< Artesunate
    DRUG_FIRST_SEV	= 0x8000,
    DRUG_NUM_SEV		= 2,		///< Number of sequenctial severe drugs listed
    DRUG_MASK		= 0xF000,
    
    ADHERENCE_FULL	= 0x00000,
    ADHERENCE_MISSED_FIRST	= 0x10000,	///< 1st or 2nd dose missed
    ADHERENCE_MISSED_LAST	= 0x20000,
    ADHERENCE_SELECTIVE	= 0x30000,	///< TODO: define exactly what this means
    ADHERENCE_MASK	= 0x30000,
    
    QUALITY_GOOD		= 0x00000,
    QUALITY_BAD		= 0x80000,
    QUALITY_MASK		= 0x80000,
    
    // Following is organised into one block with 16 potential values:
    TREATMENT_NO_HOSPITAL	= 0x0,	// no hospitalization: no seeking or preref only
    TREATMENT_HOSPITAL	= 0x100000,	// immediate hospitalization
    TREATMENT_DEL_HOSPITAL	= 0x200000,	// hospitalization one day later
    TREATMENT_UNUSED_TYPE = 0x300000,	// could be used in the future
    TREATMENT_PREREF		= 0x400000,	// combine with above for pre-referal suppository
    TREATMENT_CM_GOOD	= 0x800000,	// combine with above for good case management in hospital
    TREATMENT_NUM_TYPES	= 16,
    TREATMENT_MASK	= 0xF00000,
    
    TSDELAY_NUM_MAX	= 3,
    TSDELAY_SHIFT	= 24,
    TSDELAY_MASK	= 0x3000000,
};
