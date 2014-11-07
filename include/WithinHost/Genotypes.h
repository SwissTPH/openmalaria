/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2014 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2014 Liverpool School Of Tropical Medicine
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

#ifndef Hmod_WH_Genotypes
#define Hmod_WH_Genotypes

namespace OM { namespace WithinHost {

/** Represents infection genotypes. */
class Genotypes {
public:
    /// Represent a combination of alleles, each from a different locus
    struct Genotype{
        Genotype( uint32_t allele, double iFreq, double fit ):
            init_freq(iFreq), fitness(fit)
        {
            alleles.insert(allele);
        }
        Genotype cross( const Genotype& that )const;
        set<uint32_t> alleles;      // set of codes of all alleles
        double init_freq;
        double fitness;
    };
    
    /** Initialise with a single genotype. */
    static void initSingle();
    
    /** Initialise from XML data. Call this before other static methods are
     * used (from PK/PD code). */
    static void init( const scnXml::Scenario& scenario );
    
    /** Map a locus name and allele name to an allele code.
     * 
     * Note that two alleles from different loci will always have different
     * codes.
     * 
     * Returns max value when no match is found. */
    static uint32_t findAlleleCode( const string& locus, const string& allele );
    
    /** Get a reference to the list of all genotypes. */
    static const vector<Genotype>& getGenotypes();
    
    /** Sample the genotype using the configured approach. */
    static uint32_t sampleGenotype();
    
    /** Get the number of genotypes. Functions like sampleGenotype use values
     * from 0 to one less than this. */
    inline static size_t N(){ return N_genotypes; }
    
    /** Get the initial frequency of some genotype. */
    static double initialFreq( size_t genotype );
    
private:
    static size_t N_genotypes;
};

}
}
#endif
