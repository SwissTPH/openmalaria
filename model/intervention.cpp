#include "intervention.h"
#include "inputData.h"
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
/*

Initializing intervention data.

TODO: now that we have enums in the c code, merge these constant with the c code
Encoding them as 2powers was not consistent with the summaryOption constants.
Check the binary encoding works
*/
   int vaccineType;
  
 
   double *targetagetstep;
     int targetagetstepX;
   double *vaccineCoverage;
     int vaccineCoverageX;
   double *PEVInitialMeanEfficacy;
     int PEVInitialMeanEfficacyX;
   double *BSVInitialMeanEfficacy;
     int BSVInitialMeanEfficacyX;
   double *TBVInitialMeanEfficacy;
     int TBVInitialMeanEfficacyX;
   double PEVefficacyB;
   double BSVefficacyB;
   double TBVefficacyB;
   double PEVdecay;
   double BSVdecay;
   double TBVdecay;
   short ITN;
   double ITNdecay;
   double Pu0;
   double Pu1;
   double c;
   double z;
   short IPT;
   int numberOfIPTiDoses;
   int *iptiTargetagetstep;
     int iptiTargetagetstepX;
   double *iptiCoverage;
     int iptiCoverageX;
   double iptiEffect;
   int numberOfGenoTypes;
   double *genotypeFreq;
     int genotypeFreqX;
   double *genotypeACR;
     int genotypeACRX;
   int *genotypeProph;
     int genotypeProphX;
   int *genotypeTolPeriod;
     int genotypeTolPeriodX;
   double *genotypeAtten;
     int genotypeAttenX;



  /*!
      Commen to all vaccine types. Number of vaccine doses that are given either
      through EPI or as EPI Boosters
   */
 int _numberOfEpiDoses;

  /*!
    Number of initial efficacy values in the XML. If more than this number or
    doses are given, we assume the vaccine efficacy goes to the last defined
    initial efficacy.
  */
 int _numberOfInitEff;

// Intervention::Intervention() {

//     int isIpti;
//     isIpti=get_is_ipti();
//     if ( isIpti ==  1) {
//         IPT=1/* bool */;
//     }
//     else {
//         IPT=0/* bool */;
//     }
//     vaccineType=get_vaccine_type();
//     if ( vaccineType >  0) {
//         initVaccineParameters();
//     }
//     if ( ITN) {
//         initITNParameters();
//     }
//     if ( IPT) {
//         initIPTIParameters();
//     }
// }

void initVaccineParameters(){

    //Read in vaccine specifications
    _numberOfEpiDoses=get_number_of_epi_doses();
    _numberOfInitEff=get_number_of_init_eff();
    if ( get_vaccine_halflife_yrs(preerythrocytic_reduces_h) >  0) {
        PEVdecay=log(2.0)/((get_vaccine_halflife_yrs(preerythrocytic_reduces_h))*daysInYear/(1.0*interval));
    }
    else {
        PEVdecay=0;
    }
    if ( get_vaccine_halflife_yrs(erythrocytic_reduces_y) >  0) {
        BSVdecay=log(2.0)/((get_vaccine_halflife_yrs(erythrocytic_reduces_y))*daysInYear/(1.0*interval));
    }
    else {
        BSVdecay=0;
    }
    if ( get_vaccine_halflife_yrs(transmission_blocking_reduces_k) >  0) {
        TBVdecay=log(2.0)/((get_vaccine_halflife_yrs(transmission_blocking_reduces_k))*daysInYear/(1.0*interval));
    }
    else {
        TBVdecay=0;
    }
    PEVefficacyB=get_efficacy_b(preerythrocytic_reduces_h);
    BSVefficacyB=get_efficacy_b(erythrocytic_reduces_y);
    TBVefficacyB=get_efficacy_b(transmission_blocking_reduces_k);
    PEVInitialMeanEfficacy = (double*)malloc(((_numberOfInitEff))*sizeof(double));
    PEVInitialMeanEfficacyX = _numberOfInitEff;
    BSVInitialMeanEfficacy = (double*)malloc(((_numberOfInitEff))*sizeof(double));
    BSVInitialMeanEfficacyX = _numberOfInitEff;
    TBVInitialMeanEfficacy = (double*)malloc(((_numberOfInitEff))*sizeof(double));
    TBVInitialMeanEfficacyX = _numberOfInitEff;
    targetagetstep = (double*)malloc(((_numberOfEpiDoses))*sizeof(double));
    targetagetstepX = _numberOfEpiDoses;
    vaccineCoverage = (double*)malloc(((_numberOfEpiDoses))*sizeof(double));
    vaccineCoverageX = _numberOfEpiDoses;
    for ( int i=0;i<_numberOfInitEff; i++) {
        PEVInitialMeanEfficacy[i]=get_efficacy(preerythrocytic_reduces_h, i);
        BSVInitialMeanEfficacy[i]=get_efficacy(erythrocytic_reduces_y, i);
        TBVInitialMeanEfficacy[i]=get_efficacy(transmission_blocking_reduces_k, i);
    }
    for ( int i=0;i<_numberOfEpiDoses; i++) {
        targetagetstep[i]=floor(get_target_age_yrs(i)*daysInYear/(1.0*interval));
        vaccineCoverage[i]=get_coverage_epi_vaccine(i);
    }
}

void initITNParameters () {

    //TODO: Get ITN from the cpp code (and the other parameters)
    ITN=0/* bool */;
    if ( ITN) {
        Pu0=get_pu0();
        Pu1=get_pu1();
        c=get_sporogony_gonotrophy();
    }
}

void initIPTIParameters () {
    //TODO: 

    int i;
    numberOfGenoTypes=get_number_of_genotypes();
    iptiEffect=get_ipti_effect();
    genotypeFreq = (double*)malloc(((numberOfGenoTypes))*sizeof(double));
    genotypeFreqX = numberOfGenoTypes;
    genotypeACR = (double*)malloc(((numberOfGenoTypes))*sizeof(double));
    genotypeACRX = numberOfGenoTypes;
    genotypeProph = (int*)malloc(((numberOfGenoTypes))*sizeof(int));
    genotypeProphX = numberOfGenoTypes;
    genotypeTolPeriod = (int*)malloc(((numberOfGenoTypes))*sizeof(int));
    genotypeTolPeriodX = numberOfGenoTypes;
    genotypeAtten = (double*)malloc(((numberOfGenoTypes))*sizeof(double));
    genotypeAttenX = numberOfGenoTypes;
    for ( i=1;i<=numberOfGenoTypes; i++) {
        genotypeFreq[i - 1]=get_genotype_freq(i);
        genotypeACR[i - 1]=get_genotype_acr(i);
        genotypeProph[i - 1]=get_genotype_proph(i);
        genotypeTolPeriod[i - 1]=get_genotype_tolperiod(i);
        genotypeAtten[i - 1]=get_genotype_atten(i);
    }
    numberOfIPTiDoses=get_number_of_ipti_doses();
    iptiTargetagetstep = (int*)malloc(((numberOfIPTiDoses))*sizeof(int));
    iptiTargetagetstepX = numberOfIPTiDoses;
    iptiCoverage = (double*)malloc(((numberOfIPTiDoses))*sizeof(double));
    iptiCoverageX = numberOfIPTiDoses;
    for ( i=0;i<numberOfIPTiDoses; i++) {
        iptiTargetagetstep[i]=(int)floor(get_ipti_target_age_yrs(i)*daysInYear/(1.0*interval));
        iptiCoverage[i]=get_ipti_coverage(i);
    }
}

void clearVaccineParameters () {

    free(PEVInitialMeanEfficacy);
    free(BSVInitialMeanEfficacy);
    free(TBVInitialMeanEfficacy);
    free(targetagetstep);
    free(vaccineCoverage);
}

void clearITNParameters () {
    //TODO if so

}

void clearIPTIParameters () {

    free(genotypeFreq);
    free(genotypeACR);
    free(genotypeProph);
    free(genotypeTolPeriod);
    free(genotypeAtten);
    free(iptiTargetagetstep);
    free(iptiCoverage);
}

void initInterventionParameters () {

    int isIpti;
    isIpti=get_is_ipti();
    if ( isIpti ==  1) {
        IPT=1/* bool */;
    }
    else {
        IPT=0/* bool */;
    }
    vaccineType=get_vaccine_type();
    if ( vaccineType >  0) {
        initVaccineParameters();
    }
    if ( ITN) {
        initITNParameters();
    }
    if ( IPT) {
        initIPTIParameters();
    }
}

void clearInterventionParameters () {

    if ( vaccineType >  0) {
        clearVaccineParameters();
    }
    if ( ITN) {
        clearITNParameters();
    }
    if ( IPT) {
        clearIPTIParameters();
    }
}


// Intervention::~Intervention() {

//     if ( vaccineType >  0) {
//         clearVaccineParameters();
//     }
//     if ( ITN) {
//         clearITNParameters();
//     }
//     if ( IPT) {
//         clearIPTIParameters();
//     }
// }

