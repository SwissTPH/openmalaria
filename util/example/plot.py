# This file is part of OpenMalaria.
# 
# Copyright (C) 2005-2025 Swiss Tropical and Public Health Institute
# Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
# Copyright (C) 2020-2025 University of Basel
# Copyright (C) 2025 The Kids Research Institute Australia
#
# OpenMalaria is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

from matplotlib import pyplot as plt
import pandas as pd

mm = {
	0: 'nHost',
	1: 'nInfect',
	2: 'nExpectd',
	3: 'nPatent',
	4: 'sumLogPyrogenThres',
	5: 'sumlogDens',
	6: 'totalInfs',
	7: 'nTransmit',
	8: 'totalPatentInf',
	9: 'contrib',
	10: 'sumPyrogenThresh',
	11: 'nTreatments1',
	12: 'nTreatments2',
	13: 'nTreatments3',
	14: 'nUncomp',
	15: 'nSevere',
	16: 'nSeq',
	17: 'nHospitalDeaths',
	18: 'nIndDeaths',
	19: 'nDirDeaths',
	20: 'nEPIVaccinations',
	21: 'allCauseIMR',
	22: 'nMassVaccinations',
	23: 'nHospitalRecovs',
	24: 'nHospitalSeqs',
	25: 'nIPTDoses',
	26: 'annAvgK',
	27: 'nNMFever',
	30: 'innoculationsPerAgeGroup',
	28: 'innoculationsPerDayOfYear',
	29: 'kappaPerDayOfYear',
	31: 'Vector_Nv0',
	32: 'Vector_Nv',
	33: 'Vector_Ov',
	34: 'Vector_Sv',
    35: 'InputEIR',
    36: 'SimulatedEIR',
	39: 'Clinical_RDTs',
	40: 'Clinical_DrugUsage',
	41: 'Clinical_FirstDayDeaths',
	42: 'Clinical_HospitalFirstDayDeaths',
	43: 'nNewInfections',
	44: 'nMassITNs',
	45: 'nEPI_ITNs',
	46: 'nMassIRS',
	47: 'nMassVA',
	48: 'Clinical_Microscopy',
	49: 'Clinical_DrugUsageIV',
	50: 'nAddedToCohort',
	51: 'nRemovedFromCohort',
	52: 'nMDAs',
	53: 'nNmfDeaths',
	54: 'nAntibioticTreatments',
	55: 'nMassScreenings',
	56: 'nMassGVI',
	57: 'nCtsIRS',
	58: 'nCtsGVI',
	59: 'nCtsMDA',
	60: 'nCtsScreenings',
	61: 'nSubPopRemovalTooOld',
	62: 'nSubPopRemovalFirstEvent',
	63: 'nPQTreatments',
	64: 'nTreatDiagnostics',
	65: 'nMassRecruitOnly',
	66: 'nCtsRecruitOnly',
	67: 'nTreatDeployments',
	68: 'sumAge',
	69: 'nInfectByGenotype',
	70: 'nPatentByGenotype',
	71: 'logDensByGenotype',
	72: 'nHostDrugConcNonZero',
	73: 'sumLogDrugConcNonZero',
	74:	'expectedDirectDeaths',
	75:	'expectedHospitalDeaths',
	76:	'expectedIndirectDeaths',
	77:	'expectedSequelae',
	78:	'expectedSevere',
	79:	'innoculationsPerVector',
	80:	'nSevereWithoutComorbidities',
	81:	'expectedSevereWithoutComorbidities'	
}

mmi = {v: k for k, v in mm.items()}

df = pd.read_csv(f'output.txt', sep="\t", header=1)
df.columns = ['survey', 'age-group', 'measure', 'value']

nHosts = df[(df["measure"] == mmi['nHost'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()#.values
inputEIR = df[(df["measure"] == mmi['InputEIR'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()
simulatedEIR = df[(df["measure"] == mmi['SimulatedEIR'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()
innoculations = df[(df["measure"] == mmi['innoculationsPerAgeGroup'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()#.values
nPatent = df[(df["measure"] == mmi['nPatent'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()#.values
nUncomp = df[(df["measure"] == mmi['nUncomp'])].groupby(['survey', 'measure']).sum().value[1:].reset_index()#.values

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(inputEIR.index, inputEIR.value, marker="v", markersize=1, label=f"inputEIR")
ax.plot(simulatedEIR.index, simulatedEIR.value, marker="v", markersize=1, label=f"simulatedEIR")
ax.legend(loc="upper left")

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
prevalence = nPatent / nHosts
ax.plot(prevalence.index, prevalence.value, marker="v", markersize=1, label=f"Prevalence")
ax.legend(loc="upper left")

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ci = nUncomp / nHosts
ax.plot(ci.index, ci.value, marker="v", markersize=1, label=f"Clinical Incidence")
ax.legend(loc="upper left")

plt.show()
