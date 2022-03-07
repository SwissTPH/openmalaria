import os
from math import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pandas as pd

names = {
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
	35: 'Vector_EIR_Input',
	36: 'Vector_EIR_Simulated',
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
	79:	'innoculationsPerVector'	
}

def load(scenario):
	df = pd.read_csv(scenario, sep="\t", header=None)
	df.columns = ['survey', 'group', 'code', 'value']
	return df

def extract(df, surveys, groups, code):
	data = []

	# use only survey 2 (post intervention)
	for survey in surveys:
		y = np.array(len(groups), dtype=float)
		for group in groups:
			v = df[(df['survey'] == survey) & (df['group'] == group) & (df['code'] == code)]['value'].values
			if len(v) > 1:
				print('Too many values found for survey', survey, ' code ', code, ' and group ', group)
			elif len(v) == 0:
				print('No values found for survey', survey, ' code ', code, ' and group ', group)
			else:
				y += v[0]
		data.append(y)

	return np.array(data)

def plot(df, codes, groups, labels, title, norm_code=None):
	fig = plt.figure(figsize=(6,4))
	ax = fig.add_subplot(1,1,1)

	for i in range(len(codes)):
		surveys = np.arange(max(df['survey']))+1
		data = extract(df, surveys, groups, codes[i])

		if norm_code is not None:
			nHost = extract(df, surveys, groups, norm_code)
			data /= nHost

		ax.plot(surveys, data, marker='o', linewidth=1, markersize=4, label=labels[i])
		
	ax.set_title(title)
	ax.set_xlabel("time steps")
	ax.set_ylabel(title)
	ax.legend(loc="upper left")
	plt.tight_layout()

df = load(f'output.txt')

def plot_by_group(df, codes, groups, labels, title, norm_code=None):
	fig = plt.figure(figsize=(6,4))
	ax = fig.add_subplot(1,1,1)

	for i in range(len(codes)):
		surveys = np.arange(max(df['survey']))+1
		for j in range(len(groups)):
			data = extract(df, surveys, [groups[j]], codes[i])

			if norm_code is not None:
				nHost = extract(df, surveys, [groups[j]], norm_code)
				data /= nHost

			ax.plot(surveys, data, marker='o', linewidth=1, markersize=4, label=labels[i][j])
		
	ax.set_title(title)
	ax.set_xlabel("time steps")
	ax.set_ylabel(title)
	ax.legend(loc="upper left")
	plt.tight_layout()

def plot_cts(df, field, label):
	fig = plt.figure(figsize=(6,4))
	ax = fig.add_subplot(1,1,1)

	data = df[field].values
	ax.plot(np.arange(len(data)), data, marker='o', markersize=2, label=label)

	ax.set_title(field)
	ax.set_xlabel("time steps")
	ax.set_ylabel(field)
	ax.legend(loc="upper right")

# Normal output 
# =============

# The file is organized in 4 columns:
# survey, group, code, value
# There is one value per survey, group and code
# Groups are typically age groups (specified in xml)
# Note 1: some outputs do not have groups
# Note 2: sometimes instead of age, groups will indicate phenotype or other things depending on the output code (see wiki)

# Example scenario has 4 surveys per year and 2 age groups (< 5y and < 90y).
# Monitoring period is 20 years. A GVI interventions is deployed after 10 years
# The transmission is seasonal and you can observe a peak every year
# Simulated EIR, new infections, Malaria episodes and others are expected to decrease after the intervention is deployed
# Note: increase population size and run multiple seeds for better results
df = load(f'output.txt')

# Entomological Innoculation Rate (input and simulated), no age groups
plot(df, codes=[35,36], groups=[0], labels=['Input EIR', 'Simulated EIR'], title="EIR")

# Number of Severe Episodes per person (divided by number of hosts, which is code 0)
plot(df, codes=[15], groups=[1,2], labels=['Severe Epiosdes'], title="Severe Episodes per person per survey", norm_code=0)

# Number of Uncomplicated Episodes per person and by age group
plot_by_group(df, codes=[14], groups=[1,2], labels=[['age < 5', 'age < 90']], title="Unomplicated Episodes per person per survey", norm_code=0)

# Continous output
# ================

# The file contains one column per output (specified in the xml file)

df = pd.read_csv(f'ctsout.txt', sep="\t", header=1)

plot_cts(df, ["immunity Y"], label="immunity Y")
plot_cts(df, ["immunity h"], label="immunity h")
plot_cts(df, ["new infections"], label="new infections")

plt.show()
