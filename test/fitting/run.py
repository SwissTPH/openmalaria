import os
import sys
import math
import shutil
import numpy as np
import subprocess
import pandas as pd

FILES = [
'wu_5000_158',  'wu_5000_317',  'wu_5000_422',  'wu_5000_515',
'wu_5000_167',  'wu_5000_318',  'wu_5000_426',  'wu_5000_516',
'wu_5000_173',  'wu_5000_31',   'wu_5000_49',   'wu_5000_517',
'wu_5000_176',  'wu_5000_326',  'wu_5000_501',  'wu_5000_518',
'wu_5000_232',  'wu_5000_327',  'wu_5000_502',  'wu_5000_519',
'wu_5000_233',  'wu_5000_34',   'wu_5000_503',  'wu_5000_520',
'wu_5000_234',  'wu_5000_35',   'wu_5000_504',  'wu_5000_521',
'wu_5000_24',   'wu_5000_401',  'wu_5000_505',  'wu_5000_522',
'wu_5000_28',   'wu_5000_402',  'wu_5000_506',  'wu_5000_523',
'wu_5000_29',   'wu_5000_408',  'wu_5000_507',  'wu_5000_524',
'wu_5000_301',  'wu_5000_411',  'wu_5000_508',  'wu_5000_525',
'wu_5000_302',  'wu_5000_414',  'wu_5000_509',  'wu_5000_526',
'wu_5000_303',  'wu_5000_415',  'wu_5000_510',  'wu_5000_527',
'wu_5000_30',   'wu_5000_416',  'wu_5000_511',
'wu_5000_312',  'wu_5000_417',  'wu_5000_512',
'wu_5000_316',  'wu_5000_418',  'wu_5000_514']

def run(command, capture = False):
    result = subprocess.run(command, shell = True, capture_output = capture)#command) #, shell = True, capture_output = True)

    if result.returncode != 0:
        print(f'STDOUT: '+(result.stdout.decode() if result.stdout is not None else ''))
        print(f'STDERR: '+(result.stderr.decode() if result.stderr is not None else ''))
        return False
    else:
        return True

def runScenarios():
	currentVersion = 35
	newVersion = 41
	# newVersion = 38

	workdir = 'check/'

	parametersFile = 'parameters'
	parametersSetFile = 'parameters_default'

	shutil.rmtree(workdir, ignore_errors = True)
	os.makedirs(os.path.relpath(workdir), exist_ok=True)

	shutil.copy('../densities.csv', workdir)
	shutil.copy('../../build/openMalaria', workdir)
	shutil.copy('../../build/schema/scenario_current.xsd', workdir)
	shutil.copy(workdir+'scenario_current.xsd', workdir+'scenario_41.xsd')

	commandsFile = workdir+"commands.cmd"

	n = 0
	with open(parametersFile, "r") as fp:
		params = fp.read()
		pSet = np.loadtxt(parametersSetFile)

		for xmlFile in FILES:
			for i in range(0, len(pSet)):
				key = f'PARAM{i+1:02d}'
				val = pSet[i]
				params = params.replace(f'{key}', str(val))

			s = 1
			os.makedirs(os.path.relpath(workdir+str(s)), exist_ok=True)

			xmlParams = params.replace('@SEED@', str(s))
			with open(xmlFile+'.xml', "r") as fx:
			    xml = fx.read()
			    xml = f'{xml}{xmlParams}'
			    xml = xml.replace(f'schemaVersion="{currentVersion}"', f'schemaVersion="{newVersion}"')
			    xml = xml.replace(f'scenario_{currentVersion}', f'scenario_{newVersion}')

			    scenarioFile = 'scenario_'+xmlFile+'.xml'
			    with open(f'{workdir}{s}/{scenarioFile}', 'w') as fo:
			        fo.write(f'{xml}')

			run(f'cd {workdir} && ../../../build/openMalaria -s {s}/{scenarioFile} --ctsout {s}/ctsout_{xmlFile}.txt --output {s}/{xmlFile}.txt')

def check():
	ranges = np.loadtxt('ranges.csv', delimiter=',')
	df = pd.read_csv(f'check/1.csv')
	df = df['new weighted LF'].transpose()
	df.reset_index(drop=True, inplace=True)
	
	success = True
	for i in range(10):
		if df[i] < ranges[i][0] or df[i] > ranges[i][1]:
			print(f'LF {i} is out of range: {df[i]} [{ranges[i][0]}, {ranges[i][1]}]')
			success = False
	if success:
		print('Check successful')

# runScenarios()
run('Rscript check.R')
check()