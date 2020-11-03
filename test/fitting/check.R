source('likelihood.scicore.R')

params = list()
params$DensityBiasnonGarki = 0.177378571# 0.142722
params$DensityBiasGarki = 4.796107725 #0.142722

print(allPlots(1, 'check/', params))
